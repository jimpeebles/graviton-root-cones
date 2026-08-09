#!/usr/bin/env python3
"""Independent re-derivation of Table 1 of the draft (chamber classification)."""
from itertools import combinations, product
from independent_check import br, vsum, trees, cut_side, V, flow_vertex

def sgn(v):
    return (v > 0) - (v < 0)

def resonance_chambers(m, bound):
    props = [s for k in range(1, m) for s in combinations(range(m), k)]
    seen = {}
    for coords in product(range(-bound, bound + 1), repeat=m - 1):
        x = coords + (-sum(coords),)
        vals = [sum(x[i] for i in s) for s in props]
        if any(v == 0 for v in vals):
            continue
        seen.setdefault(tuple(v > 0 for v in vals), x)
    return seen

def realize(x):
    """Physical spinors with increasing integer slopes, Q=(S,0), S>0; None if impossible."""
    m = len(x)
    pref = [sum(x[:k]) for k in range(1, m)]
    if not any(p < 0 for p in pref):
        return None
    # choose gaps g_k >= 1; S = -sum pref_k g_k ; boost the most negative prefix's gap
    for big in range(1, 4000):
        g = [1] * (m - 1)
        k = min(range(m - 1), key=lambda i: pref[i])
        g[k] = big
        S = -sum(p * gk for p, gk in zip(pref, g))
        if S <= 0:
            continue
        r = [0]
        for gk in g:
            r.append(r[-1] + gk)
        q = tuple((xi * ri, -xi) for xi, ri in zip(x, r))
        Q = vsum(q)
        assert Q == (S, 0)
        proj = [br(qi, Q) for qi in q]
        if [sgn(p) for p in proj] != [sgn(v) * 1 for v in [S * xi for xi in x]]:
            continue
        # confirm same chamber
        props = [s for kk in range(1, m) for s in combinations(range(m), kk)]
        if all(sgn(sum(proj[i] for i in s)) == sgn(sum(x[i] for i in s)) for s in props):
            # genericity for the vertex: no zero pair bracket
            if all(br(q[i], q[j]) != 0 for i, j in combinations(range(m), 2)):
                return q
    return None

def feasible_family(q):
    m = len(q)
    Q = vsum(q)
    x = [br(qi, Q) for qi in q]
    fam = []
    for T in trees(m):
        ok = True
        for e in T:
            u, v = e
            A = cut_side(T, e, u, m)
            sA = sum(x[i] for i in A)
            duv = br(q[u], q[v])
            exits = (duv < 0)
            if (sA > 0) != exits:
                ok = False
                break
        if ok:
            fam.append(T)
    return frozenset(fam)

def orient(q, e):
    a, b = e
    return (a, b) if br(q[a], q[b]) < 0 else (b, a)

def in_arb_toward(q, T, root, m):
    """Is retarded-oriented tree T an in-arborescence toward root?"""
    # BFS from root; every tree edge must point from child side to parent side
    adj = {i: set() for i in range(m)}
    for a, b in T:
        adj[a].add(b)
        adj[b].add(a)
    parent = {root: None}
    order = [root]
    for c in order:
        for x in adj[c]:
            if x not in parent:
                parent[x] = c
                order.append(x)
    for node in order[1:]:
        if orient(q, tuple(sorted((node, parent[node])))) != (node, parent[node]):
            return False
    return True

def bareiss(mat):
    n = len(mat)
    if n == 0:
        return 1
    a = [row[:] for row in mat]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            sw = next((r for r in range(k + 1, n) if a[r][k]), None)
            if sw is None:
                return 0
            a[k], a[sw] = a[sw], a[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[-1][-1]

def classify(m, bound):
    chambers = resonance_chambers(m, bound)
    fams = {}
    stats = dict(physical=0, empty=0, arb_nonempty=0, arb_nontrivial=0, det_ok=0, det_bad=0, maxfam=0)
    for sig, x in chambers.items():
        q = realize(x)
        if q is None:
            continue
        stats["physical"] += 1
        fam = feasible_family(q)
        fams[fam] = fams.get(fam, 0) + 1
        stats["maxfam"] = max(stats["maxfam"], len(fam))
        if not fam:
            stats["empty"] += 1
            continue
        union = frozenset(e for T in fam for e in T)
        arb_roots = []
        for root in range(m):
            gen = frozenset(
                T for T in trees(m)
                if all(e in union for e in T) and in_arb_toward(q, T, root, m)
            )
            if gen == fam:
                arb_roots.append(root)
        if arb_roots:
            stats["arb_nonempty"] += 1
            if len(fam) > 1:
                stats["arb_nontrivial"] += 1
            expected = V(list(q))
            for root in arb_roots:
                L = [[0] * m for _ in range(m)]
                for e in union:
                    t, h = orient(q, e)
                    w = abs(br(q[t], q[h]))
                    L[t][t] += w
                    L[t][h] -= w
                minor = [[L[r][c] for c in range(m) if c != root] for r in range(m) if r != root]
                d = bareiss(minor)
                stats["det_ok" if d == expected else "det_bad"] += 1
    print(f"m={m}: chambers={len(chambers)}, physical={stats['physical']}, "
          f"distinct families={len(fams)}, V=0 chambers={stats['empty']}, "
          f"arborescent nonempty={stats['arb_nonempty']}, nontrivial={stats['arb_nontrivial']}, "
          f"max family={stats['maxfam']}, det checks ok/bad={stats['det_ok']}/{stats['det_bad']}")

if __name__ == "__main__":
    classify(4, 6)
    classify(5, 8)