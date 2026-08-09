#!/usr/bin/env python3
"""Referee's independent verification, written from scratch (no imports from the bundle).

Checks:
  A. Own implementation of the published vertex V / Vbar (arXiv:2603.04330 eq for V)
     and the Berends-Giele recursion, with BOTH candidate top-level conventions
     (sum over partitions with A>=2 vs A>=3 blocks).
  B. Decay-chamber anchor: compare both conventions with the published closed form
     M_{1..5}|_{R_{5,4}} = ([21]+[31]+[41])([32]+[42])[43]   (eq 39)
     and n=4 analogue prod_{i=1}^{n-2} sum_{j>i} [ji]        (eq 26).
  C. S_n permutation invariance of both conventions at n=5.
  D. Draft's C* chamber formula (eq newM5) vs both conventions, random chamber points.
  E. Symbolic check: Table 2 signs alone determine every step function for the
     16 four-block Cayley trees and the triples; recover the claimed tree families.
  F. Resonance chamber counts m=4,5; physical (some negative prefix) counts;
     all-positive-prefix counts vs R/m.
  G. Flow-cone equivalence (own implementation of Theorem 3.1) vs step functions.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
import random

# ---------- basic linear algebra on 2-spinors ----------
def br(a, b):
    return a[0] * b[1] - a[1] * b[0]

def vsum(vs):
    return (sum(v[0] for v in vs), sum(v[1] for v in vs))

# ---------- all labelled trees on n nodes (recursive edge growth) ----------
def all_trees(n):
    if n == 1:
        return [frozenset()]
    trees = []
    nodes = list(range(n))
    # grow trees by adding node k connected to any node j<k is wrong (only gives increasing trees);
    # use Cayley via parent arrays checked for connectivity instead: enumerate all functions
    # parent: {1..n-1} -> {0..n-1} that yield a tree on labels {0..n-1} rooted anywhere.
    # Simplest correct: enumerate all edge subsets of size n-1 and keep spanning trees.
    edges = list(combinations(nodes, 2))
    for subset in combinations(edges, n - 1):
        # connectivity check
        adj = {i: set() for i in nodes}
        for a, b in subset:
            adj[a].add(b); adj[b].add(a)
        seen = {0}; stack = [0]
        while stack:
            c = stack.pop()
            for x in adj[c]:
                if x not in seen:
                    seen.add(x); stack.append(x)
        if len(seen) == n:
            trees.append(frozenset(subset))
    return trees

TREE_CACHE = {}
def trees(n):
    if n not in TREE_CACHE:
        TREE_CACHE[n] = all_trees(n)
    return TREE_CACHE[n]

def cut_side(tree, edge, start, n):
    adj = {i: set() for i in range(n)}
    for a, b in tree:
        if (a, b) == edge:
            continue
        adj[a].add(b); adj[b].add(a)
    seen = {start}; stack = [start]
    while stack:
        c = stack.pop()
        for x in adj[c]:
            if x not in seen:
                seen.add(x); stack.append(x)
    return seen

# ---------- published vertex ----------
def V(qs, advanced=False):
    n = len(qs)
    if n == 1:
        return 1
    total = 0
    for T in trees(n):
        w = 1
        ok = True
        for e in T:
            u, v = e
            A = cut_side(T, e, u, n)      # u in A
            qA = vsum([qs[i] for i in A])
            qB = vsum([qs[i] for i in range(n) if i not in A])
            num = br(qA, qB)
            den = br(qs[u], qs[v])
            if num == 0 or den == 0:
                raise ZeroDivisionError("wall")
            arg = -Fraction(num, den)
            if advanced:
                arg = -arg
            if arg < 0:
                ok = False
                break
            w *= abs(den)
        if ok:
            total += w
    return total

# ---------- set partitions ----------
def partitions(items):
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p

# ---------- recursion ----------
def Mbar(qs, labels, memo):
    key = tuple(sorted(labels))
    if key in memo:
        return memo[key]
    if len(key) == 1:
        return 1
    if len(key) == 2:
        return 0
    total = 0
    for p in partitions(key):
        if len(p) < 3:
            continue
        blocks = [vsum([qs[i] for i in blk]) for blk in p]
        f = V(tuple(blocks))
        if f:
            for blk in p:
                f *= Mbar(qs, blk, memo)
        total -= f
    memo[key] = total
    return total

def M(qs_positive, min_blocks):
    """Stripped amplitude; minus leg = -sum(qs_positive). min_blocks = 2 or 3."""
    memo = {}
    labels = list(range(len(qs_positive)))
    total = 0
    for p in partitions(labels):
        if len(p) < min_blocks:
            continue
        blocks = tuple(vsum([qs_positive[i] for i in blk]) for blk in p)
        kern = V(blocks) - V(blocks, advanced=True)
        if kern:
            for blk in p:
                kern *= Mbar(qs_positive, blk, memo)
        total -= kern
    return total

# ================= B. decay chamber anchor =================
def decay_sample(n, rng):
    """omega_a>0 a=1..n-1; ztilde_1<...<z_{n-2}<z_n<z_{n-1}; z_n = weighted mean."""
    while True:
        w = [rng.randint(1, 6) for _ in range(n - 2)] + [rng.randint(20, 40)]
        z = sorted(rng.sample(range(-30, 30), n - 1))
        # assign z_1..z_{n-2} the smallest, z_{n-1} the largest
        zs = z[:n - 2] + [z[-1]]
        q = [(wa, wa * za) for wa, za in zip(w, zs)]  # q_i = omega_i*(1, z_i) -> store (omega, omega*z)
        qn = tuple(-c for c in vsum(q))
        # qn = omega_n*(1, z_n): omega_n = qn[0] (<0), z_n = qn[1]/qn[0]
        if qn[0] == 0:
            continue
        zn = Fraction(qn[1], qn[0])
        if not (zs[n - 3] < zn < zs[-1]):
            continue
        allq = [tuple(x) for x in q] + [qn]
        if any(br(a, b) == 0 for a, b in combinations(allq, 2)):
            continue
        return [tuple(x) for x in q], allq

def decay_product(allq, n):
    # prod_{i=1}^{n-2} sum_{j=i+1}^{n-1} [ji]
    total = 1
    for i in range(1, n - 1):
        s = sum(br(allq[j - 1], allq[i - 1]) for j in range(i + 1, n))
        total *= s
    return total

def check_decay(n, samples=25, seed=11):
    rng = random.Random(seed)
    agree2 = agree3 = agree2neg = agree3neg = 0
    for _ in range(samples):
        qpos, allq = decay_sample(n, rng)
        try:
            m2 = M(qpos, 2)
            m3 = M(qpos, 3)
        except ZeroDivisionError:
            continue
        target = decay_product(allq, n)
        agree2 += (m2 == target); agree2neg += (m2 == -target)
        agree3 += (m3 == target); agree3neg += (m3 == -target)
    print(f"[B] n={n} decay-chamber: A>=2 matches eq(26): +{agree2}/-{agree2neg}; "
          f"A>=3 matches: +{agree3}/-{agree3neg}  (of {samples})")

# ================= C. permutation invariance =================
def check_perm(n, samples=6, seed=7):
    rng = random.Random(seed)
    bad2 = bad3 = 0; tested = 0
    while tested < samples:
        qpos = [(rng.randint(-8, 8) or 3, rng.randint(-8, 8) or -2) for _ in range(n - 1)]
        qn = tuple(-c for c in vsum(qpos))
        allq = [tuple(x) for x in qpos] + [qn]
        if any(br(a, b) == 0 for a, b in combinations(allq, 2)):
            continue
        try:
            base2 = M(qpos, 2); base3 = M(qpos, 3)
            perms = list(permutations(range(n)))
            rng.shuffle(perms)
            for pi in perms[:20]:
                arr = [allq[i] for i in pi]
                if M(arr[:-1], 2) != base2:
                    bad2 += 1; break
            for pi in perms[:20]:
                arr = [allq[i] for i in pi]
                if M(arr[:-1], 3) != base3:
                    bad3 += 1; break
        except ZeroDivisionError:
            continue
        tested += 1
    print(f"[C] n={n} S_n invariance over {tested} pts: A>=2 failures={bad2}, A>=3 failures={bad3}")

# ================= D. draft C* chamber formula =================
CSTAR_ROOT = (-1, 1, 1, -1)          # signs of x_i
CSTAR_PAIRSUM = {(0,1): -1, (0,2): -1, (0,3): -1, (1,2): 1, (1,3): 1, (2,3): 1}
CSTAR_PAIRBR = {(0,1): -1, (0,2): -1, (0,3): 1, (1,2): 1, (1,3): -1, (2,3): -1}
CSTAR_TRIPLE = {(0,1,2): (-1,1,1), (0,1,3): (1,1,-1), (0,2,3): (-1,1,-1), (1,2,3): (1,-1,1)}

def sgn(v):
    return (v > 0) - (v < 0)

def in_cstar(q):
    Q = vsum(q)
    x = [br(qi, Q) for qi in q]
    if tuple(sgn(v) for v in x) != CSTAR_ROOT:
        return False
    for (i, j), s in CSTAR_PAIRSUM.items():
        if sgn(x[i] + x[j]) != s:
            return False
    for (i, j), s in CSTAR_PAIRBR.items():
        if sgn(br(q[i], q[j])) != s:
            return False
    for labs, ss in CSTAR_TRIPLE.items():
        qS = vsum([q[i] for i in labs])
        if tuple(sgn(br(q[i], qS)) for i in labs) != ss:
            return False
    return True

def draft_formula(q):
    Q = vsum(q)
    q5 = tuple(-c for c in Q)
    allq = list(q) + [q5]
    def w(i, j):
        return abs(br(allq[i - 1], allq[j - 1]))
    return (w(1,3)*w(1,4)*(w(1,2) + w(2,3))
            - w(1,4)*w(2,4)*(w(2,3) + w(3,4))
            - w(2,5)*w(1,4)*w(3,4)
            - w(3,5)*w(1,2)*w(2,4))

def check_cstar(samples=40, seed=99):
    rng = random.Random(seed)
    n2 = n3 = tot = 0
    while tot < samples:
        # random x in root chamber with some negative prefix, random increasing slopes
        a, b_, c = rng.randint(1,9), rng.randint(1,9), rng.randint(1,9)
        d = b_ + c - a
        if d <= 0:
            continue
        x = (-a, b_, c, -d)
        g = [rng.randint(1,9) for _ in range(3)]
        r = [0, g[0], g[0]+g[1], g[0]+g[1]+g[2]]
        q = tuple((xi*ri, -xi) for xi, ri in zip(x, r))
        if vsum(q)[0] <= 0 or not in_cstar(q):
            continue
        try:
            m2 = M(list(q), 2); m3 = M(list(q), 3)
        except ZeroDivisionError:
            continue
        f = draft_formula(q)
        n2 += (m2 == f); n3 += (m3 == f); tot += 1
    print(f"[D] C* chamber: draft eq(newM5) == M(A>=2): {n2}/{tot}; == M(A>=3): {n3}/{tot}")

# ================= E. symbolic determination in C* =================
def check_symbolic():
    """Resolve every step function for all 16 trees using ONLY the C* sign table."""
    xs = CSTAR_ROOT
    def subset_sign(A):
        # sign of sum_{i in A} x_i from table (singles/pairs/triples via complement)
        A = tuple(sorted(A))
        if len(A) == 1: return xs[A[0]]
        if len(A) == 2: return CSTAR_PAIRSUM[A]
        if len(A) == 3:
            comp = tuple(i for i in range(4) if i not in A)[0]
            return -xs[comp]
        raise ValueError
    surv_V, surv_Vb = [], []
    for T in trees(4):
        okV = okVb = True
        for e in T:
            u, v = e
            A = cut_side(T, e, u, 4)
            sA = subset_sign(A)
            sden = CSTAR_PAIRBR[e]
            argV = -sA * sden       # sign of -[A,B]/[uv]
            if argV < 0: okV = False
            if -argV < 0: okVb = False
        if okV: surv_V.append(sorted(tuple(sorted((a+1,b+1))) for a,b in T))
        if okVb: surv_Vb.append(sorted(tuple(sorted((a+1,b+1))) for a,b in T))
    print(f"[E] F(V_1234) from signs alone: {surv_V}")
    print(f"[E] F(Vbar_1234) from signs alone: {surv_Vb}")
    # triples
    for labs, ss in CSTAR_TRIPLE.items():
        names = {labs[i]: ss[i] for i in range(3)}
        surv = []
        for T in trees(3):
            ok = True
            for e in T:
                u, v = e
                lu, lv = labs[u], labs[v]
                A = cut_side(T, e, u, 3)
                # cut sign: [q_A, q_S - q_A] = [q_A, q_S]; for singleton A={u}: y_u
                if len(A) == 1:
                    sA = names[labs[next(iter(A))]] if False else names[lu]
                else:
                    other = next(i for i in range(3) if i not in A)
                    sA = -names[labs[other]]
                sden = CSTAR_PAIRBR[tuple(sorted((lu, lv)))]
                if -sA * sden < 0:
                    ok = False; break
            if ok:
                surv.append(sorted(tuple(sorted((labs[a]+1, labs[b]+1))) for a, b in T))
        print(f"[E] F(V_{{{''.join(str(l+1) for l in labs)}}}) = {surv}")

# ================= F. chamber counts =================
def chamber_counts(m):
    subsets = [s for k in range(1, m) for s in combinations(range(m), k) if 0 in s or True]
    # use all nonempty proper subsets; dedupe sign vectors via canonical half
    props = [s for k in range(1, m) for s in combinations(range(m), k)]
    props = [s for s in props if len(s) < m]
    seen = {}
    bound = 6 if m == 4 else 8
    for coords in product(range(-bound, bound + 1), repeat=m - 1):
        x = coords + (-sum(coords),)
        vals = [sum(x[i] for i in s) for s in props]
        if any(v == 0 for v in vals):
            continue
        sig = tuple(v > 0 for v in vals)
        seen.setdefault(sig, x)
    total = len(seen)
    phys = allpos = 0
    for sig, x in seen.items():
        pref = [sum(x[:k]) for k in range(1, m)]
        if any(p < 0 for p in pref):
            phys += 1
        if all(p > 0 for p in pref):
            allpos += 1
    print(f"[F] m={m}: chambers={total}, some-neg-prefix (physical)={phys}, "
          f"all-pos-prefix={allpos}, R/m={total}/{m}={total//m if total%m==0 else total/m}")

# ================= G. flow-cone equivalence =================
def flow_vertex(qs, advanced=False):
    n = len(qs)
    Q = vsum(qs)
    x = [br(qi, Q) for qi in qs]
    total = 0
    for T in trees(n):
        ok = True
        w = 1
        for e in T:
            u, v = e
            A = cut_side(T, e, u, n)
            sA = sum(x[i] for i in A)
            duv = br(qs[u], qs[v])
            if sA == 0 or duv == 0:
                raise ZeroDivisionError("wall")
            # retarded orientation u->v iff [uv]<0; edge exits A (u in A) iff u->v
            exits = (duv < 0)
            if advanced:
                exits = not exits
            if (sA > 0) != exits:
                ok = False; break
            w *= abs(duv)
        if ok:
            total += w
    return total

def check_flow(samples=60, seed=3):
    rng = random.Random(seed)
    for n in range(2, 7):
        done = 0
        while done < samples:
            qs = [(rng.randint(-9, 9) or 2, rng.randint(-9, 9) or -3) for _ in range(n)]
            try:
                a, b_ = V(qs), V(qs, advanced=True)
                c, d = flow_vertex(qs), flow_vertex(qs, advanced=True)
            except ZeroDivisionError:
                continue
            assert (a, b_) == (c, d), (n, qs, a, b_, c, d)
            done += 1
        print(f"[G] flow-cone == step-function vertex, n={n}: {done}/{samples}")

if __name__ == "__main__":
    check_decay(4)
    check_decay(5)
    check_perm(5)
    check_cstar()
    check_symbolic()
    chamber_counts(4)
    chamber_counts(5)
    check_flow()