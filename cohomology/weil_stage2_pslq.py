"""Stage 2: PSLQ structure hunt on the Weil-form prime side P_mn (and arch
side A_mn when weil_arch_Amn.pkl exists).

Inputs: weil_prime_Pmn.pkl (120-digit strings), weil_arch_Amn.pkl (130-digit).
Discipline (memory: Hahn/Wilson lessons): maxcoeff >= 1e10; ~12 digits per
PSLQ input is the spurious threshold, so at dps 120 keep vectors <= 8 long;
report every hit as UNVERIFIED until recomputed at 2x dps.

Tests, in order of strength:
  1. P_mn vs tier-1 {A, B} and tier-2 {A, B, Alog2, Blog2, 1} lattices.
  2. Same for arch entries A_mn (the dps-40 attempt was spurious).
  3. Q-relations among the P_mn themselves (incremental rank, like Wilson).
  4. P_mn vs {zeta'/zeta(3), zeta'/zeta(5), zeta'/zeta(7)} -- do a few zp/z
     terms exhaust the values? (They shouldn't exactly, but coefficients may
     be structured.)
"""
import os
import pickle

import mpmath as mp

mp.mp.dps = 115          # slightly under the 120-digit data
MAXCOEFF = int(1e10)
MAXSTEPS = 100000


def load(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


def verify(rel, vec):
    resid = mp.fsum(r * v for r, v in zip(rel, vec))
    scale = max(abs(v) for v in vec) * max(abs(r) for r in rel)
    return abs(resid) < scale * mp.mpf(10) ** (-(mp.mp.dps - 20))


def pslq_report(name, v, cands):
    vec = [v] + [c for _, c in cands]
    rel = mp.pslq(vec, maxcoeff=MAXCOEFF, maxsteps=MAXSTEPS)
    if rel is not None and rel[0] != 0 and verify(rel, vec):
        terms = " + ".join(f"{c}*{nm}" for c, (nm, _) in zip(rel[1:], cands) if c)
        print(f"    {name}: HIT (unverified@2x): {rel[0]}*v + {terms} = 0")
        return True
    print(f"    {name}: no relation")
    return False


def rank_test(values, label):
    print(f"\n  Rank test ({label}):")
    basis = []
    for name, v in values:
        if abs(v) < mp.mpf("1e-40"):
            continue
        if not basis:
            basis.append((name, v))
            continue
        vec = [v] + [b for _, b in basis]
        if len(vec) > 9:          # stay under spurious threshold
            basis.append((name, v))
            continue
        rel = mp.pslq(vec, maxcoeff=MAXCOEFF, maxsteps=MAXSTEPS)
        if rel is not None and rel[0] != 0 and verify(rel, vec):
            names = [n for n, _ in basis]
            print(f"    {name}: dependent: {rel[0]}*v = -({', '.join(f'{c}*{n}' for c, n in zip(rel[1:], names) if c)})")
        else:
            basis.append((name, v))
    print(f"    RANK({label}) >= {len(basis)}  basis head: {[n for n, _ in basis[:6]]}")


def main():
    A_ = mp.sqrt(2 * mp.e) * mp.besselk(0, mp.mpf(1) / 2)
    B_ = mp.sqrt(2 * mp.e) * mp.besselk(1, mp.mpf(1) / 2)
    tier1 = [("A", A_), ("B", B_)]
    tier2 = tier1 + [("Alog2", A_ * mp.log(2)), ("Blog2", B_ * mp.log(2)),
                     ("1", mp.mpf(1))]
    zp = [(f"zpz{2 * k + 1}", -mp.zeta(2 * k + 1, 1, 1) / mp.zeta(2 * k + 1))
          for k in range(1, 4)]

    P = {eval(k) if isinstance(k, str) else k: mp.mpf(v)
         for k, v in load("weil_prime_Pmn.pkl").items()}
    keys = sorted(k for k in P if k[0] <= k[1])

    print("=" * 68)
    print("TEST 1+4: prime entries P_mn vs lattices")
    print("=" * 68)
    for k in keys:
        print(f"  P[{k[0]},{k[1]}] = {mp.nstr(P[k], 20)}")
        pslq_report("tier1 {A,B}", P[k], tier1)
        pslq_report("tier2 {A,B,log2,1}", P[k], tier2)
        pslq_report("zeta'/zeta(3,5,7)", P[k], zp)

    print("=" * 68)
    print("TEST 3: Q-rank of the P_mn themselves")
    print("=" * 68)
    rank_test([(f"P{m}{n}", P[(m, n)]) for m, n in keys], "P entries")

    arch = load("weil_arch_Amn.pkl")
    if arch is None:
        print("\nweil_arch_Amn.pkl not ready -- rerun for arch tests.")
        return
    Am = {eval(k) if isinstance(k, str) else k: mp.mpf(v[0])
          for k, v in arch.items()}
    print("=" * 68)
    print("TEST 2: arch entries A_mn vs lattices (redo of spurious dps-40)")
    print("=" * 68)
    for k in sorted(kk for kk in Am if kk[0] <= kk[1]):
        print(f"  A[{k[0]},{k[1]}] = {mp.nstr(Am[k], 20)}")
        pslq_report("tier1 {A,B}", Am[k], tier1)
        pslq_report("tier2 {A,B,log2,1}", Am[k], tier2)
    rank_test([(f"A{m}{n}", Am[(m, n)]) for m, n in sorted(Am) if m <= n],
              "arch entries")


if __name__ == "__main__":
    main()
