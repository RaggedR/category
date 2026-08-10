"""Wilson lattice rank test — the falsification run for the D-module frame.

Frame prediction: lattice rank of the matching coefficients = D-module rank
of the kernel.  MP kernels (Bessel/Tricomi) are rank 2 -> ZA + ZB (verified).
Hahn kernel sqrt(z)U(1,1,z) is rank 2 reducible -> QA + Q*sqrt(pi) (verified
2026-08-06, after raising pslq maxcoeff).  Wilson's kernel should be rank >= 3,
so its dual-split coefficients c_mk should span a Q-vector space of dim >= 3.

Method: recompute c_mk = (1/pi^2) int_0^inf p_m p_k w_plus W dx  for m,k <= 5
at dps=80 (old pkl stores only float64 -> useless for PSLQ), then do
incremental PSLQ rank detection on the set {Re c_mk, Im c_mk}:
maintain a basis; for each value, ask pslq([v] + basis) with maxcoeff=1e12;
if a relation with nonzero leading coefficient appears (and verifies to
1e-50), v is dependent; else it joins the basis.  Rank = final basis size.

LESSON (Hahn 2026-08-06): mpmath pslq default maxcoeff=1000 hides lattices
whose integer coordinates reach 1e8.  Always maxcoeff >= 1e10.
"""
import os
import mpmath as mp
import pickle
import wilson_matching as wm

mp.mp.dps = int(os.environ.get("WRT_DPS", 80))

NMAX = 5
CACHE = os.environ.get("WRT_CACHE", "/Users/robin/git/category/wilson_rank_test.pkl")
MAXCOEFF = int(float(os.environ.get("WRT_MAXCOEFF", 1e12)))
MAXSTEPS = int(os.environ.get("WRT_MAXSTEPS", 20000))


def compute_C():
    H0 = mp.pi**2
    C = {}
    for m in range(NMAX + 1):
        for k in range(m, NMAX + 1):
            val = mp.quad(
                lambda x, m=m, k=k: wm.eval_wilson(m, x) * wm.eval_wilson(k, x)
                * wm.w_plus(x) * wm.W(x) / H0,
                wm.NODES)
            C[(m, k)] = val
            print(f"c_{m},{k} = {mp.nstr(val, 30)}", flush=True)
    return C


def verify_relation(rel, vec):
    resid = mp.fsum(r * v for r, v in zip(rel, vec))
    scale = max(abs(v) for v in vec) * max(abs(r) for r in rel)
    return abs(resid) < scale * mp.mpf(10) ** (-(mp.mp.dps - 25))


def rank_test(values, label):
    """Incremental Q-linear rank of a list of (name, mpf) pairs."""
    print(f"\n=== Rank test: {label} ===")
    basis = []          # list of (name, value)
    coords = {}         # name -> relation over basis (for reporting)
    for name, v in values:
        if abs(v) < mp.mpf('1e-40'):
            continue
        if not basis:
            basis.append((name, v))
            print(f"  {name}: NEW basis element #1")
            continue
        vec = [v] + [b for _, b in basis]
        try:
            rel = mp.pslq(vec, maxcoeff=MAXCOEFF, maxsteps=MAXSTEPS)
        except Exception as e:
            rel = None
            print(f"  {name}: pslq error {e}")
        if rel is not None and rel[0] != 0 and verify_relation(rel, vec):
            coords[name] = (rel, [n for n, _ in basis])
            print(f"  {name}: dependent, {rel[0]}*v = "
                  f"-({', '.join(f'{c}*{n}' for c, n in zip(rel[1:], [n for n,_ in basis]))})")
        else:
            basis.append((name, v))
            print(f"  {name}: NEW basis element #{len(basis)}"
                  + ("" if rel is None else f"  (pslq gave leading-0 or unverified rel)"))
    print(f"\n  RANK({label}) = {len(basis)}   basis: {[n for n, _ in basis]}")
    return basis, coords


def known_constant_scan(basis):
    """See whether basis elements relate to famous constants."""
    G14 = mp.gamma(mp.mpf(1)/4)
    G34 = mp.gamma(mp.mpf(3)/4)
    cands = [
        ('pi', mp.pi), ('pi^2', mp.pi**2), ('sqrt(pi)', mp.sqrt(mp.pi)),
        ('G14^2', G14**2), ('G34^2', G34**2), ('G14^4/pi', G14**4/mp.pi),
        ('e*E1(1)', mp.e * mp.e1(1)), ('e^2*E1(2)', mp.e**2 * mp.e1(2)),
        ('sqrt(pi)e^2E1(2)', mp.sqrt(mp.pi)*mp.e**2*mp.e1(2)),
        ('K0(1/2)', mp.besselk(0, mp.mpf(1)/2)),
        ('K1(1/2)', mp.besselk(1, mp.mpf(1)/2)),
    ]
    print("\n=== Known-constant scan for basis elements (pairs+triples) ===")
    for bn, bv in basis:
        found = False
        for i in range(len(cands)):
            for j in range(i + 1, len(cands)):
                vec = [bv, cands[i][1], cands[j][1]]
                rel = mp.pslq(vec, maxcoeff=MAXCOEFF, maxsteps=MAXSTEPS)
                if rel is not None and rel[0] != 0 and verify_relation(rel, vec):
                    print(f"  {bn}: {rel[0]}*v + {rel[1]}*{cands[i][0]} "
                          f"+ {rel[2]}*{cands[j][0]} = 0")
                    found = True
                    break
            if found:
                break
        if not found:
            print(f"  {bn} = {mp.nstr(bv, 25)}: no pair relation found")


def main():
    print(f"dps = {mp.mp.dps}, maxcoeff = {MAXCOEFF:.0e}")
    C = compute_C()
    with open(CACHE, 'wb') as f:
        pickle.dump({str(k): (mp.nstr(v.real, 70), mp.nstr(v.imag, 70))
                     for k, v in C.items()}, f)
    print(f"\nSaved high-precision strings to {CACHE}")

    # Order values by descending magnitude so basis elements are well-scaled
    reals, imags = [], []
    for (m, k), v in sorted(C.items()):
        reals.append((f"Re c{m}{k}", v.real))
        imags.append((f"Im c{m}{k}", v.imag))

    b_re, _ = rank_test(reals, "Re parts")
    b_im, _ = rank_test(imags, "Im parts")
    b_all, _ = rank_test(reals + imags, "all parts combined")

    known_constant_scan(b_all)

    print("\nDone. Frame prediction was rank >= 3 (per Re/Im sector).")


if __name__ == "__main__":
    main()
