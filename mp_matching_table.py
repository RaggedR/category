"""Full table of MP matching coefficients: the MP analogue of the Hermite
matching matrix D_mk.

b_{mk} = <Phi_m, Phat_k>_W = int Phat_m(x) Phat_k(x) |Gamma(1/4+ix)|^3 dx

where Phat_n = P_n / sqrt(h_n) are orthonormal MP polynomials
(lambda = 1/4, phi = pi/2) and Phi_m = Phat_m * sqrt(W) are the MP functions.

Also computes the raw integrals J_{mk} = int P_m P_k W^{3/2} dx and probes
the ratios J_{mk}/J_{00} for algebraic structure (the Hermite analogue had
I_mk / I_00 * 3^{(m+k)/2} integer).
"""
import mpmath as mp
import pickle, os
from mp_matching import W, evalP, quadW, NODES

mp.mp.dps = 40
MMAX = 11
CACHE = "mp_matching_table.pkl"


def main():
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            nrm, J = pickle.load(f)
    else:
        nrm = [quadW(lambda x, n=n: evalP(n, x)**2 * W(x))
               for n in range(MMAX + 1)]
        J = {}
        for m in range(MMAX + 1):
            for k in range(m, MMAX + 1):
                if (m + k) % 2:
                    continue
                J[(m, k)] = quadW(lambda x, m=m, k=k:
                                  evalP(m, x) * evalP(k, x) * W(x)**mp.mpf(1.5))
        with open(CACHE, "wb") as f:
            pickle.dump((nrm, J), f)

    def b(m, k):
        mm, kk = min(m, k), max(m, k)
        return J[(mm, kk)] / mp.sqrt(nrm[m] * nrm[k])

    # ---- orthonormal table, two parity blocks ----
    for par, name in [(0, "EVEN"), (1, "ODD")]:
        idx = list(range(par, MMAX + 1, 2))
        print(f"\n=== b_mk ({name} block), lambda=1/4, phi=pi/2 ===")
        print("      " + "".join(f"k={k:<11d}" for k in idx))
        for m in idx:
            print(f"m={m:<3d} " + "".join(f"{float(b(m,k)): .5e} " for k in idx))

    # ---- raw ratios: hunt for structure ----
    print("\n=== raw ratios  J_mk / J_00  (Hermite analogue: rational * 3^{-(m+k)/2}) ===")
    J00 = J[(0, 0)]
    for par in [0, 1]:
        idx = list(range(par, MMAX + 1, 2))
        for m in idx:
            for k in idx:
                if k < m:
                    continue
                r = J[(m, k)] / J00
                print(f"  J_{m},{k}/J_00 = {mp.nstr(r, 20)}"
                      f"   identify: {mp.identify(r, ['sqrt(2)', 'sqrt(3)', 'sqrt(5)', 'pi'])}")

    # ---- decay along row 0 (predicted rho from balance law) ----
    print("\n=== decay along m=0 row ===")
    idx = list(range(0, MMAX + 1, 2))
    for i in range(1, len(idx)):
        k0, k1 = idx[i - 1], idx[i]
        sl = (mp.log(abs(b(0, k1))) - mp.log(abs(b(0, k0)))) / (k1 - k0)
        print(f"  k={k1:2d}  ln|b_0k| = {float(mp.log(abs(b(0,k1)))):8.3f}"
              f"   slope per k = {float(sl):8.4f}   rho_local = {float(mp.e**sl):.4f}")

    # ---- J00 itself: int |Gamma(1/4+ix)|^3 dx ----
    print(f"\nJ_00 = int |Gamma(1/4+ix)|^3 dx = {mp.nstr(J00, 25)}")
    print("  identify:", mp.identify(J00, ['pi', 'gamma(1/4)', 'gamma(3/4)', 'sqrt(pi)', 'sqrt(2)']))


if __name__ == "__main__":
    main()
