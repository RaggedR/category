"""Dual pair of MP functions and their polynomial approximations.

The weight factors analytically:  W(x) = Gamma(1/4+ix) * Gamma(1/4-ix).
Distribute one whole factor to each side of a BILINEAR pairing:

    F_m(x) = Gamma(1/4+ix) Phat_m(x)      (poles above, at x = i(1/4+n))
    G_m(x) = Gamma(1/4-ix) Phat_m(x)      (poles below)
    <F_m, G_k> = int F_m G_k dx = int Phat_m Phat_k W dx = delta_mk.

Approximating F_m by MP polynomials in L2(W):

    c_mk = <F_m, Phat_k>_W = int Phat_m Phat_k Gamma(1/4+ix)^2 Gamma(1/4-ix) dx

(G_m gives the conjugates).  Three whole Gamma factors -- no branch points.

Predictions tested here:
  * x -> -x  =>  c_mk = (-1)^{m+k} conj(c_mk):
    m+k even -> real entries, m+k odd -> purely imaginary (no parity vanishing!)
  * closed form of the base integral via Kummer U and Bessel K:
    int Gamma(1/4+ix)^2 Gamma(1/4-ix) dx
      = 2 pi Gamma(1/2)^2 U(1/2, 1, 1) = 2 pi^{3/2} sqrt(e) K_0(1/2).
  * decay of |c_0k| vs the k^{-3/4} of the sqrt(W) matching matrix.
"""
import mpmath as mp
import pickle, os
from mp_matching import evalP

mp.mp.dps = 40
MMAX = 11
CACHE = "mp_dual_matching.pkl"
NODES = [-50, -40, -30, -20, -12, -8, -5, -3, -1.5, 0,
         1.5, 3, 5, 8, 12, 20, 30, 40, 50]


def K(x):
    """The three-Gamma kernel Gamma(1/4+ix)^2 Gamma(1/4-ix)."""
    return mp.gamma(mp.mpf(1)/4 + 1j*x)**2 * mp.gamma(mp.mpf(1)/4 - 1j*x)


def W(x):
    g = mp.gamma(mp.mpf(1)/4 + 1j*x)
    return (g * mp.conj(g)).real


def main():
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            nrm, C = pickle.load(f)
    else:
        nrm = [mp.quad(lambda x, n=n: evalP(n, x)**2 * W(x), NODES)
               for n in range(MMAX + 1)]
        C = {}
        for m in range(MMAX + 1):
            for k in range(m, MMAX + 1):
                J = mp.quad(lambda x, m=m, k=k:
                            evalP(m, x) * evalP(k, x) * K(x), NODES)
                C[(m, k)] = J / mp.sqrt(nrm[m] * nrm[k])
        with open(CACHE, "wb") as f:
            pickle.dump((nrm, C), f)

    def c(m, k):
        return C[(min(m, k), max(m, k))]

    # ---- prediction 1: real / imaginary checkerboard ----
    worst = max(abs(mp.im(c(m, k))) if (m + k) % 2 == 0 else abs(mp.re(c(m, k)))
                for m in range(MMAX + 1) for k in range(MMAX + 1))
    print(f"checkerboard violation (should be ~0): {mp.nstr(worst, 3)}")

    # ---- prediction 2: closed form of the base integral ----
    base = C[(0, 0)] * nrm[0]          # = int K(x) dx
    closed = 2 * mp.pi**mp.mpf(1.5) * mp.sqrt(mp.e) * mp.besselk(0, mp.mpf(1)/2)
    print(f"int G^2 Gbar dx = {mp.nstr(base, 20)}")
    print(f"2 pi^3/2 e^1/2 K0(1/2) = {mp.nstr(closed, 20)}")

    # ---- the table ----
    print("\n=== c_mk = <F_m, Phat_k>_W  (m+k even: real; m+k odd: imaginary) ===")
    print("     " + "".join(f"k={k:<11d}" for k in range(MMAX + 1)))
    for m in range(MMAX + 1):
        row = []
        for k in range(MMAX + 1):
            v = c(m, k)
            row.append(f"{float(mp.re(v)): .4f}      " if (m + k) % 2 == 0
                       else f"{float(mp.im(v)): .4f}i     ")
        print(f"m={m:<3d} " + "".join(row))

    # ---- decay of row 0 ----
    print("\n=== decay of |c_0k| (sqrt(W) matching had ~ k^{-3/4}) ===")
    for k in range(1, MMAX + 1):
        v0, v1 = abs(c(0, k - 1)), abs(c(0, k))
        print(f"  k={k:2d}  |c_0k| = {mp.nstr(v1, 6)}   ratio = {mp.nstr(v1/v0, 6)}")


if __name__ == "__main__":
    main()
