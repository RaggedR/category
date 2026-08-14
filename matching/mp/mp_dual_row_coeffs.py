"""Extend rows m=0,4 of the dual matching coefficients to k=40 -- via the GF.

Instead of quadrature, use the master generating function

    sum_{m,k} J_mk w^m u^k = 2 pi^{3/2} e^{zeta/2} K_0(zeta/2)
                             / sqrt((1-iw)(1-iu)),
    zeta = phi(w) phi(u),   phi(t) = (1+it)/(1-it),

and extract Taylor coefficients by discrete Fourier sums on circles
|w| = |u| = 1/2 (singularity at |t| = 1; aliasing error ~ (1/2)^M).
Then c_mk = J_mk / sqrt(h_m h_k) with h_n = sqrt(2) pi^{3/2} C(2n,n)/4^n.

Cache: mp_dual_row_coeffs.pkl  ->  {(m, k): c_mk} for m in {0, 4}, k <= 40.
"""
import pickle, os
import mpmath as mp

mp.mp.dps = 60
KMAX = 40
CACHE = "mp_dual_row_coeffs.pkl"


def h(n):
    return mp.sqrt(2) * mp.pi**mp.mpf(1.5) * mp.binomial(2*n, n) / 4**n


def G(w, u):
    zeta = (1 + 1j*w) * (1 + 1j*u) / ((1 - 1j*w) * (1 - 1j*u))
    return (2 * mp.pi**mp.mpf(1.5) * mp.exp(zeta/2) * mp.besselk(0, zeta/2)
            / mp.sqrt((1 - 1j*w) * (1 - 1j*u)))


def compute():
    Mw, Mu = 64, 128            # aliasing: (1/2)^64 ~ 5e-20, (1/2)^128 ~ 3e-39
    r = mp.mpf(1) / 2
    ws = [r * mp.expjpi(2 * mp.mpf(a) / Mw) for a in range(Mw)]
    us = [r * mp.expjpi(2 * mp.mpf(b) / Mu) for b in range(Mu)]
    vals = [[G(ws[a], us[b]) for b in range(Mu)] for a in range(Mw)]

    C = {}
    for m in (0, 4):
        # inner Fourier sum over w for fixed u-node
        inner = [sum(vals[a][b] * mp.expjpi(-2 * mp.mpf(m * a) / Mw)
                     for a in range(Mw)) / (Mw * r**m) for b in range(Mu)]
        for k in range(KMAX + 1):
            J = sum(inner[b] * mp.expjpi(-2 * mp.mpf(k * b) / Mu)
                    for b in range(Mu)) / (Mu * r**k)
            C[(m, k)] = J / mp.sqrt(h(m) * h(k))
    return C


if __name__ == "__main__":
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            C = pickle.load(f)
    else:
        C = compute()
        with open(CACHE, "wb") as f:
            pickle.dump(C, f)

    # sanity: compare with the 12x12 quadrature cache
    with open("mp_dual_matching.pkl", "rb") as f:
        nrm, Cq = pickle.load(f)
    worst = max(abs(C[(m, k)] - Cq[(min(m, k), max(m, k))])
                for m in (0, 4) for k in range(12))
    print(f"max |GF - quadrature| over rows 0,4 (k<=11): {mp.nstr(worst, 3)}")

    print("\ndecay of |c_0k|:")
    for k in range(0, KMAX + 1, 4):
        print(f"  k={k:2d}  |c_0k| = {mp.nstr(abs(C[(0, k)]), 8)}")
