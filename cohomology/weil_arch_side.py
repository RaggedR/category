"""Archimedean side of the Weil form in the MP-1/4 basis, high precision.

Companion to weil_prime_side.py (prime side, 120 digits via zeta'/zeta series).
The old weil_form_mp.py arch entries were dps-30 oscillatory-node quadrature
(~1e-6); here the integrand is NON-oscillatory (Re psi smooth, weight decays
like e^{-pi r/2}), so plain mp.quad reaches ~dps digits.

Computes, for m,n <= 6, m = n (mod 2):
  A_mn  = (1/pi) int_0^inf h_mn(r) [Re psi(1/4+ir/2) - log pi] dr   (dps 130)
  pole_mn = 2 sqrt(pi) [PhPh(i/4) + PhPh(-i/4)]                      (dps 130)
  Z_mn  = 2 sum_j h_mn(gamma_j)   over zeros to height ~85           (dps 50)

Identity check (per validated convention): pole + arch - prime = Z.
Prime values are read from weil_prime_Pmn.pkl when it exists (Stage 1 of
weil_prime_side.py); the check runs at dps 50 (zero-side cost limits it).

Caches: weil_arch_Amn.pkl (130-digit strings), weil_zero_Zmn.pkl (50-digit).
"""
import os
import pickle

import mpmath as mp

MMAX = 6
A_CACHE = "weil_arch_Amn.pkl"
Z_CACHE = "weil_zero_Zmn.pkl"
P_CACHE = "weil_prime_Pmn.pkl"


def mp_poly(m, x):
    p0, p1 = mp.mpf(1), 2 * x
    if m == 0:
        return p0
    for j in range(1, m):
        p0, p1 = p1, (2 * x * p1 - (j - mp.mpf(1) / 2) * p0) / (j + 1)
    return p1


def h_norm(m):
    return 2 * mp.pi * mp.gamma(m + mp.mpf(1) / 2) / (mp.sqrt(2) * mp.factorial(m))


def h_mn(m, n, r):
    x = r / 2
    g = mp.gamma(mp.mpf(1) / 4 + 1j * x)
    return ((r**2 + mp.mpf(1) / 4) * mp_poly(m, x) * mp_poly(n, x)
            / mp.sqrt(h_norm(m) * h_norm(n)) * (g * mp.conj(g)).real)


def pairs():
    return [(m, n) for m in range(MMAX + 1) for n in range(m, MMAX + 1)
            if (m - n) % 2 == 0]


def arch_entry(m, n):
    kern = lambda r: h_mn(m, n, r) * (mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * r / 2))
                                      - mp.log(mp.pi))
    # non-oscillatory, weight ~ e^{-pi r/2}: 1e-140 by r ~ 205
    return (1 / mp.pi) * mp.quad(kern, [0, 2, 5, 10, 20, 40, 80, 140, 210])


def pole_entry(m, n):
    nrm = mp.sqrt(h_norm(m) * h_norm(n))
    z = mp_poly(m, mp.mpc(0, 0.25)) * mp_poly(n, mp.mpc(0, 0.25)) / nrm
    zb = mp_poly(m, mp.mpc(0, -0.25)) * mp_poly(n, mp.mpc(0, -0.25)) / nrm
    return (2 * mp.sqrt(mp.pi) * (z + zb)).real


def stage_arch(dps=130):
    print("=" * 68)
    print(f"ARCH + POLE at dps={dps}, m,n <= {MMAX}")
    print("=" * 68)
    mp.mp.dps = dps
    A, POLE = {}, {}
    for (m, n) in pairs():
        a = arch_entry(m, n)
        p = pole_entry(m, n)
        A[(m, n)] = A[(n, m)] = a
        POLE[(m, n)] = POLE[(n, m)] = p
        print(f"  ({m},{n}): arch = {mp.nstr(a, 25)}   pole = {mp.nstr(p, 25)}",
              flush=True)
    with open(A_CACHE, "wb") as f:
        pickle.dump({k: (mp.nstr(A[k], dps), mp.nstr(POLE[k], dps))
                     for k in A}, f)
    print(f"  cached -> {A_CACHE}")
    return A, POLE


def stage_zeros(dps=50, height=85):
    print("=" * 68)
    print(f"ZERO SIDE at dps={dps}, zeros to height ~{height}")
    print("=" * 68)
    mp.mp.dps = dps
    gammas = []
    j = 1
    while True:
        g = mp.im(mp.zetazero(j))
        if g > height:
            break
        gammas.append(g)
        j += 1
    print(f"  {len(gammas)} zeros, weight at top ~ e^(-pi*{height}/2) ~ "
          f"{mp.nstr(mp.exp(-mp.pi * height / 2), 3)}")
    Z = {}
    for (m, n) in pairs():
        Z[(m, n)] = 2 * mp.fsum(h_mn(m, n, g) for g in gammas)
        print(f"  Z[{m},{n}] = {mp.nstr(Z[(m, n)], 25)}", flush=True)
    with open(Z_CACHE, "wb") as f:
        pickle.dump({k: mp.nstr(v, dps) for k, v in Z.items()}, f)
    print(f"  cached -> {Z_CACHE}")
    return Z


def identity_check(A, POLE, Z):
    print("=" * 68)
    print("IDENTITY CHECK: pole + arch - prime  vs  zero sum")
    print("=" * 68)
    if not os.path.exists(P_CACHE):
        print(f"  {P_CACHE} not found yet -- rerun after weil_prime_side.py finishes")
        return
    mp.mp.dps = 60
    with open(P_CACHE, "rb") as f:
        P = {k: mp.mpf(v) for k, v in pickle.load(f).items()}
    for (m, n) in pairs():
        w = mp.mpf(mp.nstr(POLE[(m, n)], 55)) + mp.mpf(mp.nstr(A[(m, n)], 55)) \
            - P[(m, n)]
        z = Z[(m, n)]
        print(f"  ({m},{n}): W = {mp.nstr(w, 20):>25}   Z = {mp.nstr(z, 20):>25}"
              f"   diff = {mp.nstr(w - z, 3)}", flush=True)


if __name__ == "__main__":
    A, POLE = stage_arch()
    Z = stage_zeros()
    identity_check(A, POLE, Z)
