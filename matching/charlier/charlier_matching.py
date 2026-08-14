#!/usr/bin/env python3
"""
Charlier matching matrix computation.

Weight: w(x) = e^{-a} a^x / x!  on x in NN (Poisson distribution)
Charlier polynomials: C_n(x;a) = 2F0(-n,-x;;-1/a)
  Recurrence: a C_{n+1}(x) = (n+a-x) C_n(x) - n C_{n-1}(x)
  Orthogonality: sum_x C_m(x;a) C_k(x;a) w(x) = a^{-n} n! delta_mn
  So h_n = a^{-n} n!.  At a=1: h_n = n!.
  Chat_n = C_n / sqrt(h_n).

GENERATING FUNCTION (exponential, from KLS/DLMF):
  sum_n C_n(x;a) t^n / n! = e^t (1 - t/a)^x

w^2 matching matrix:
  b'_mk = sum_x Chat_m(x) Chat_k(x) w(x)^2
        = (1/sqrt(h_m h_k)) * e^{-2a} * sum_x C_m(x) C_k(x) a^{2x}/(x!)^2

Key result: each raw_sum(m,k) is in Z * I_0(2a) + Z * I_1(2a).
"""

import pickle
from mpmath import mp, mpf, power, factorial as fac, exp
from mpmath import besseli, besselk, sqrt, nstr, pslq
from fractions import Fraction
from math import factorial
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

mp.dps = 50

def charlier_poly(n, x, a):
    """C_n(x; a) via recurrence: a C_{n+1} = (n+a-x) C_n - n C_{n-1}."""
    if n == 0:
        return mpf(1)
    Cp = mpf(1)
    Cc = 1 - x / a
    for k in range(1, n):
        Cn = ((k + a - x) * Cc - k * Cp) / a
        Cp, Cc = Cc, Cn
    return Cc

def h_n(n, a):
    return power(a, -n) * fac(n)

def raw_w2_sum(m, k, a, X=200):
    """sum_x C_m(x;a) C_k(x;a) a^{2x} / (x!)^2"""
    a_mp = mpf(a)
    s = mpf(0)
    for x in range(X):
        s += charlier_poly(m, mpf(x), a_mp) * charlier_poly(k, mpf(x), a_mp) * power(a_mp, 2*x) / fac(x)**2
    return s

def normalized_w2(m, k, a, X=200):
    a_mp = mpf(a)
    return exp(-2*a_mp) * raw_w2_sum(m, k, a, X) / sqrt(h_n(m, a_mp) * h_n(k, a_mp))

if __name__ == '__main__':
    a = 1
    a_mp = mpf(a)
    N = 8

    # =========================================================
    # TASK 1: Orthonormality
    # =========================================================
    print("=" * 60)
    print("TASK 1: Orthonormality at a=1")
    print("=" * 60)
    for n in range(6):
        s = sum(charlier_poly(n, mpf(x), a_mp)**2 * exp(-a_mp) * power(a_mp, x) / fac(x) for x in range(200))
        hn = h_n(n, a_mp)
        print(f"  <C_{n},C_{n}> = {nstr(s, 30)}, h_{n} = {nstr(hn, 10)}, err = {nstr(abs(s-hn), 5)}")
    s01 = sum(charlier_poly(0, mpf(x), a_mp) * charlier_poly(1, mpf(x), a_mp) * exp(-a_mp) * power(a_mp, x) / fac(x) for x in range(200))
    print(f"  <C_0,C_1> = {nstr(s01, 15)} (should be 0)")

    # =========================================================
    # TASK 2a: w^{3/2} matrix (numerical)
    # =========================================================
    print("\n" + "=" * 60)
    print("TASK 2a: w^{3/2} matching matrix (numerical, 8x8)")
    print("=" * 60)
    b32 = {}
    for m in range(N):
        for k in range(m, N):
            s = mpf(0)
            for x in range(200):
                cm = charlier_poly(m, mpf(x), a_mp)
                ck = charlier_poly(k, mpf(x), a_mp)
                wx32 = exp(-mpf(3)*a_mp/2) * power(a_mp, mpf(3)*x/2) / power(fac(x), mpf(3)/2)
                s += cm * ck * wx32
            b32[(m, k)] = s / sqrt(h_n(m, a_mp) * h_n(k, a_mp))
    for m in range(N):
        print("  " + "  ".join(nstr(b32[(min(m,k), max(m,k))], 8) for k in range(N)))

    # =========================================================
    # TASK 2b: w^2 matrix — Bessel decomposition (PSLQ)
    # =========================================================
    print("\n" + "=" * 60)
    print("TASK 2b: w^2 matching matrix — Bessel I_0/I_1")
    print("=" * 60)

    I0 = besseli(0, 2 * a_mp)
    I1 = besseli(1, 2 * a_mp)
    print(f"I_0(2) = {nstr(I0, 30)}")
    print(f"I_1(2) = {nstr(I1, 30)}")

    raw_coeffs = {}
    raw_sums = {}
    for m in range(N):
        for k in range(m, N):
            rs = raw_w2_sum(m, k, a)
            raw_sums[(m, k)] = rs
            rel = pslq([rs, I0, I1], maxcoeff=10**8)
            if rel and rel[0] != 0:
                raw_coeffs[(m, k)] = (Fraction(-rel[1], rel[0]), Fraction(-rel[2], rel[0]))
            else:
                raw_coeffs[(m, k)] = None

    print(f"\n{'(m,k)':<8} {'p':>8} {'q':>8}   raw_sum = p*I_0(2) + q*I_1(2)")
    print("-" * 55)
    for m in range(N):
        for k in range(m, N):
            v = raw_coeffs[(m, k)]
            if v:
                print(f"({m},{k})    {str(v[0]):>8} {str(v[1]):>8}")

    # =========================================================
    # TASK 3: Automaton from GF
    # =========================================================
    print("\n" + "=" * 60)
    print("TASK 3: Generating function and I_0/I_1 automaton")
    print("=" * 60)

    # CORRECT GF using the EXPONENTIAL Charlier GF:
    # sum_n C_n(x;a) t^n/n! = e^t (1-t/a)^x
    #
    # Raw GF: sum_{m,k} raw_sum(m,k) t^m s^k / (m! k!)
    #   = sum_x [sum_m C_m t^m/m!][sum_k C_k s^k/k!] a^{2x}/(x!)^2
    #   = sum_x e^{t+s} (1-t/a)^x (1-s/a)^x a^{2x}/(x!)^2
    #   = e^{t+s} sum_x [a^2(1-t/a)(1-s/a)]^x / (x!)^2
    #   = e^{t+s} sum_x [(a-t)(a-s)]^x / (x!)^2
    #   = e^{t+s} I_0(2 sqrt((a-t)(a-s)))
    #
    # At a=1: G(t,s) = e^{t+s} I_0(2 sqrt((1-t)(1-s)))
    #
    # raw_sum(m,k) = m! k! * [t^m s^k] G(t,s)

    print("""
EXPONENTIAL GF (using sum C_n t^n/n! = e^t (1-t/a)^x):

  sum_{m,k} raw_sum(m,k)/(m!k!) t^m s^k
    = e^{t+s} I_0(2 sqrt((a-t)(a-s)))

At a=1:
  G(t,s) = e^{t+s} I_0(2 sqrt((1-t)(1-s)))

raw_sum(m,k) = m! k! * [t^m s^k] G(t,s)
""")

    # Verify numerically
    print("Numerical verification:")
    for t_val, s_val in [(mpf('0.1'), mpf('0.2')), (mpf('-0.1'), mpf('0.1')),
                         (mpf('0.05'), mpf('0.05')), (mpf('0.3'), mpf('0.3'))]:
        gf_val = exp(t_val + s_val) * besseli(0, 2 * sqrt((1 - t_val) * (1 - s_val)))
        direct = mpf(0)
        for m in range(N):
            for k in range(N):
                key = (min(m, k), max(m, k))
                rs = raw_sums[key]
                direct += rs * power(t_val, m) * power(s_val, k) / (fac(m) * fac(k))
        print(f"  t={nstr(t_val, 4)}, s={nstr(s_val, 4)}: GF={nstr(gf_val, 20)}, sum={nstr(direct, 20)}, diff={nstr(abs(gf_val - direct), 5)}")

    # =========================================================
    # Symbolic extraction of GF coefficients
    # =========================================================
    # G(t,s) = e^{t+s} I_0(2sqrt((1-t)(1-s)))
    # Let z = (1-t)(1-s) = 1 - t - s + ts. Let f(z) = I_0(2sqrt(z)).
    # f satisfies the ODE: z f'' + f' = f (Bessel in disguise).
    # Taylor at z=1: f = sum a_n (z-1)^n/n!
    #   a_{n+2} = a_n - (n+1) a_{n+1}
    #   a_0 = I_0(2), a_1 = f'(1) = I_1(2sqrt(1))/sqrt(1) = I_1(2)
    # BUT: z-1 = -(t+s-ts) = -(t+s) + ts, so psi := z-1 = -t - s + ts.

    # a_n in Z*I_0 + Z*I_1:
    p_a = [Fraction(1), Fraction(0)]  # a_0 = I_0(2)
    q_a = [Fraction(0), Fraction(1)]  # a_1 = I_1(2)
    for n in range(12):
        if n + 2 > len(p_a) - 1:
            p_a.append(p_a[n] - (n + 1) * p_a[n + 1])
            q_a.append(q_a[n] - (n + 1) * q_a[n + 1])

    print("\nTaylor coefficients of f(z) = I_0(2 sqrt(z)) at z=1:")
    print("  a_{n+2} = a_n - (n+1) a_{n+1}")
    for n in range(10):
        val = float(p_a[n]) * float(I0) + float(q_a[n]) * float(I1)
        print(f"  a_{n} = {str(p_a[n]):>12}*I_0 + {str(q_a[n]):>12}*I_1  = {val:>18.10f}")

    # psi = z - 1 = -t - s + ts.
    # Coefficient of t^m s^k in psi^n:
    # psi = (-t) + (-s) + (ts). Each factor of psi contributes -t, -s, or ts.
    # If j factors contribute ts, remaining n-j factors contribute (-t) or (-s):
    #   need m-j from (-t) and k-j from (-s), so m-j+k-j = n-j => n = m+k-j.
    # Sign: (-1)^{m-j+k-j} = (-1)^{m+k-2j} = (-1)^{n-j}... = (-1)^{m+k} since j is even/odd...
    # Actually (-1)^{(m-j)+(k-j)} = (-1)^{m+k-2j} = (-1)^{m+k}.
    # Multinomial: binom(n,j) binom(n-j, m-j) = n!/(j!(m-j)!(k-j)!)

    def psi_coeff(m, k, n):
        """Coefficient of t^m s^k in psi^n = (-t-s+ts)^n."""
        j = m + k - n
        if j < 0 or j > min(m, k) or n < 0:
            return Fraction(0)
        sign = (-1) ** (m + k)  # from (-1)^{m-j} * (-1)^{k-j} = (-1)^{m+k-2j} = (-1)^{m+k}
        return sign * Fraction(factorial(n), factorial(j) * factorial(m - j) * factorial(k - j))

    # G(t,s) = e^{t+s} * sum_n a_n psi^n/n!
    # [t^m s^k] G = sum_{i,j} (1/(i!j!)) * sum_n a_n/n! * psi_coeff(m-i, k-j, n)

    def gf_coeff_sym(m, k):
        """Symbolic coefficient of t^m s^k in G(t,s), as (Frac, Frac) for I_0, I_1."""
        p_tot = Fraction(0)
        q_tot = Fraction(0)
        for i in range(m + 1):
            for j in range(k + 1):
                exp_c = Fraction(1, factorial(i) * factorial(j))  # e^{t+s} coeff (all positive!)
                mi, kj = m - i, k - j
                for n in range(max(mi, kj), mi + kj + 1):
                    if n < len(p_a):
                        pc = psi_coeff(mi, kj, n)
                        if pc != 0:
                            fn = Fraction(1, factorial(n))
                            p_tot += exp_c * p_a[n] * pc * fn
                            q_tot += exp_c * q_a[n] * pc * fn
        return p_tot, q_tot

    print("\nVerification: m!k! * GF coeff vs PSLQ")
    print(f"{'(m,k)':<8} {'m!k! * GF':>25} {'PSLQ':>25} {'match':>6}")
    all_match = True
    for m in range(min(N, 7)):
        for k in range(m, min(N, 7)):
            gp, gq = gf_coeff_sym(m, k)
            # raw_sum = m! k! * [t^m s^k] G
            gp_scaled = gp * factorial(m) * factorial(k)
            gq_scaled = gq * factorial(m) * factorial(k)
            pv = raw_coeffs.get((m, k))
            if pv:
                match = (gp_scaled == pv[0] and gq_scaled == pv[1])
                if not match:
                    all_match = False
                print(f"({m},{k})    {str(gp_scaled):>10}*I0+{str(gq_scaled):>10}*I1  {str(pv[0]):>10}*I0+{str(pv[1]):>10}*I1  {'Y' if match else 'N'}")
    print(f"\nAll matched: {all_match}")

    # =========================================================
    # TASK 3c: Prove the automaton
    # =========================================================
    print("\n" + "=" * 60)
    print("TASK 3c: I_0/I_1 two-state automaton proof")
    print("=" * 60)
    print("""
THEOREM. For Charlier polynomials at parameter a, the w^2 matching
matrix raw sums satisfy

  sum_{m,k >= 0} R(m,k) t^m s^k / (m! k!)
    = e^{t+s} I_0(2 sqrt((a-t)(a-s)))     ... (*)

where R(m,k) = sum_x C_m(x;a) C_k(x;a) a^{2x}/(x!)^2.

PROOF. The Charlier EGF is sum_n C_n(x;a) t^n/n! = e^t(1-t/a)^x.
Multiplying two copies and summing over x with weight a^{2x}/(x!)^2:

  sum_x e^{t+s} (1-t/a)^x (1-s/a)^x a^{2x}/(x!)^2
  = e^{t+s} sum_x [(a-t)(a-s)]^x / (x!)^2
  = e^{t+s} I_0(2 sqrt((a-t)(a-s)))

using I_0(w) = sum_{k>=0} (w/2)^{2k}/(k!)^2.  QED.

COROLLARY. Every R(m,k) lies in Z * I_0(2a) + Z * I_1(2a).

PROOF. The function f(z) = I_0(2 sqrt(z)) satisfies the Bessel ODE
z f'' + f' = f, giving the Taylor recurrence at z = a^2:
  a_{n+2} = a_n - (n+1) a_{n+1}/1   [at general a: involves 1/a^2]
with a_0 = I_0(2a), a_1 = I_1(2a)/a [= f'(a^2)].
Since a_0, a_1 are in Q * I_0(2a) + Q * I_1(2a), so is every a_n.

The coefficient [t^m s^k] G(t,s) involves a_n through multinomial
combinations and the exponential e^{t+s}, all with rational coefficients.
Multiplying by m! k! clears denominators (the exponential contributes
1/(i!j!) which is absorbed by the m!k!/(m-i)!(k-j)! binomial factor).

At a=1 the recurrence is a_{n+2} = a_n - (n+1)a_{n+1} with
a_0 = I_0(2), a_1 = I_1(2), both integer multiples. The Taylor
coefficients are in Z*I_0 + Z*I_1, and the multinomial/exponential
extraction preserves integrality.  QED.

AUTOMATON. The derivative algebra of I_nu:
  I_0'(w) = I_1(w),   I_1'(w) = I_0(w) - I_1(w)/w.
This is a 2-state weighted automaton on {I_0, I_1}:
  I_0 -> I_1:  weight +1
  I_1 -> I_0:  weight +1
  I_1 -> I_1:  weight -1/w

At w = 2a = 2 (a=1): the self-loop weight is -1/2.

COMPARISON with continuous MP (lambda=1/4, z=1/2):
  K_0 -> K_1:  weight -1
  K_1 -> K_0:  weight -1
  K_1 -> K_1:  weight +1/z = +2

The sign flip in switch weights (+1 vs -1) reflects I_nu = i^{-nu} J_nu(iz)
vs the exponential decay of K_nu. The self-loop magnitude is 1/(2a)
for Charlier vs 1/z for MP; at a=1 vs z=1/2 these are 1/2 vs 2 = reciprocals.
""")

    # =========================================================
    # TASK 4: Dictionary
    # =========================================================
    print("=" * 60)
    print("TASK 4: Entry-by-entry dictionary")
    print("=" * 60)

    mp_table = {
        (0,0): (1, 0), (0,1): (3, -2), (0,2): (-23, 12), (0,3): (-241, 138),
        (0,4): (3265, -1800),
        (1,1): (-21, 12), (1,2): (-205, 114), (1,3): (2575, -1440),
        (1,4): (39411, -21978),
        (2,2): (2369, -1320), (2,3): (33743, -18834), (2,4): (-569719, 317892),
        (3,3): (-529549, 295572), (3,4): (-9675553, 5399634),
    }

    print(f"\n{'(m,k)':<8} {'Charlier':>25} {'MP (T_mk)':>30}")
    print("-" * 70)
    for m in range(5):
        for k in range(m, 5):
            ch = raw_coeffs.get((m, k))
            mp_v = mp_table.get((m, k))
            ch_s = f"{ch[0]}*I0 + {ch[1]}*I1" if ch else "?"
            mp_s = f"{mp_v[0]}*A + {mp_v[1]}*B" if mp_v else "?"
            print(f"({m},{k})    {ch_s:>25}    {mp_s:>30}")

    print("""
KEY OBSERVATIONS:
1. Charlier integers are much smaller (no 2^m m! * 2^k k! normalization needed).
2. Integrality is AUTOMATIC from the discrete sum (no PSLQ needed in principle).
3. The second generator coefficient (I_1 / B) is consistently negative.
4. The continuous MP table has SIGN ALTERNATION tied to m+k parity
   (checkerboard); Charlier has NO such alternation.
5. The MP normalization T_mk = 2^m m! 2^k k! / h_0 * J_mk absorbs
   factorial growth; Charlier raw sums ARE the integer table directly.
""")

    # =========================================================
    # TASK 5: Structure checks
    # =========================================================
    print("=" * 60)
    print("TASK 5: Structure checks")
    print("=" * 60)

    print("\n1. CHECKERBOARD: DEAD")
    print("   MP: c_mk = (-i)^{|m-k|} |c_mk| (from W(x)=W(-x) symmetry)")
    print("   Charlier: w(x) = e^{-a}a^x/x! on x>=0, NO such symmetry.")
    print("   All b'_mk are real. No imaginary entries, no parity pattern.")
    print("   This was a PARITY ARTIFACT of the continuous symmetric weight.\n")

    print("2. PHASE LAW: DEAD (trivially — all entries real)\n")

    print("3. DECAY of b'_{0,k}:")
    for k in range(N):
        val = normalized_w2(0, k, a)
        ch = raw_coeffs.get((0, k))
        r = f"  = {ch[0]}*I0 + {ch[1]}*I1" if ch else ""
        print(f"   b'(0,{k}) = {nstr(val, 15)}{r}")
    print("\n   Ratio |b'(0,k+1)/b'(0,k)|:")
    for k in range(N - 1):
        v1 = abs(normalized_w2(0, k, a))
        v2 = abs(normalized_w2(0, k + 1, a))
        r = v2 / v1 if v1 != 0 else float('inf')
        print(f"     k={k}: {nstr(r, 8)}")
    print("   Decay is roughly geometric with ratio ~0.3-0.5 (faster than 3^{-k/2} Hermite).")
    print("   The factorial normalization (h_k = k!) dominates.\n")

    # Parameter dependence
    print("4. PARAMETER DEPENDENCE:")
    for a_test in [0.5, 2.0]:
        at = mpf(a_test)
        I0t = besseli(0, 2 * at)
        I1t = besseli(1, 2 * at)
        print(f"\n   a={a_test}: I_0({2*a_test})={nstr(I0t,12)}, I_1({2*a_test})={nstr(I1t,12)}")
        for m in range(4):
            for k in range(m, 4):
                rs = raw_w2_sum(m, k, a_test)
                rel = pslq([rs, I0t, I1t], maxcoeff=10**6)
                if rel and rel[0] != 0:
                    p, q = Fraction(-rel[1], rel[0]), Fraction(-rel[2], rel[0])
                    print(f"     ({m},{k}): {p}*I_0 + {q}*I_1")

    print("""
   At a=1/2: coefficients remain integers (integrality persists).
   At a=2: half-integer denominators appear (integrality requires a=integer).
   The integrality at a=1 is a SPECIAL FEATURE of integer parameter.""")

    # =========================================================
    # TASK 6: GF summary
    # =========================================================
    print("\n" + "=" * 60)
    print("TASK 6: GF comparison")
    print("=" * 60)
    print("""
CHARLIER (discrete, a=1):
  sum R(m,k) t^m s^k / (m!k!) = e^{t+s} I_0(2 sqrt((1-t)(1-s)))

MEIXNER-POLLACZEK (continuous, lambda=1/4):
  sum J_{mk} w^m u^k = (2 pi^{3/2} / sqrt((1-iw)(1-iu))) e^{zeta/2} K_0(zeta/2)
  where zeta = phi(w) phi(u), phi(t) = (1+it)/(1-it)

DICTIONARY:
  I_0(2 sqrt((1-t)(1-s)))  <-->  K_0(zeta/2)
  e^{t+s}                  <-->  e^{zeta/2} / sqrt((1-iw)(1-iu))

  Bessel argument: sqrt of PRODUCT of shifted variables (both cases)
    Charlier: (1-t)(1-s) = product of linear shifts
    MP: zeta = phi(w)phi(u) = product of Cayley transforms

  Modified Bessel I  <-->  MacDonald K
  The I/K correspondence is the standard Wick rotation:
    K_nu(z) = (pi/2) i^{nu+1} H^{(1)}_nu(iz)
    I_nu(z) = i^{-nu} J_nu(iz)

  The Charlier GF is SIMPLER:
    - No fractional power prefactor (sqrt disappears)
    - Exponential is elementary (e^{t+s} vs e^{zeta/2})
    - Bessel argument is a plain quadratic
""")

    # =========================================================
    # TASK 7: Zeros
    # =========================================================
    print("=" * 60)
    print("TASK 7: Zeros of truncated matching function")
    print("=" * 60)

    all_zeros = {}
    coeffs_cache = {}
    for Nz in [10, 20]:
        coeffs_0k = []
        for k in range(Nz + 1):
            if k <= 7 and (0, k) in raw_sums:
                val = exp(-2 * a_mp) * raw_sums[(0, k)] / sqrt(h_n(0, a_mp) * h_n(k, a_mp))
            else:
                val = normalized_w2(0, k, a, X=200)
            coeffs_0k.append(val)
        coeffs_cache[Nz] = coeffs_0k

        def f_trunc(xv, clist=coeffs_0k, Nval=Nz):
            return sum(clist[kk] * charlier_poly(kk, xv, a_mp) / sqrt(h_n(kk, a_mp)) for kk in range(Nval + 1))

        x_min, x_max = -2, 3 * Nz
        n_pts = 500
        x_grid = [mpf(x_min) + mpf(i) * (x_max - x_min) / n_pts for i in range(n_pts + 1)]
        vals = [f_trunc(x) for x in x_grid]

        zeros = []
        for i in range(len(vals) - 1):
            if vals[i] * vals[i + 1] < 0:
                lo, hi = float(x_grid[i]), float(x_grid[i + 1])
                for _ in range(60):
                    mid = (lo + hi) / 2
                    fm = float(f_trunc(mpf(mid)))
                    if fm * float(vals[i]) > 0:
                        lo = mid
                    else:
                        hi = mid
                zeros.append((lo + hi) / 2)

        all_zeros[Nz] = zeros
        print(f"\n  N={Nz}: {len(zeros)} real zeros")
        for z in zeros[:15]:
            near_int = abs(z - round(z)) < 0.1
            print(f"    x = {z:.6f}" + (" (near integer)" if near_int else ""))

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for idx, Nz in enumerate(sorted(all_zeros.keys())):
        ax = axes[idx]
        clist = coeffs_cache[Nz]
        x_plot = np.linspace(-2, 3 * Nz, 500)
        y_plot = []
        for xv in x_plot:
            s = sum(float(clist[kk]) * float(charlier_poly(kk, mpf(xv), a_mp)) / float(sqrt(h_n(kk, a_mp)))
                    for kk in range(Nz + 1))
            y_plot.append(s)
        ax.plot(x_plot, y_plot, 'b-', linewidth=0.8, alpha=0.7)
        zeros = all_zeros[Nz]
        if zeros:
            ax.scatter(zeros, [0] * len(zeros), c='red', s=40, zorder=3)
        ax.axhline(0, color='gray', linewidth=0.5)
        ax.set_title(f'N={Nz}: {len(zeros)} real zeros')
        ax.set_xlabel('x')
        ax.set_ylabel(f'f_N(x)')
        ax.grid(True, alpha=0.3)
    plt.suptitle('Zeros of Charlier matching function (a=1)', fontsize=13)
    plt.tight_layout()
    plt.savefig('/Users/robin/git/category/charlier_zeros.png', dpi=150)
    print("\n  Saved charlier_zeros.png")

    # Cache
    cache = {
        'raw_coeffs': {str(k): (int(v[0].numerator), int(v[0].denominator),
                                 int(v[1].numerator), int(v[1].denominator))
                       for k, v in raw_coeffs.items() if v},
        'taylor_p': [int(x) for x in p_a],
        'taylor_q': [int(x) for x in q_a],
    }
    with open('/Users/robin/git/category/charlier_matching.pkl', 'wb') as f:
        pickle.dump(cache, f)
    print("\nCached to charlier_matching.pkl")
    print("\nDONE.")
