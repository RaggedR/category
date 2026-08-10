"""Joint lambda=1/4 x lambda=3/4 Meixner-Pollaczek cross-matching.

Two MP families at phi=pi/2:
  lambda=1/4: (n+1)P_{n+1} = 2xP_n - (n-1/2)P_{n-1}, weight W_1/4 = |Gamma(1/4+ix)|^2
  lambda=3/4: (n+1)Q_{n+1} = 2xQ_n - (n+1/2)Q_{n-1}, weight W_3/4 = |Gamma(3/4+ix)|^2

Key identity (reflection formula):
  W_1/4(x) * W_3/4(x) = 2pi^2 / cosh(2pi x)

Elementary cross pairing:
  E_mk = int P_hat_m(x) Q_hat_k(x) * 2pi^2 sech(2pi x) dx

Since P_hat_m Q_hat_k is a polynomial times normalization, E_mk reduces to
Q-linear combinations of sech moments, which are Euler numbers.
"""
import mpmath as mp
import pickle, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction
from math import factorial

mp.mp.dps = 50
MMAX = 11  # 12x12 table
CACHE = "joint_sech_matching.pkl"
NODES = [-50, -40, -30, -20, -12, -8, -5, -3, -1.5, 0,
         1.5, 3, 5, 8, 12, 20, 30, 40, 50]


# ================================================================
# Polynomial evaluation
# ================================================================

def evalP(n, x):
    """Evaluate lambda=1/4 MP polynomial P_n(x)."""
    p0, p1 = mp.mpf(1), 2*x
    if n == 0: return p0
    for j in range(1, n):
        p0, p1 = p1, (2*x*p1 - (j - mp.mpf(1)/2)*p0) / (j + 1)
    return p1


def evalQ(n, x):
    """Evaluate lambda=3/4 MP polynomial Q_n(x)."""
    q0, q1 = mp.mpf(1), 2*x
    if n == 0: return q0
    for j in range(1, n):
        q0, q1 = q1, (2*x*q1 - (j + mp.mpf(1)/2)*q0) / (j + 1)
    return q1


def W14(x):
    g = mp.gamma(mp.mpf(1)/4 + 1j*x)
    return (g * mp.conj(g)).real

def W34(x):
    g = mp.gamma(mp.mpf(3)/4 + 1j*x)
    return (g * mp.conj(g)).real

def sech_kernel(x):
    return 2 * mp.pi**2 * mp.sech(2 * mp.pi * x)

def sqrtsech_kernel(x):
    return mp.pi * mp.sqrt(2) * mp.sqrt(mp.sech(2 * mp.pi * x))


# ================================================================
# Sech moments (exact via Euler numbers)
# ================================================================

def sech_moment_rational(twon):
    """mu_{2n} / pi^2 = |E_{2n}| / 16^n as a Fraction.
    Here twon is the exponent of x (must be even)."""
    n = twon // 2
    En = abs(int(mp.eulernum(twon)))
    return Fraction(En, 16**n)


# ================================================================
# Exact polynomial coefficients
# ================================================================

def mp_poly_coeffs(nmax, lam_frac):
    """Monomial coefficients of MP polys for given lambda (as Fraction)."""
    polys = [[Fraction(1)], [Fraction(0), Fraction(2)]]
    for n in range(1, nmax):
        pn, pm = polys[n], polys[n-1]
        nxt = [Fraction(0)] * (n + 2)
        for j, c in enumerate(pn):
            nxt[j+1] += 2 * c
        for j, c in enumerate(pm):
            nxt[j] -= (Fraction(n) + 2*lam_frac - 1) * c
        polys.append([c / (n + 1) for c in nxt])
    return polys


def norm_ratio(n, lam_frac):
    """h_n/h_0 for MP with parameter lambda (as Fraction)."""
    r = Fraction(1)
    for j in range(1, n+1):
        r *= (Fraction(j) + 2*lam_frac - 1) / Fraction(j)
    return r


# ================================================================
# Main computation
# ================================================================

def main():
    lam14 = Fraction(1, 4)
    lam34 = Fraction(3, 4)

    if os.path.exists(CACHE):
        print("Loading from cache...")
        with open(CACHE, "rb") as f:
            data = pickle.load(f)
        nrm14 = data["nrm14"]
        nrm34 = data["nrm34"]
        E_quad = data["E"]
        M_quad = data["M"]
    else:
        print("Computing norms for lambda=1/4...")
        nrm14 = [mp.quad(lambda x, n=n: evalP(n, x)**2 * W14(x), NODES)
                  for n in range(MMAX + 1)]
        print("Computing norms for lambda=3/4...")
        nrm34 = [mp.quad(lambda x, n=n: evalQ(n, x)**2 * W34(x), NODES)
                  for n in range(MMAX + 1)]

        # Verify reflection identity
        print("\nReflection identity W_1/4 * W_3/4 = 2pi^2 sech(2pi x):")
        for xv in [mp.mpf(0), mp.mpf(1), mp.mpf(2)]:
            lhs = W14(xv) * W34(xv)
            rhs = sech_kernel(xv)
            print(f"  x={float(xv)}: diff = {mp.nstr(abs(lhs-rhs), 4)}")

        print("\nComputing E_mk (sech kernel)...")
        E_quad = {}
        for m in range(MMAX + 1):
            for k in range(MMAX + 1):
                val = mp.quad(lambda x, m=m, k=k:
                              evalP(m, x) * evalQ(k, x) * sech_kernel(x), NODES)
                E_quad[(m, k)] = val / mp.sqrt(nrm14[m] * nrm34[k])
            print(f"  row {m} done", flush=True)

        print("Computing M_mk (sqrt-sech kernel)...")
        M_quad = {}
        for m in range(MMAX + 1):
            for k in range(MMAX + 1):
                val = mp.quad(lambda x, m=m, k=k:
                              evalP(m, x) * evalQ(k, x) * sqrtsech_kernel(x), NODES)
                M_quad[(m, k)] = val / mp.sqrt(nrm14[m] * nrm34[k])
            print(f"  row {m} done", flush=True)

        with open(CACHE, "wb") as f:
            pickle.dump({"nrm14": nrm14, "nrm34": nrm34, "E": E_quad, "M": M_quad}, f)

    # ================================================================
    # Exact evaluation
    # ================================================================
    print("\n" + "="*70)
    print("EXACT EVALUATION VIA EULER-NUMBER SECH MOMENTS")
    print("="*70)

    # Verify sech moments
    print("\nSech moments mu_{2n}/pi^2 = |E_{2n}|/16^n:")
    for n in range(8):
        r = sech_moment_rational(2*n)
        exact = r * mp.pi**2
        numer = mp.quad(lambda x, n2=2*n: x**(n2) * sech_kernel(x), NODES)
        print(f"  n={n}: {r} = {float(r):.10f}, quad diff = {mp.nstr(abs(exact-numer), 4)}")

    # Polynomial coefficients
    P_coeffs = mp_poly_coeffs(MMAX + 1, lam14)
    Q_coeffs = mp_poly_coeffs(MMAX + 1, lam34)

    # Compute exact r_mk = E_mk_unnorm / pi^2 as Fraction
    r_exact = {}
    for m in range(MMAX + 1):
        for k in range(MMAX + 1):
            deg = m + k
            prod = [Fraction(0)] * (deg + 1)
            for i, a in enumerate(P_coeffs[m]):
                for j, b in enumerate(Q_coeffs[k]):
                    prod[i + j] += a * b
            val = Fraction(0)
            for j in range(0, deg + 1, 2):
                val += prod[j] * sech_moment_rational(j)
            r_exact[(m, k)] = val

    # Parity: m+k odd => r_mk = 0 (only even powers survive)
    print("\nParity vanishing (m+k odd entries):")
    all_zero = True
    for m in range(8):
        for k in range(8):
            if (m + k) % 2 == 1 and r_exact[(m,k)] != 0:
                print(f"  NONZERO: r({m},{k}) = {r_exact[(m,k)]}")
                all_zero = False
    if all_zero:
        print("  All m+k odd entries are exactly ZERO. The matrix is block-diagonal by parity.")

    # Verify exact vs quadrature
    print("\nExact vs quadrature verification:")
    max_diff = mp.mpf(0)
    for m in range(MMAX + 1):
        for k in range(MMAX + 1):
            exact_val = float(r_exact[(m,k)]) * mp.pi**2 / mp.sqrt(nrm14[m] * nrm34[k])
            quad_val = mp.re(E_quad[(m,k)])
            diff = abs(exact_val - quad_val)
            if diff > max_diff:
                max_diff = diff
    print(f"  Max |exact - quadrature| = {mp.nstr(max_diff, 4)}")

    # ================================================================
    # Rational table
    # ================================================================
    print("\n=== Exact r_mk = E_mk_unnorm / pi^2 (rational) ===")
    for m in range(8):
        row = []
        for k in range(8):
            r = r_exact[(m, k)]
            if r == 0:
                row.append("0")
            else:
                row.append(f"{r.numerator}/{r.denominator}" if r.denominator != 1 else str(r.numerator))
        print(f"  m={m}: " + ", ".join(f"{s:>10s}" for s in row))

    # ================================================================
    # Integrality search
    # ================================================================
    print("\n" + "="*70)
    print("INTEGRALITY")
    print("="*70)

    # Look at denominators: they are powers of 2.
    # For the lambda=1/4 single family, the normalisation was 2^m m! 2^k k! / h_0.
    # Here the unnormalized r_mk already absorbs the polynomial coefficients
    # and the Euler-number denominators 16^n.
    # P_m coefficients have denominators that are products of integers up to m+1
    # (from the recurrence division by (n+1)), so roughly m!.
    # Similarly Q_k coefficients have denominators ~ k!.

    # Try T_mk = 4^m m! * 4^k k! * r_mk (the 4 from both the MP coeff denominators
    # and the 16^n in the sech moments)
    print("\nTrying T_mk = 4^m m! * 4^k k! * r_mk:")
    T1 = {}
    for m in range(MMAX + 1):
        for k in range(MMAX + 1):
            T1[(m,k)] = Fraction(4**m * factorial(m) * 4**k * factorial(k)) * r_exact[(m,k)]

    all_int = True
    for m in range(8):
        row = []
        for k in range(8):
            T = T1[(m,k)]
            if T.denominator != 1:
                all_int = False
                row.append(str(T))
            else:
                row.append(str(T.numerator))
        print(f"  m={m}: " + ", ".join(f"{s:>12s}" for s in row))

    if not all_int:
        # Try other normalizations
        print("\n4^m m! 4^k k! didn't give all integers. Trying 2^m m! * 2^k k!:")
        T2 = {}
        for m in range(MMAX + 1):
            for k in range(MMAX + 1):
                T2[(m,k)] = Fraction(2**m * factorial(m) * 2**k * factorial(k)) * r_exact[(m,k)]
        for m in range(8):
            row = []
            for k in range(8):
                T = T2[(m,k)]
                if T.denominator != 1:
                    row.append(str(T))
                else:
                    row.append(str(T.numerator))
            print(f"  m={m}: " + ", ".join(f"{s:>12s}" for s in row))

    # Find minimal normalization: for each (m,k), what's the denominator of r_mk?
    print("\nDenominators of r_mk (should be powers of 2):")
    for m in range(8):
        dens = []
        for k in range(8):
            if r_exact[(m,k)] != 0:
                dens.append(r_exact[(m,k)].denominator)
            else:
                dens.append(1)
        print(f"  m={m}: " + ", ".join(f"{d:>8d}" for d in dens))

    # Factor out powers of 2
    print("\nlog2 of denominators:")
    for m in range(8):
        logs = []
        for k in range(8):
            if r_exact[(m,k)] != 0:
                d = r_exact[(m,k)].denominator
                l = 0
                while d % 2 == 0:
                    d //= 2
                    l += 1
                if d != 1:
                    logs.append(f"{l}+{d}")
                else:
                    logs.append(str(l))
            else:
                logs.append("-")
        print(f"  m={m}: " + ", ".join(f"{s:>6s}" for s in logs))

    # From the output we saw:
    # (0,0)->0, (0,2)->3, (0,4)->7, (0,6)->11 => 4n-1 for (0,2n)
    # (1,1)->2, (1,3)->5, (1,5)->9 => 4n-2? n=1->2, n=2->5=4*2-3, n=3->9=4*3-3
    # Actually (m,k) with m+k=s: den ~ 2^{2s-1} for (0,s), 2^{2s-2} for (1,s-1)?
    # Let me just find the right normalization by trial.

    # Try: N_mk = (2m)!! * (2k)!! * r_mk   where (2m)!! = 2^m m!
    # That's 2^m m! * 2^k k! = same as T2. Already tried.

    # Try: (2m+1)!! * (2k+1)!! * r_mk? Or mixed: 2^m m! * (2k+1)!! ?
    # Given the two different lambda values, maybe the normalisation differs.
    # For lambda=1/4, the h_n ratio involves (2j-1)/(2j), denom = 2^n after n steps
    # For lambda=3/4, the h_n ratio involves (2j+1)/(2j), denom = 2^n after n steps
    # So sqrt(h_m * h_k) has denom 2^{(m+k)/2} ish.

    # Actually, let me just compute what power of 2 makes r_mk integer:
    print("\nMinimal power of 2 to integralize r_mk:")
    min_pow = {}
    for m in range(MMAX + 1):
        for k in range(MMAX + 1):
            if r_exact[(m,k)] == 0:
                min_pow[(m,k)] = 0
                continue
            r = r_exact[(m,k)]
            d = r.denominator
            p = 0
            while d % 2 == 0:
                d //= 2
                p += 1
            assert d == 1, f"Non-power-of-2 denominator at ({m},{k}): {r.denominator}"
            min_pow[(m,k)] = p

    for m in range(10):
        row = []
        for k in range(10):
            if (m+k) % 2 == 1:
                row.append("-")
            else:
                row.append(str(min_pow[(m,k)]))
        print(f"  m={m}: " + ", ".join(f"{s:>3s}" for s in row))

    # Define T_mk with the minimal integralization
    print("\n=== Integer table T_mk = 2^{p(m,k)} * r_mk ===")
    T_int = {}
    for m in range(MMAX + 1):
        for k in range(MMAX + 1):
            p = min_pow[(m,k)]
            T_int[(m,k)] = r_exact[(m,k)] * 2**p

    for m in range(8):
        row = []
        for k in range(8):
            if (m+k) % 2 == 1:
                row.append("0")
            else:
                T = T_int[(m,k)]
                row.append(str(T.numerator))
        print(f"  m={m}: " + ", ".join(f"{s:>8s}" for s in row))

    # ================================================================
    # Better: uniform normalization. Let's try T_mk = 2^{2(m+k)} * r_mk
    # ================================================================
    print("\n=== T_mk = 4^{m+k} * r_mk ===")
    for m in range(8):
        row = []
        for k in range(8):
            T = Fraction(4**(m+k)) * r_exact[(m,k)]
            if T.denominator == 1:
                row.append(str(T.numerator))
            else:
                row.append(str(T))
        print(f"  m={m}: " + ", ".join(f"{s:>12s}" for s in row))

    # Try T_mk = 2^{2(m+k)-1} * r_mk for m+k >= 1
    print("\n=== T_mk = 2 * 4^{m+k} * r_mk ===")
    for m in range(8):
        row = []
        for k in range(8):
            T = Fraction(2 * 4**(m+k)) * r_exact[(m,k)]
            if T.denominator == 1:
                row.append(str(T.numerator))
            else:
                row.append(str(T))
        print(f"  m={m}: " + ", ".join(f"{s:>12s}" for s in row))

    # ================================================================
    # Generating function
    # ================================================================
    print("\n" + "="*70)
    print("GENERATING FUNCTION")
    print("="*70)
    # sum_{m,k} I_mk w^m u^k where I_mk = int P_m Q_k * 2pi^2 sech(2pi x) dx
    # = int (1-iw)^{-1/4+ix}(1+iw)^{-1/4-ix} * (1-iu)^{-3/4+ix}(1+iu)^{-3/4-ix} * 2pi^2 sech(2pi x) dx
    # = (1+w^2)^{-1/4}(1+u^2)^{-3/4} * int e^{-isx} 2pi^2 sech(2pi x) dx
    #   where e^{is} = phi(w)phi(u), phi(t) = (1+it)/(1-it)
    # int e^{-isx} 2pi^2 sech(2pi x) dx = pi^2 sech(s/4)  [Fourier transform of sech]
    #
    # GF(w,u) = pi^2 (1+w^2)^{-1/4}(1+u^2)^{-3/4} sech(ln(phi(w)phi(u))/4)

    def GF(w, u):
        phi_w = (1 + 1j*w) / (1 - 1j*w)
        phi_u = (1 + 1j*u) / (1 - 1j*u)
        zeta = phi_w * phi_u
        s = mp.log(zeta)
        return mp.pi**2 * (1 + w**2)**(-mp.mpf(1)/4) * (1 + u**2)**(-mp.mpf(3)/4) * mp.sech(s/4)

    # Verify at w=u=0: should give I_{00} = pi^2 * r_{00} = pi^2
    print(f"GF(0,0) = {mp.nstr(GF(0,0), 15)}, pi^2 = {mp.nstr(mp.pi**2, 15)}")

    # Verify a few Taylor coefficients by finite differences
    h = mp.mpf(10)**(-8)
    checks = [
        ((0,0), GF(0,0)),
        ((1,0), (GF(h,0) - GF(-h,0))/(2*h)),
        ((0,1), (GF(0,h) - GF(0,-h))/(2*h)),
        ((2,0), (GF(h,0) - 2*GF(0,0) + GF(-h,0))/h**2 / 2),
        ((0,2), (GF(0,h) - 2*GF(0,0) + GF(0,-h))/h**2 / 2),
        ((1,1), (GF(h,h) - GF(h,-h) - GF(-h,h) + GF(-h,-h))/(4*h**2)),
    ]
    print("\nGF Taylor coefficient verification:")
    for (m,k), gf_val in checks:
        exact_val = r_exact[(m,k)] * mp.pi**2
        print(f"  I({m},{k}): GF={mp.nstr(gf_val, 12)}, exact={mp.nstr(exact_val, 12)}, "
              f"diff={mp.nstr(abs(gf_val - exact_val), 4)}")

    # Simplification of the GF
    # sech(ln(zeta)/4) = 2/(zeta^{1/4} + zeta^{-1/4})
    # With zeta = phi(w)phi(u):
    # GF = 2pi^2 / ((1+w^2)^{1/4}(1+u^2)^{3/4}(zeta^{1/4}+zeta^{-1/4}))
    #
    # zeta^{1/4} = [(1+iw)(1+iu)/((1-iw)(1-iu))]^{1/4}
    # Let alpha = [(1+iw)(1+iu)]^{1/4}, beta = [(1-iw)(1-iu)]^{1/4}
    # alpha*beta = [(1+w^2)(1+u^2)]^{1/4}
    # zeta^{1/4} = alpha/beta, zeta^{-1/4} = beta/alpha
    # zeta^{1/4}+zeta^{-1/4} = (alpha^2+beta^2)/(alpha*beta)
    #
    # GF = 2pi^2 * alpha*beta / ((1+w^2)^{1/4}(1+u^2)^{3/4}(alpha^2+beta^2))
    #    = 2pi^2 * [(1+w^2)(1+u^2)]^{1/4} / ((1+w^2)^{1/4}(1+u^2)^{3/4}(alpha^2+beta^2))
    #    = 2pi^2 / ((1+u^2)^{1/2}(alpha^2+beta^2))
    #
    # where alpha^2 = [(1+iw)(1+iu)]^{1/2}, beta^2 = [(1-iw)(1-iu)]^{1/2}
    #
    # GF(w,u) = 2pi^2 / (sqrt(1+u^2) * (sqrt((1+iw)(1+iu)) + sqrt((1-iw)(1-iu))))

    def GF2(w, u):
        A = mp.sqrt((1+1j*w)*(1+1j*u))
        B = mp.sqrt((1-1j*w)*(1-1j*u))
        return 2*mp.pi**2 / (mp.sqrt(1+u**2) * (A + B))

    print(f"\nGF2(0,0) = {mp.nstr(GF2(0,0), 15)} (should be pi^2 = {mp.nstr(mp.pi**2, 15)})")

    # Note the asymmetry: (1+u^2)^{1/2} not (1+w^2)^{1/2}!
    # This reflects lambda=1/4 vs lambda=3/4. The total exponent is
    # -1/4-3/4+1/4 = -3/4 for u and -1/4-(-1/4) = 0 ... no wait.
    # Actually (1+w^2)^{1/4} cancels in the numerator's (alpha*beta) contribution.
    # The lambda=1/4 family contributes -1/4 and the sqrt contributes +1/4, net 0.
    # The lambda=3/4 family contributes -3/4 and the sqrt contributes +1/4, net -1/2.

    print("\nGF formula: GF(w,u) = 2pi^2 / (sqrt(1+u^2)(sqrt((1+iw)(1+iu))+sqrt((1-iw)(1-iu))))")
    print("The asymmetry sqrt(1+u^2) vs sqrt(1+w^2) reflects lambda_P=1/4 vs lambda_Q=3/4.")

    # ================================================================
    # Zeros of truncations
    # ================================================================
    print("\n" + "="*70)
    print("ZEROS OF TRUNCATIONS sum_{k<=N} E_{0k} Q_hat_k(z)")
    print("="*70)

    Q_all_coeffs = mp_poly_coeffs(max(42, MMAX + 1), lam34)
    zeros_data = {}
    for N in [10, 20, 40]:
        # Build polynomial in z: sum_k E_{0k}^{norm} * Q_k(z) / sqrt(h_k)
        # E_{0k}^{norm} = r_exact[(0,k)] * pi^2 / sqrt(h0_14 * h0_34 * nr14(0) * nr34(k))
        # But we just need to get the polynomial. Let's just use unnormalized:
        # f_N(z) proportional to sum_k r_exact[(0,k)] / sqrt(nr34(k)) * Q_k(z)
        # (the zeros don't depend on overall scale)
        # Actually for correct shape, use normalized: E_{0k} * Q_hat_k = E_{0k}/sqrt(h_k^{3/4}) * Q_k
        # E_{0k} = pi^2 r_{0k} / sqrt(h_0^{1/4} h_k^{3/4})
        # So coefficient of Q_k in the sum = pi^2 r_{0k} / sqrt(h_0^{1/4} h_k^{3/4}) * 1/sqrt(h_k^{3/4})
        # = pi^2 r_{0k} / (sqrt(h_0^{1/4}) * h_k^{3/4})
        # The norms h_k^{3/4} = h_0^{3/4} * norm_ratio_34(k).
        # Zeros only depend on relative coefficients, so use c_k = r_{0k} / norm_ratio_34(k).

        N_use = min(N, MMAX)
        poly = [mp.mpf(0)] * (N_use + 1)
        for k in range(N_use + 1):
            rr = r_exact[(0, k)]
            if rr == 0:
                continue
            ck = mp.mpf(rr.numerator) / mp.mpf(rr.denominator)
            nr = norm_ratio(k, lam34)
            ck /= mp.mpf(nr.numerator) / mp.mpf(nr.denominator)
            for j, c in enumerate(Q_all_coeffs[k]):
                cf = mp.mpf(c.numerator) / mp.mpf(c.denominator)
                poly[j] += ck * cf

        while len(poly) > 1 and abs(poly[-1]) < mp.mpf(10)**(-30):
            poly.pop()

        try:
            mp.mp.dps = 80
            roots = mp.polyroots(poly[::-1], maxsteps=300, extraprec=200)
            mp.mp.dps = 50
            zeros_data[N] = roots
            real_r = sorted([float(z.real) for z in roots if abs(z.imag) < 0.1])
            cmplx = sorted([(float(z.real), float(z.imag)) for z in roots if z.imag > 0.1],
                           key=lambda z: z[0])
            print(f"N={N_use}: {len(real_r)} real, {len(cmplx)} complex pairs")
            if real_r:
                print(f"  Real: {', '.join('%.3f' % r for r in real_r[:12])}")
            if cmplx:
                print(f"  Complex: {', '.join('(%.2f,%.2f)' % z for z in cmplx[:8])}")
        except Exception as e:
            mp.mp.dps = 50
            print(f"N={N_use}: root-finding failed: {e}")

    # Plot zeros
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    colors = {10: '#2166ac', 20: '#b2182b', 40: '#4dac26'}
    markers = {10: 'o', 20: 's', 40: 'D'}
    for N in [10, 20, 40]:
        N_use = min(N, MMAX)
        if N_use in zeros_data or N in zeros_data:
            roots = zeros_data.get(N, zeros_data.get(N_use, []))
            xs = [float(z.real) for z in roots]
            ys = [float(z.imag) for z in roots]
            ax.scatter(xs, ys, s=25, c=colors[N], marker=markers[N],
                      alpha=0.7, label=f'N={N_use}', edgecolors='k', linewidths=0.3)
    ax.set_xlabel('Re(z)', fontsize=12)
    ax.set_ylabel('Im(z)', fontsize=12)
    ax.set_title(r'Zeros of $\sum_{k\leq N} E_{0k}\,\hat{Q}_k(z)$ (joint $\frac{1}{4}\otimes\frac{3}{4}$ cross-matching)', fontsize=13)
    ax.legend(fontsize=11)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0, color='gray', linewidth=0.5)
    ax.grid(True, alpha=0.3)
    fig.savefig("/Users/robin/git/category/joint_sech_zeros.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("\nSaved joint_sech_zeros.png")

    # ================================================================
    # Final summary tables for LaTeX
    # ================================================================
    print("\n" + "="*70)
    print("FINAL TABLES FOR LATEX")
    print("="*70)

    # Table 1: Exact r_mk (nonzero entries only: m+k even)
    print("\n=== r_mk = E_mk^{unnorm} / pi^2 (m+k even only) ===")
    for m in range(8):
        for k in range(m % 2, 8, 2):
            r = r_exact[(m,k)]
            print(f"  ({m},{k}): {r}")

    # Table 2: Normalized E_mk values
    print("\n=== Normalized E_mk table (8x8) ===")
    print("     " + "".join(f"{'k='+str(k):>10s}" for k in range(8)))
    for m in range(8):
        row = []
        for k in range(8):
            if (m+k) % 2 == 1:
                row.append(f"{'0':>10s}")
            else:
                v = float(r_exact[(m,k)]) * mp.pi**2
                v /= float(mp.sqrt(nrm14[m] * nrm34[k]))
                row.append(f"{float(v):10.5f}")
        print(f"  m={m}: " + "".join(row))

    # Table 3: M_mk (sqrt-sech) table
    print("\n=== Sqrt-sech M_mk table (8x8) ===")
    print("     " + "".join(f"{'k='+str(k):>10s}" for k in range(8)))
    for m in range(8):
        row = []
        for k in range(8):
            v = mp.re(M_quad[(m,k)])
            row.append(f"{float(v):10.5f}")
        print(f"  m={m}: " + "".join(row))

    # Row 0 sequence for OEIS
    print("\n=== Row 0 numerators (r_{0,2n} with common denominator extracted) ===")
    row0_nums = []
    for k in range(0, MMAX+1, 2):
        r = r_exact[(0,k)]
        # express as num/2^p
        p = 0
        d = r.denominator
        while d % 2 == 0:
            d //= 2
            p += 1
        row0_nums.append(r.numerator)  # signed numerator
    print("Numerators of r_{0,2n}:", row0_nums)
    print("Denominators of r_{0,2n}:", [r_exact[(0,2*n)].denominator for n in range(len(row0_nums))])

    # Diagonal
    print("\nDiagonal r_{n,n}:", [r_exact[(n,n)] for n in range(min(8, MMAX+1)) if n % 2 == 0])

    # ================================================================
    # Combinatorics: decompose into Euler numbers
    # ================================================================
    print("\n" + "="*70)
    print("COMBINATORIAL DECOMPOSITION INTO EULER NUMBERS")
    print("="*70)
    # r_mk = sum_j (P_m * Q_k)_{2j} * |E_{2j}| / 16^j
    # Show the decomposition for small entries
    for m in range(4):
        for k in range(m % 2, 4, 2):
            deg = m + k
            prod = [Fraction(0)] * (deg + 1)
            for i, a in enumerate(P_coeffs[m]):
                for jj, b in enumerate(Q_coeffs[k]):
                    prod[i + jj] += a * b
            terms = []
            for j in range(0, deg + 1, 2):
                if prod[j] != 0:
                    En = abs(int(mp.eulernum(j)))
                    terms.append(f"{prod[j]}*|E_{j}|/16^{j//2}")
            print(f"  r({m},{k}) = " + " + ".join(terms) + f" = {r_exact[(m,k)]}")


if __name__ == "__main__":
    main()
