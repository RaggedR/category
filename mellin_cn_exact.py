"""Exact evaluation of the Mellin constants c_n in

    M[psi_n(sqrt2 .)](1/2 + 2ix) = c_n * F_m(x),

where psi_n = H_n e^{-x^2/2}/sqrt(2^n n! sqrt(pi)) (Hermite functions),
F_m = Gamma(1/4+ix) MPhat_m (even n = 2m, lambda=1/4) or
F_m = Gamma(3/4+ix) MPhat_m (odd n = 2m+1, lambda=3/4).

Key structural fact (bump_ng.py): for even n,
    M psi_n(s) = norm_n * 2^{s/2-1} Gamma(s/2) * P_n^BN(s),
and the sqrt2 dilation contributes 2^{-s/2}, killing the 2-power.  So the
identity forces  P_n^BN(1/2 + 2ix) = rho * P_m^MP(x)  as POLYNOMIALS,
with rho a constant.  If sympy confirms proportionality, then

    even n=2m:  c = (1/2)      * rho * sqrt(h_m^{1/4}) * norm_n
    odd  n=2m+1: c = 2^{-1/2}  * rho * sqrt(h_m^{3/4}) * norm_n

with h_m^lambda = 2 pi Gamma(m+2lambda) / (2^{2lambda} m!)   (MP norm, phi=pi/2)
and  norm_n = 1/sqrt(2^n n! sqrt(pi)).

Numerical cross-check against the paper's ratios:
n=4: -0.613,  n=6: -0.672i,  n=5: -0.549,  n=7: -0.508i.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
x = sp.symbols('x')
s = sp.symbols('s')
I = sp.I


def bump_ng_poly(n):
    """P_n^BN(s) with the 2^{k/2} factors, as in bump_ng.py."""
    H = sp.hermite(n, x)
    coeffs = sp.Poly(H, x).all_coeffs()[::-1]
    par = n % 2
    P = sp.Integer(0)
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        j = (k - par) // 2
        rising = sp.prod([(s + par) / 2 + i for i in range(j)])
        P += c * sp.sqrt(2)**k * rising
    return sp.expand(P)


def mp_poly(m, lam):
    """MP polynomial P_m^lambda(x), phi=pi/2:
    (j+1)P_{j+1} = 2x P_j - (j + 2lam - 1) P_{j-1}, P_0=1, P_1=2x."""
    p0, p1 = sp.Integer(1), 2*x
    if m == 0:
        return p0
    for j in range(1, m):
        p0, p1 = p1, sp.expand((2*x*p1 - (j + 2*lam - 1)*p0) / (j + 1))
    return p1


def h_mp(m, lam):
    """MP norm h_m = 2 pi Gamma(m+2lam) / (2^{2lam} m!)."""
    return 2 * sp.pi * sp.gamma(m + 2*lam) / (2**(2*lam) * sp.factorial(m))


def main():
    print("=" * 72)
    print("STEP 1: proportionality  P_n^BN(1/2+2ix)  vs  P_m^MP(x)")
    print("=" * 72)
    results = {}
    for n in range(0, 9):
        par = n % 2
        m_idx = (n - par) // 2
        lam = sp.Rational(1, 4) if par == 0 else sp.Rational(3, 4)
        PBN = bump_ng_poly(n).subs(s, sp.Rational(1, 2) + 2*I*x)
        PMP = mp_poly(m_idx, lam)
        ratio = sp.simplify(sp.cancel(sp.expand(PBN) / PMP))
        is_const = x not in ratio.free_symbols
        print(f"n={n} (MP index {m_idx}, lambda={lam}): "
              f"ratio = {sp.nsimplify(ratio)}   constant: {is_const}")
        if not is_const:
            print(f"   PBN = {sp.expand(PBN)}")
            print(f"   PMP = {sp.expand(PMP)}")
            continue
        results[n] = (m_idx, lam, sp.simplify(ratio))

    print()
    print("=" * 72)
    print("STEP 2: exact c_n and numerical comparison")
    print("=" * 72)
    paper = {4: '-0.613', 6: '-0.672i', 5: '-0.549', 7: '-0.508i'}
    exact = {}
    for n, (m_idx, lam, rho) in results.items():
        norm_n = 1 / sp.sqrt(2**n * sp.factorial(n) * sp.sqrt(sp.pi))
        h = h_mp(m_idx, lam)
        # Mellin prefactor is 2^{(s+par)/2-1} * 2^{-s/2} = 2^{par/2 - 1};
        # bump_ng_poly's 2^{k/2} already carries the 2^{par/2}, so a clean
        # factor 1/2 remains for BOTH parities.
        pref = sp.Rational(1, 2)
        c = sp.simplify(pref * rho * sp.sqrt(h) * norm_n)
        exact[n] = c
        cval = complex(c.evalf(30))
        note = f"   [paper: {paper[n]}]" if n in paper else ""
        print(f"n={n}: c_n = {c}\n      = {cval:.15g}{note}")

    print()
    print("=" * 72)
    print("STEP 3: direct numerical verification of the full identity")
    print("=" * 72)
    # M[psi_n(sqrt2 .)](1/2+2ix0) vs c_n * Gamma(lam+ix0) * MPhat_m(x0)
    x0 = mp.mpf('0.3')
    s0 = mp.mpf('0.5') + 2j * x0
    for n, (m_idx, lam, rho) in sorted(results.items()):
        lam_mp = mp.mpf(1)/4 if n % 2 == 0 else mp.mpf(3)/4

        def psi(t, n=n):
            return (mp.hermite(n, t) * mp.exp(-t**2/2)
                    / mp.sqrt(2**n * mp.factorial(n) * mp.sqrt(mp.pi)))

        lhs = mp.quad(lambda t: psi(mp.sqrt(2)*t) * t**(s0-1),
                      [0, 1, 3, 6, 10, 20, 40])

        # MPhat_m at x0
        Pm = mp_poly(m_idx, sp.Rational(1, 4) if n % 2 == 0 else sp.Rational(3, 4))
        Pm_val = mp.mpc(complex(Pm.subs(x, sp.Float(0.3, 30)).evalf(35)))
        h_val = mp.mpc(complex(h_mp(m_idx,
                       sp.Rational(1, 4) if n % 2 == 0 else sp.Rational(3, 4)).evalf(35)))
        F_val = mp.gamma(lam_mp + 1j*x0) * Pm_val / mp.sqrt(h_val)

        c_val = mp.mpc(complex(exact[n].evalf(35)))
        rhs = c_val * F_val
        ok = abs(lhs - rhs) / abs(lhs) < mp.mpf('1e-14')
        print(f"n={n}: LHS = {mp.nstr(lhs, 20)}")
        print(f"      RHS = {mp.nstr(rhs, 20)}   match: {ok}")


if __name__ == "__main__":
    main()
