"""Mellin transforms of Hermite functions: Bump-Ng polynomials and their zeros.

M psi_n(s) = int_0^infty psi_n(x) x^{s-1} dx,  psi_n = c_n H_n(x) e^{-x^2/2}.

Monomial identity: int_0^infty x^{k+s-1} e^{-x^2/2} dx = 2^{(s+k)/2 - 1} Gamma((s+k)/2).

For even n only even k appear -> factor Gamma(s/2) * 2^{s/2-1} * P_n(s).
For odd n only odd k appear  -> factor Gamma((s+1)/2) * 2^{(s+1)/2-1} * P_n(s).
Gamma has no zeros, so zeros of M psi_n = zeros of the polynomial P_n.
"""
import sympy as sp

s = sp.symbols('s')
x = sp.symbols('x')

for n in range(0, 7):
    H = sp.hermite(n, x)  # physicist's Hermite polynomial
    coeffs = sp.Poly(H, x).all_coeffs()[::-1]  # coeffs[k] * x^k
    par = n % 2
    # write Gamma((s+k)/2) = Gamma((s+par)/2) * rising factorial:
    # for k = par + 2j:  Gamma((s+par)/2 + j) = Gamma((s+par)/2) * prod_{i=0}^{j-1} ((s+par)/2 + i)
    P = sp.Integer(0)
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        assert k % 2 == par
        j = (k - par) // 2
        rising = sp.prod([(s + par) / 2 + i for i in range(j)])
        P += c * sp.Rational(2) ** sp.Rational(k, 2) * rising  # 2^{k/2} from 2^{(s+k)/2}
    P = sp.expand(P)
    Pmonic = sp.Poly(P, s)
    print(f"n={n}:  M psi_n(s) = c_{n} * 2^((s+{par})/2 - 1) * Gamma((s+{par})/2) * P_{n}(s)")
    print(f"   P_{n}(s) = {sp.factor(P)}")
    if Pmonic.degree() > 0:
        roots = sp.nroots(Pmonic)
        for r in roots:
            print(f"      zero: s = {sp.re(r):.6f} + {sp.im(r):.6f}i   (Re = {sp.re(r)})")
    # sanity check n=0,1 against direct integration
    print()
