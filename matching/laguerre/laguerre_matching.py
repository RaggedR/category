"""
Laguerre Matching Matrix
========================
For weight W(x) = x^alpha * e^{-x} on (0, inf), orthonormal Laguerre polynomials L̂_n^(alpha),
the matching matrix is:

  b_mk = ∫_0^inf L̂_m(x) L̂_k(x) x^{3*alpha/2} e^{-3x/2} dx

Classical closed form (alpha=0):
  ∫_0^inf L_m(x) L_k(x) e^{-p*x} dx  (p > 1)
  = (p-1)^{m+k} / p^{m+k+1}  * sum ...
  Actually: ∫_0^inf e^{-px} L_m(x) L_k(x) dx = (-1)^{m+k} (p-1)^{m+k}/p^{m+k+1} * P_k^{m-k}((2p-1)...

  More precisely (DLMF 18.17.34 or Gradshteyn & Ryzhik):
  ∫_0^inf e^{-sx} L_m(x) L_k(x) dx = (s-1)^{m+k} / s^{m+k+1} * C(m+k, m) * 2F1(-m,-k;1;s^2/(s-1)^2)

  Actually the exact formula is:
  ∫_0^inf e^{-sx} L_m(x) L_k(x) dx = C(m+k,m) * (s-1)^{m+k} / s^{m+k+1}   [for s > 1, m,k >=0 integers]

  No wait. Let me use the generating function approach directly.
  GF: sum_{m,k} ∫ L_m L_k e^{-3x/2} w^m u^k dx = ∫_0^inf sum_m L_m w^m * sum_k L_k u^k * e^{-3x/2} dx
     = ∫_0^inf e^{-x*w/(1-w)} / (1-w) * e^{-x*u/(1-u)} / (1-u) * e^{-3x/2} dx
     = 1/((1-w)(1-u)) * ∫_0^inf exp(-x * [3/2 + w/(1-w) + u/(1-u)]) dx
     = 1/((1-w)(1-u)) * 1 / [3/2 + w/(1-w) + u/(1-u)]
     = 1/((1-w)(1-u)) * (1-w)(1-u) / [3(1-w)(1-u)/2 + w(1-u) + u(1-w)]
     = 1 / [3(1-w)(1-u)/2 + w(1-u) + u(1-w)]

  Numerator: 1
  Denominator: 3/2 - 3w/2 - 3u/2 + 3wu/2 + w - wu + u - wu
             = 3/2 - w/2 - u/2 - wu/2
             = (3 - w - u - wu) / 2

  So GF = 2 / (3 - w - u - wu)
         = 2 / (3 - w - u - wu)

  This is a RATIONAL generating function!

  For alpha=0: b_mk = ∫ L̂_m L̂_k e^{-3x/2} dx where L̂_n = L_n (since norms = 1 for alpha=0)

  Actually norm of L_n^(0) is Gamma(n+1)/n! = 1, so L̂_n^(0) = L_n^(0).

  b_mk = coefficient of w^m u^k in GF * (normalization factors)
       Since L̂_m = L_m for alpha=0, b_mk = [w^m u^k] * 2/(3-w-u-wu)

For alpha != 0:
  Laguerre GF: sum_n L_n^(alpha)(x) t^n = (1-t)^{-alpha-1} exp(-xt/(1-t))

  ∫_0^inf x^{3alpha/2} e^{-3x/2} * (1/(1-w))^{alpha+1} exp(-xw/(1-w)) * (1/(1-u))^{alpha+1} exp(-xu/(1-u)) dx
  = 1/((1-w)^{alpha+1}(1-u)^{alpha+1}) * ∫_0^inf x^{3alpha/2} exp(-x[3/2 + w/(1-w) + u/(1-u)]) dx
  = 1/((1-w)^{alpha+1}(1-u)^{alpha+1}) * Gamma(3alpha/2 + 1) / [3/2 + w/(1-w) + u/(1-u)]^{3alpha/2+1}
  = Gamma(3alpha/2+1) / ((1-w)^{alpha+1}(1-u)^{alpha+1}) * ((1-w)(1-u))^{3alpha/2+1} / ((3-w-u-wu)/2)^{3alpha/2+1}
  = Gamma(3alpha/2+1) * 2^{3alpha/2+1} * (1-w)^{3alpha/2-alpha}(1-u)^{3alpha/2-alpha} / (3-w-u-wu)^{3alpha/2+1}
  = Gamma(3alpha/2+1) * 2^{3alpha/2+1} * (1-w)^{alpha/2}(1-u)^{alpha/2} / (3-w-u-wu)^{3alpha/2+1}

  This needs to be multiplied by the norm factors: L̂_n^(alpha) = L_n^(alpha) / sqrt(Gamma(n+alpha+1)/n!)

  So: sum_{m,k} ∫ L_m^(alpha) L_k^(alpha) x^{3alpha/2} e^{-3x/2} dx * w^m u^k
    = Gamma(3alpha/2+1) * 2^{3alpha/2+1} * (1-w)^{alpha/2}(1-u)^{alpha/2} / (3-w-u-wu)^{3alpha/2+1}

  For alpha=0: Gamma(1) * 2 * 1 / (3-w-u-wu)^1 = 2/(3-w-u-wu) ✓
"""

import mpmath
import sympy as sp
import numpy as np
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import product
from fractions import Fraction

mpmath.mp.dps = 50

# ============================================================
# PART 1: Alpha=0 exact rational table via GF
# ============================================================

def gf_coefficients_alpha0(N):
    """
    Compute b_mk = [w^m u^k] in 2/(3-w-u-wu) for 0 <= m,k <= N-1.
    Use exact rational arithmetic.
    """
    # We need Taylor coefficients of 2/(3 - w - u - wu)
    # Let f(w,u) = 2/(3-w-u-wu)
    # f = (2/3) * 1/(1 - (w+u+wu)/3)
    # = (2/3) * sum_{n>=0} ((w+u+wu)/3)^n

    # Coefficient of w^m u^k in ((w+u+wu)/3)^n:
    # (w+u+wu)^n = sum_{a+b+c=n} n!/(a!b!c!) * w^a * u^b * (wu)^c
    #            = sum_{a+b+c=n} n!/(a!b!c!) * w^{a+c} * u^{b+c}
    # We want a+c = m, b+c = k, so a = m-c, b = k-c, and a,b,c >= 0
    # => c ranges from 0 to min(m,k), m-c >= 0, k-c >= 0
    # Then n = a+b+c = m+k-c

    # So [w^m u^k] in (w+u+wu)^n:
    # = sum_{c: n=m+k-c, 0<=c<=min(m,k)} n!/(m-c)!(k-c)!c!
    # = (for fixed n, this is nonzero only if n <= m+k)

    # [w^m u^k] f(w,u) = (2/3) * sum_{n>=0} 1/3^n * [w^m u^k](w+u+wu)^n
    # Only n = m+k-c contributes for c = m+k-n, i.e., one c per n.
    # Actually: [w^m u^k] in ((w+u+wu)/3)^n is nonzero only for specific n.
    # For fixed m,k:
    # [w^m u^k] f = (2/3) * sum_{c=0}^{min(m,k)} 1/3^{m+k-c} * (m+k-c)!/((m-c)!(k-c)!c!)
    # = (2/3) * (1/3^{m+k}) * sum_{c=0}^{min(m,k)} 3^c * (m+k-c)!/((m-c)!(k-c)!c!)

    table = [[Fraction(0)] * N for _ in range(N)]

    for m in range(N):
        for k in range(N):
            val = Fraction(0)
            for c in range(min(m, k) + 1):
                n = m + k - c
                from math import factorial
                coeff = Fraction(factorial(n), factorial(m-c) * factorial(k-c) * factorial(c))
                val += coeff * Fraction(3**c, 3**(m+k))
            val *= Fraction(2, 3)
            table[m][k] = val

    return table

# ============================================================
# PART 2: Numerical verification and normalization
# ============================================================

def laguerre_poly(n, x, alpha=0):
    """Laguerre polynomial L_n^(alpha)(x) using recurrence."""
    if n == 0:
        return mpmath.mpf(1)
    elif n == 1:
        return 1 + alpha - x
    else:
        L_prev2 = mpmath.mpf(1)
        L_prev1 = 1 + alpha - x
        for k in range(2, n+1):
            L_curr = ((2*k-1+alpha-x)*L_prev1 - (k-1+alpha)*L_prev2) / k
            L_prev2 = L_prev1
            L_prev1 = L_curr
        return L_curr

def laguerre_norm_sq(n, alpha=0):
    """||L_n^(alpha)||^2 = Gamma(n+alpha+1)/n! under weight x^alpha e^{-x}"""
    return mpmath.gamma(n + alpha + 1) / mpmath.fac(n)

def compute_bmk_numeric(m, k, alpha=0, p=1.5, npts=200):
    """
    Numerical integration of b_mk = ∫_0^inf L̂_m L̂_k x^{3*alpha/2} e^{-3x/2} dx
    """
    norm_m = mpmath.sqrt(laguerre_norm_sq(m, alpha))
    norm_k = mpmath.sqrt(laguerre_norm_sq(k, alpha))

    def integrand(x):
        Lm = laguerre_poly(m, x, alpha)
        Lk = laguerre_poly(k, x, alpha)
        return Lm * Lk * x**(3*alpha/2) * mpmath.exp(-p*x) / (norm_m * norm_k)

    result = mpmath.quad(integrand, [0, mpmath.inf], maxdegree=8, error=True)
    return result[0]

# ============================================================
# PART 3: Generating function expansion for general alpha
# ============================================================

def gf_coefficients_alpha_sympy(alpha_val, N):
    """
    GF for general alpha (unnormalized):
    G(w,u) = Gamma(3*alpha/2+1) * 2^{3*alpha/2+1} * (1-w)^{alpha/2} * (1-u)^{alpha/2} / (3-w-u-wu)^{3*alpha/2+1}

    Extract Taylor coefficients [w^m u^k] G(w,u) and divide by norm factors.
    """
    import sympy as sp
    w, u = sp.symbols('w u')
    a = sp.Rational(alpha_val)  # alpha as rational

    # Compute GF symbolically
    exponent_top = sp.Rational(3*alpha_val, 2)  # 3*alpha/2

    prefactor = sp.gamma(exponent_top + 1) * 2**(exponent_top + 1)

    # Expand (1-w)^{alpha/2} * (1-u)^{alpha/2} around w=u=0
    # = sum_{a,b} C(alpha/2, a) C(alpha/2, b) (-w)^a (-u)^b
    # (3-w-u-wu)^{-(3*alpha/2+1)} needs expansion too

    # For alpha=0: straightforward. For alpha=±1/2: irrational exponents.
    # Use mpmath for numerical extraction.
    pass

def extract_gf_coeffs_numeric(alpha_val, N, dps=50):
    """
    Numerically extract Taylor coefficients of the GF:
    G(w,u) = Gamma(3*alpha/2+1) * 2^{3*alpha/2+1} * (1-w)^{alpha/2} * (1-u)^{alpha/2} / (3-w-u-wu)^{3*alpha/2+1}

    These are [w^m u^k] G(w,u) = unnormalized integral of L_m L_k x^{3*alpha/2} e^{-3x/2}
    """
    mpmath.mp.dps = dps
    a = mpmath.mpf(alpha_val)

    table = np.zeros((N, N), dtype=object)

    prefactor = mpmath.gamma(3*a/2 + 1) * mpmath.power(2, 3*a/2 + 1)

    # Use numerical differentiation via contour integrals (Cauchy)
    # Or: expand using series/multinomial

    # For rational alpha, we can use symbolic expansion
    # For alpha=0: exact rational
    # For alpha=±1/2: extract numerically

    def G_func(w_val, u_val):
        return prefactor * (1-w_val)**(a/2) * (1-u_val)**(a/2) / (3-w_val-u_val-w_val*u_val)**(3*a/2+1)

    # Extract via 2D Taylor series using finite differences
    # Better: use univariate series extraction
    # [w^m] G(w, u) is itself a function of u; then [u^k] of that.

    # Use mpmath's taylor series
    results = {}

    for m in range(N):
        for k in range(N):
            def f_mk(eps_w, eps_u, m=m, k=k):
                return G_func(eps_w, eps_u)

            # Use numerical diff: Taylor coefficient via contour
            r = mpmath.mpf('0.3')  # radius
            def integrand_m(theta_m, k=k, m=m):
                w_val = r * mpmath.exp(mpmath.j * theta_m)
                def integrand_k(theta_k):
                    u_val = r * mpmath.exp(mpmath.j * theta_k)
                    return G_func(w_val, u_val) / (w_val**(m) * u_val**(k)) / (4 * mpmath.pi**2)
                return mpmath.quad(integrand_k, [0, 2*mpmath.pi])

            val = mpmath.quad(integrand_m, [0, 2*mpmath.pi])
            results[(m, k)] = val

    return results

# ============================================================
# PART 4: PSLQ for alpha = ±1/2
# ============================================================

def pslq_identify(val, constants, names, tol=1e-30, max_coeff=10000):
    """Try to identify val as integer linear combination of constants."""
    vec = [val] + list(constants)
    try:
        rel = mpmath.identify(val, tol=tol)
        if rel:
            return rel
    except:
        pass

    try:
        coeffs = mpmath.pslq(vec, tol=tol, maxcoeff=max_coeff)
        if coeffs is not None:
            # val * coeffs[0] + sum(c_i * const_i) = 0
            # val = -sum(c_i * const_i) / coeffs[0]
            result = {names[i]: -coeffs[i+1] for i in range(len(constants))}
            denom = coeffs[0]
            return result, denom
    except:
        pass
    return None

# ============================================================
# PART 5: Zeros of truncations
# ============================================================

def compute_zeros(row, alpha=0, N_trunc=20):
    """
    Compute zeros of truncation sum_{k<=N} b_{row,k} L̂_k(x)
    by discretizing on a grid and using polynomial root finding.
    """
    # Actually, express as a polynomial via Laguerre polynomial coefficients
    # L_n(x) = sum_{j=0}^n C(n,j) (-x)^j / j!
    # So sum_k b_k L̂_k(x) is a polynomial in x (of degree N_trunc)

    mpmath.mp.dps = 40

    # Get b_{row, k} coefficients
    if alpha == 0:
        table = gf_coefficients_alpha0(N_trunc + 1)
        bmk_coeffs = [float(table[row][k]) for k in range(N_trunc + 1)]
    else:
        # Numerical
        bmk_coeffs = []
        for k in range(N_trunc + 1):
            v = compute_bmk_numeric(row, k, alpha=alpha)
            bmk_coeffs.append(float(v))

    # Build polynomial: sum_k b_k * L_k(x)
    # Represent as polynomial in x
    from math import comb, factorial

    poly = np.zeros(N_trunc + 1)

    for k in range(N_trunc + 1):
        bk = bmk_coeffs[k]
        # L_k(x) = sum_{j=0}^k C(k,j) (-1)^j / j! * x^j
        norm_sq = float(mpmath.gamma(k + alpha + 1) / mpmath.fac(k))
        norm_k = np.sqrt(norm_sq)
        for j in range(k + 1):
            coeff = comb(k, j) * ((-1)**j) / float(factorial(j))
            poly[j] += bk * coeff / norm_k

    # Find roots of the polynomial
    poly_coeffs = poly[::-1]  # highest degree first for numpy
    roots = np.roots(poly_coeffs)

    return roots, bmk_coeffs

# ============================================================
# MAIN COMPUTATION
# ============================================================

def main():
    print("=" * 60)
    print("LAGUERRE MATCHING MATRIX COMPUTATION")
    print("=" * 60)

    N = 12  # table size

    # --------------------------------------------------------
    # TASK 1 & 2: Alpha=0 exact rational table + GF
    # --------------------------------------------------------
    print("\n--- TASK 1: Exact rational table, alpha=0 ---")

    table_alpha0 = gf_coefficients_alpha0(N)

    print("\nb_mk (rational) for alpha=0, first 8x8:")
    for m in range(8):
        row = []
        for k in range(8):
            row.append(str(table_alpha0[m][k]))
        print(f"  m={m}: {row}")

    # Verify numerically
    print("\nNumerical verification (first 4x4):")
    for m in range(4):
        for k in range(4):
            exact = float(table_alpha0[m][k])
            numerical = float(compute_bmk_numeric(m, k, alpha=0))
            print(f"  b_{m}{k}: exact={exact:.6f}, numerical={numerical:.6f}, match={abs(exact-numerical)<1e-10}")

    # --------------------------------------------------------
    # TASK 2: Identify decay rate
    # --------------------------------------------------------
    print("\n--- TASK 2: Decay rate analysis ---")

    # GF = 2/(3-w-u-wu)
    # b_{0,k} = [w^0 u^k] 2/(3-w-u-wu) = [u^k] 2/(3-u) = 2 * (1/3)^{k+1}... wait
    # Actually at w=0: 2/(3-0-u-0*u) = 2/(3-u) = (2/3) * 1/(1-u/3) = (2/3) sum (u/3)^k
    # So b_{0,k} = (2/3) * (1/3)^k = 2/3^{k+1}

    print("\nb_{0,k} = 2/3^{k+1}:")
    for k in range(8):
        exact = table_alpha0[0][k]
        predicted = Fraction(2, 3**(k+1))
        print(f"  k={k}: {exact} = {predicted}, match={exact==predicted}")

    # Ratio b_{0,k}/b_{0,k-1}:
    print("\nDecay ratio b_{0,k}/b_{0,k-1}:")
    for k in range(1, 8):
        if table_alpha0[0][k-1] != 0:
            ratio = table_alpha0[0][k] / table_alpha0[0][k-1]
            print(f"  k={k}: ratio = {ratio} = {float(ratio):.6f}")

    # --------------------------------------------------------
    # TASK 3: Integer normalization
    # --------------------------------------------------------
    print("\n--- TASK 3: Integer structure ---")

    # Multiply b_mk by 3^{m+k+1}/2 to get integers?
    print("\n3^{m+k+1}/2 * b_mk (should be integer):")
    for m in range(6):
        row = []
        for k in range(6):
            val = table_alpha0[m][k] * Fraction(3**(m+k+1), 2)
            row.append(str(val))
        print(f"  m={m}: {row}")

    # Alternative: look at b_mk * 3^{m+k+1}/2
    # b_mk = 2/3^{m+k+1} * [coeff from formula]
    # Let's check: is b_mk = 2/(m+k choose m) * C(m+k+1, m+1) * 1/3^{m+k+1}?
    # Actually from GF: b_mk = [w^m u^k] 2/(3-w-u-wu)
    # = 2 * [w^m u^k] 1/3 * 1/(1-(w+u+wu)/3)
    # = (2/3) * sum_{n>=0} 1/3^n [w^m u^k](w+u+wu)^n
    # = (2/3) * 1/3^{m+k} * sum_{c=0}^{min(m,k)} 3^c * (m+k-c)!/((m-c)!(k-c)!c!)

    # Define E_mk = 3^{m+k+1}/2 * b_mk = sum_{c=0}^{min(m,k)} 3^c * (m+k-c)!/((m-c)!(k-c)!c!)
    # This should be an integer!

    from math import factorial
    print("\nE_mk = 3^{m+k+1}/2 * b_mk (integers):")
    E_table = []
    for m in range(N):
        row = []
        for k in range(N):
            E = Fraction(0)
            for c in range(min(m,k)+1):
                E += Fraction(3**c * factorial(m+k-c), factorial(m-c)*factorial(k-c)*factorial(c))
            assert E.denominator == 1, f"Not integer at m={m},k={k}: {E}"
            row.append(int(E))
        E_table.append(row)

    print("E_mk table (first 8x8):")
    for m in range(8):
        print(f"  m={m}: {E_table[m][:8]}")

    # Diagonal: E_nn = sum_{c=0}^n 3^c * (2n-c)!/((n-c)!^2 * c!)
    print("\nDiagonal E_nn:")
    for n in range(10):
        print(f"  n={n}: E_{n}{n} = {E_table[n][n]}")

    # Check combinatorial formula: Lah numbers involved?
    # E_mk is a Jacobi polynomial or Legendre at a specific point?
    # E_mk = sum_c C(m,c)*C(k,c)*c! * (m+k-2c choose m-c) * ... hmm
    # Actually E_mk = sum_c 3^c * (m+k-c)!/((m-c)!(k-c)!c!)
    #              = C(m+k,m) * sum_c C(m,c)*C(k,c)*c! * 3^c / C(m+k,m)?
    # Let's factor out C(m+k,m):
    from math import comb
    print("\nE_mk / C(m+k,m):")
    for m in range(5):
        for k in range(5):
            c_mk = comb(m+k, m)
            ratio = Fraction(E_table[m][k], c_mk)
            print(f"  ({m},{k}): {ratio}")

    # --------------------------------------------------------
    # TASK 4: Alpha = ±1/2
    # --------------------------------------------------------
    print("\n--- TASK 4: Alpha = ±1/2 ---")

    for alpha_val in [-0.5, 0.5]:
        print(f"\nalpha = {alpha_val}:")
        print(f"GF factor: Gamma({3*alpha_val/2}+1) * 2^{3*alpha_val/2+1} * (1-w)^{alpha_val/2} * (1-u)^{alpha_val/2} / (3-w-u-wu)^{3*alpha_val/2+1}")

        # Compute a small table numerically
        N_small = 6
        bmk_alpha = np.zeros((N_small, N_small))
        for m in range(N_small):
            for k in range(N_small):
                v = compute_bmk_numeric(m, k, alpha=alpha_val)
                bmk_alpha[m][k] = float(v)

        print(f"Numerical b_mk for alpha={alpha_val} (first 6x6):")
        for m in range(N_small):
            row_str = [f"{bmk_alpha[m][k]:.6f}" for k in range(N_small)]
            print(f"  m={m}: {row_str}")

        # PSLQ to find structure
        # For alpha=1/2: involves Gamma(7/4) type constants?
        # sqrt(pi), Gamma(1/4), etc.
        sqrt_pi = float(mpmath.sqrt(mpmath.pi))
        gamma_34 = float(mpmath.gamma(mpmath.mpf('3/4')))
        gamma_14 = float(mpmath.gamma(mpmath.mpf('1/4')))

        print(f"\nDecay ratio b_{{0,k}}/b_{{0,k-1}} for alpha={alpha_val}:")
        for k in range(1, min(6, N_small)):
            if abs(bmk_alpha[0][k-1]) > 1e-30:
                ratio = bmk_alpha[0][k] / bmk_alpha[0][k-1]
                print(f"  k={k}: ratio = {ratio:.8f}")

        # Try to identify b_{0,0} via PSLQ
        v00 = mpmath.mpf(bmk_alpha[0][0])
        print(f"\nPSLQ on b_{{0,0}} = {float(v00):.10f}:")
        # b_{0,0} for alpha=1/2: ∫_0^inf x^{3/4} e^{-3x/2} / (Gamma(3/2)/0!)...
        # = Gamma(7/4) / (3/2)^{7/4} / norm^2
        # norm^2 for n=0, alpha=1/2: Gamma(3/2) = sqrt(pi)/2
        # b_{0,0} = ∫_0^inf L̂_0^(1/2) L̂_0^(1/2) x^{3/4} e^{-3x/2} dx
        #         = ∫_0^inf x^{3/4} e^{-3x/2} dx / Gamma(3/2)
        #         = Gamma(7/4) / (3/2)^{7/4} / (sqrt(pi)/2)
        predicted_00_half = float(mpmath.gamma(mpmath.mpf('7/4')) / (mpmath.mpf('3/2')**mpmath.mpf('7/4')) / mpmath.gamma(mpmath.mpf('3/2')))
        print(f"  Predicted (Gamma(7/4)/(3/2)^(7/4)/Gamma(3/2)): {predicted_00_half:.10f}")
        print(f"  Matches: {abs(float(v00) - predicted_00_half) < 1e-8}")

    # --------------------------------------------------------
    # TASK 5: Consistency with Hermite (alpha=0 vs alpha=±1/2)
    # --------------------------------------------------------
    print("\n--- TASK 5: Hermite consistency ---")

    # Under x = t^2, L_n^(-1/2)(t^2) ~ H_{2n}(t) (even Hermite)
    # and L_n^(1/2)(t^2) ~ H_{2n+1}(t) / t (odd Hermite)

    # Hermite b_mk from GF: sqrt(2pi/3) exp(-(s^2-4su+u^2)/3)
    # I_{mk}^H = m!k! [s^m u^k] sqrt(2pi/3) exp(-(s^2-4su+u^2)/3)

    print("\nHermite matching matrix (even subblock, 3x3):")
    # D_{2m,2k} from Hermite = ?
    # I_{2m,2k} = sqrt(2pi/3) * sum_{a+b+c: 2m=2a+b, 2k=b+2c} (-1/3)^a (4/3)^b (-1/3)^c * (2m)!(2k)!/(a!b!c!)

    def hermite_Imk(m, k):
        """Compute I_mk for Hermite by GF expansion."""
        from math import factorial
        total = Fraction(0)
        for b in range(min(m,k)*2+1):
            if (m - b) % 2 == 0 and (k - b) % 2 == 0:
                a = (m - b) // 2
                c = (k - b) // 2
                if a >= 0 and c >= 0:
                    coeff = Fraction((-1)**(a+c) * 4**b * factorial(m) * factorial(k),
                                    3**(a+b+c) * factorial(a) * factorial(b) * factorial(c))
                    total += coeff
        return total

    print("  Hermite D_{2m,2k} (ratio to sqrt(2pi/3), scaled by 3^{m+k}):")
    for m in range(3):
        for k in range(3):
            D = 3**(m+k) * hermite_Imk(2*m, 2*k)
            print(f"    D_H({2*m},{2*k}) = {D} (float: {float(D):.4f})")

    print("\nLaguerre b_{m,k} for alpha=-1/2 (should relate to even Hermite):")
    for m in range(3):
        for k in range(3):
            v = compute_bmk_numeric(m, k, alpha=-0.5)
            print(f"    b_L({m},{k}) = {float(v):.6f}")

    # --------------------------------------------------------
    # TASK 6: Zeros of truncations
    # --------------------------------------------------------
    print("\n--- TASK 6: Zeros of truncations ---")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Zeros of Laguerre Matching Truncations', fontsize=14)

    trunc_sizes = [10, 20, 40]
    rows_to_plot = [0, 4]

    all_zeros_data = {}

    for ax_row, row_idx in enumerate(rows_to_plot):
        for ax_col, N_trunc in enumerate(trunc_sizes):
            ax = axes[ax_row][ax_col]

            try:
                roots, coeffs = compute_zeros(row_idx, alpha=0, N_trunc=N_trunc)
                all_zeros_data[(row_idx, N_trunc)] = roots

                real_roots = [r.real for r in roots]
                imag_roots = [r.imag for r in roots]

                # Color by whether they're near real axis
                colors = ['blue' if abs(r.imag) < 0.1 else 'red' for r in roots]

                ax.scatter(real_roots, imag_roots, c=colors, s=20, alpha=0.7)
                ax.axhline(y=0, color='k', linewidth=0.5)
                ax.axvline(x=0, color='k', linewidth=0.5)
                ax.set_title(f'Row {row_idx}, N={N_trunc}')
                ax.set_xlabel('Re(z)')
                ax.set_ylabel('Im(z)')
                ax.grid(True, alpha=0.3)

                # Mark positive real axis
                real_pos = [r.real for r in roots if abs(r.imag) < 0.1 and r.real > 0]
                print(f"  Row {row_idx}, N={N_trunc}: {len(real_pos)} real positive zeros, {len(roots)} total")

            except Exception as e:
                print(f"  Error for row={row_idx}, N={N_trunc}: {e}")
                ax.set_title(f'Row {row_idx}, N={N_trunc} (error)')

    plt.tight_layout()
    plt.savefig('/Users/robin/git/category/laguerre_zeros.png', dpi=120, bbox_inches='tight')
    plt.close()
    print("\nSaved: laguerre_zeros.png")

    # --------------------------------------------------------
    # Summary and caching
    # --------------------------------------------------------
    print("\n--- SUMMARY ---")
    print(f"GF for alpha=0: 2/(3-w-u-wu) [RATIONAL!]")
    print(f"Decay ratio: b_{{0,k}}/b_{{0,k-1}} = 1/3 (exact)")
    print(f"b_{{0,k}} = 2/3^{{k+1}}")

    # Identify E_mk combinatorially
    # E_mk = sum_{c=0}^{min(m,k)} 3^c * (m+k-c)!/((m-c)!(k-c)!c!)
    # = sum_c C(m,c)*C(k,c)*c! * (m+k-2c choose m-c) * ... no
    # Let's compare with Lah numbers and rook numbers
    # Lah number L(n,k) = C(n-1,k-1)*n!/k! counts partitions of [n] into k lists
    # But our formula: (m+k-c)!/((m-c)!(k-c)!c!) = C(m+k-c, m-c, k-c, c) (multinomial)
    #                = C(m+k-c, c) * C(m+k-2c, m-c)

    # Alternative: E_mk = P_k^{(m-k,m-k)}(something) ? Jacobi poly?
    # At w=1, u=0: sum_m E_mm / 3^{2m+1}/2 ...

    # The GF 2/(3-w-u-wu) factors as:
    # 2/(3-w-u-wu) = 2/((1+w)(3-u) - u(1+w) + u) hmm
    # = 2/((3-u)(1+w) - u) no
    # = 2/((1+u)(3-w) - 2u) no
    # Actually: 3-w-u-wu = 3-(1+u)w-u = 3-u-(1+u)w
    # So [w^m u^k] 2/(3-w-u-wu):
    # = 2 * [w^m u^k] 1/((3-u) - (1+u)w)
    # = 2/(3-u) * [w^m u^k] 1/(1 - (1+u)/(3-u) * w)
    # = 2/(3-u) * (1+u)^m/(3-u)^m * [u^k]
    # = 2 * [u^k] (1+u)^m / (3-u)^{m+1}

    print("\nAlternate factorization: b_{m,k} = 2 * [u^k] (1+u)^m / (3-u)^{m+1}")
    print("Verify:")
    # (1+u)^m = sum_j C(m,j) u^j
    # 1/(3-u)^{m+1} = 1/3^{m+1} * sum_l C(m+l, l) (u/3)^l
    # [u^k] = sum_{j+l=k} C(m,j) * C(m+l,l) / 3^{m+1+l}
    #       = sum_{j=0}^{min(m,k)} C(m,j) * C(m+k-j, k-j) / 3^{m+k-j+1}
    for m in range(4):
        for k in range(4):
            val = Fraction(0)
            for j in range(min(m,k)+1):
                val += Fraction(comb(m,j) * comb(m+k-j, k-j), 3**(m+k-j+1))
            val *= 2
            exact = table_alpha0[m][k]
            print(f"  b_{m}{k}: formula={float(val):.6f}, exact={float(exact):.6f}, match={val==exact}")

    # So E_mk = 3^{m+k+1}/2 * b_mk = sum_{j=0}^{min(m,k)} C(m,j)*C(m+k-j,k-j)*3^j
    # = sum_{j=0}^{min(m,k)} C(m,j)*C(m+k-j,m)*3^j
    print("\nE_mk = sum_{j=0}^{min(m,k)} C(m,j)*C(m+k-j,m)*3^j")
    print("Verify:")
    for m in range(4):
        for k in range(4):
            val = sum(comb(m,j)*comb(m+k-j,m)*3**j for j in range(min(m,k)+1))
            print(f"  E_{m}{k}: formula={val}, table={E_table[m][k]}, match={val==E_table[m][k]}")

    # This is a Jacobi polynomial!
    # sum_{j=0}^m C(m,j)*C(m+k-j,m)*3^j = P_m^{(k,0)}(1+2*3) = P_m^{(k,0)}(7)?
    # Actually: C(m+k-j,m) * C(m,j) = C(m+k,m) * C(m+k-j,k) / C(m+k,k) hmm
    # Jacobi: P_m^{(alpha,beta)}(x) = sum C(m+alpha,m-j)*C(m+beta,j) ((x-1)/2)^j ((x+1)/2)^{m-j}

    # Actually E_{m,k} = P_k^{(m,0)}(2) * m! / [something]? Let's check Jacobi at 2:
    # P_k^{(m,0)}(2) = C(m+k,k) * 3F2 or by formula
    import sympy as sp
    x = sp.Symbol('x')
    print("\nJacobi polynomial check P_k^(m,0)(x=2):")
    for m in range(3):
        for k in range(3):
            # P_k^{(m,0)}(2) = C(m+k,k) * (something)
            # Standard: P_k^{(a,b)}(1) = C(k+a, k)
            # P_k^{(a,b)}(-1) = (-1)^k C(k+b,k)
            # At x=2: use the formula
            # P_k^{(m,0)}(x) = sum_{j=0}^k C(k+m,k-j)*C(k,j)*((x-1)/2)^j
            val_jacobi = sum(comb(k+m,k-j)*comb(k,j)*(1/2)**j for j in range(k+1))  # at x=2: (x-1)/2=1/2
            print(f"  P_{k}^({m},0)(2) = {val_jacobi:.4f}, E_{m}{k}={E_table[m][k]}, ratio={E_table[m][k]/val_jacobi if val_jacobi != 0 else 'inf':.4f}")

    # Hmm, let me check the relation more carefully
    # E_mk = sum_{j=0}^{min(m,k)} C(m,j)*C(m+k-j,m)*3^j
    # Compare: P_m^{(k,0)}(2) = sum_{j=0}^m C(m+k,m-j)*C(m,j)*(1/2)^j (using x=2, so (x-1)/2=1/2, (x+1)/2=3/2)
    # Hmm let me just check if E_mk = C(m+k,m) * sum ... or something cleaner

    print("\nE_mk / C(m+k,k):")
    for m in range(4):
        for k in range(4):
            ratio = E_table[m][k] / comb(m+k,k)
            print(f"  E_{m}{k}/C({m+k},{k}) = {ratio:.6f}")

    # Save results
    results_data = {
        'table_alpha0': table_alpha0,
        'E_table': E_table,
        'N': N,
    }
    with open('/Users/robin/git/category/laguerre_matching.pkl', 'wb') as f:
        pickle.dump(results_data, f)
    print("\nSaved: laguerre_matching.pkl")

    return table_alpha0, E_table, all_zeros_data

if __name__ == '__main__':
    table_alpha0, E_table, zeros_data = main()
