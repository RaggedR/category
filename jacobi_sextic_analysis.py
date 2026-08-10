"""Analysis pass for the sextic test (JACOBI-PLAN.md Task 1, part 2).

Frame: primitive integer polynomials PZ_m = c_m * pi_m in Z[x] (content 1).
Entries E_mk = int PZ_m PZ_k W^{3/2} dx / nu0  in Q.
Hypothesis: den(E_mk) uses ONLY thirds-class primes (from (2/3)_n), so an
s-dependent rescale by B(s) = prod_{j<=s}(3j-1) (Chebyshev's Q(s) analog,
s=(m+k)/2) makes entries integral.

Also verifies the Gegenbauer identification
  beta_n = r_n/r_{n-1} = 9n(9n-10)/((18n-1)(18n-19)),  lambda = -1/18.
"""

import sympy as sp

x = sp.symbols('x')
alpha = sp.Rational(-5, 9)
NMAX = 8

def poch(a, n):
    r = sp.Integer(1)
    for j in range(n):
        r *= (a + j)
    return r

def moment_ratio(expnt, j):
    return poch(sp.Rational(1, 2), j) / poch(expnt + sp.Rational(3, 2), j)

MU = {2 * j: moment_ratio(alpha, j) for j in range(0, 4 * NMAX + 2)}
NU = {2 * j: moment_ratio(sp.Rational(3, 2) * alpha, j) for j in range(0, 4 * NMAX + 2)}

def functional(p, M):
    p = sp.Poly(sp.expand(p), x)
    return sum(c * M[m[0]] for m, c in zip(p.monoms(), p.coeffs()) if m[0] % 2 == 0)

# monic Gram-Schmidt
P = [sp.Integer(1)]
for n in range(1, NMAX + 1):
    pn = x ** n
    for j in range(n):
        pn -= functional(x ** n * P[j], MU) / functional(P[j] * P[j], MU) * P[j]
    P.append(sp.expand(pn))

# Gegenbauer beta check
print("beta_n = r_n/r_{n-1} vs 9n(9n-10)/((18n-1)(18n-19)):")
r = [functional(P[m] * P[m], MU) for m in range(NMAX + 1)]
for n in range(1, NMAX + 1):
    pred = sp.Rational(9 * n * (9 * n - 10), (18 * n - 1) * (18 * n - 19))
    print(f"  n={n}: {r[n] / r[n - 1]} vs {pred}  {'OK' if r[n] / r[n - 1] == pred else 'MISMATCH'}")

# primitive integer polynomials
PZ, content = [], []
for n in range(NMAX + 1):
    den = sp.lcm([sp.Rational(c).q for c in sp.Poly(P[n], x).coeffs()])
    pz = sp.expand(P[n] * den)
    g = sp.gcd([sp.Integer(c) for c in sp.Poly(pz, x).coeffs()])
    PZ.append(sp.expand(pz / g))
    content.append(den / g)
    print(f"PZ_{n} = {PZ[n]}   (c_{n} = {content[n]} = {sp.factorrat(sp.Rational(content[n]))})")

# entries in the primitive-Z frame
print("\nE_mk = <PZ_m, PZ_k>_{W^{3/2}} / nu0, factored:")
E = [[functional(PZ[m] * PZ[k], NU) for k in range(NMAX + 1)] for m in range(NMAX + 1)]
for m in range(NMAX + 1):
    for k in range(m, NMAX + 1):
        if (m + k) % 2 == 1 or E[m][k] == 0:
            continue
        print(f"  E[{m}][{k}] = {E[m][k]} = {sp.factorrat(sp.Rational(E[m][k]))}")

# denominator analysis: thirds-class only?
print("\nden(E_mk) factored; thirds hypothesis check (den | 2^a 3^b prod(3j-1), j<=s):")
def Bclass(s):  # prod_{j=1..s} (3j-1) = 2*5*8*11*...
    r = sp.Integer(1)
    for j in range(1, s + 1):
        r *= (3 * j - 1)
    return r

all_ok = True
for m in range(NMAX + 1):
    for k in range(m, NMAX + 1):
        if (m + k) % 2 == 1 or E[m][k] == 0:
            continue
        s = (m + k) // 2
        den = sp.Rational(E[m][k]).q
        # strip 2s and 3s
        rem = den
        for p_ in (2, 3):
            while rem % p_ == 0:
                rem //= p_
        ok = (Bclass(s + 2) % rem == 0)  # allow slack of 2 in s
        if not ok:
            all_ok = False
        print(f"  ({m},{k}) s={s}: den={sp.factorint(den)}  odd-part-divides-B(s+2)? {ok}")
print(f"\nthirds hypothesis (with slack 2): {'HOLDS' if all_ok else 'FAILS'}")

# rescale search: D_mk = E_mk * 3^? * B(s+c) minimal clearing
print("\nrescaled D_mk = E_mk * (minimal 2^a 3^b B-product):  row-0/diag/edge sequences")
for m in range(NMAX + 1):
    row = []
    for k in range(NMAX + 1):
        if (m + k) % 2 == 1:
            continue
        s = (m + k) // 2
        v = E[m][k] * Bclass(s)
        row.append((k, v))
    print(f"  row {m} (xB(s)): {row}")
