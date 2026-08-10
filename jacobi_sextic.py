"""Sextic test for the general-Jacobi matching conjecture (JACOBI-PLAN.md Task 1).

Weight W = (1-x^2)^alpha with alpha = -5/9, polynomials = Gegenbauer
C^(lambda), lambda = alpha + 1/2 = -1/18 (we work with the MONIC family,
exact rational Gram-Schmidt from moments).

Forward matching coefficients  b_mk = int P^_m P^_k W^{3/2} dx.

Theta substitution x = cos(theta):
  entries use Wallis power p  = 3*alpha + 1 = -2/3
  norms   use Wallis power p' = 2*alpha + 1 = -1/9

Wallis (WITH the phase factor -- pitfall 1):
  I(p,q) = cos(pi q/2) * pi * Gamma(p+1) / (2^p G((p+q+2)/2) G((p-q+2)/2))

Predicted ratio, a = (p+2)/2:
  I(p,2n)/I(p,0) = (1-a)_n / (a)_n
  p=-2/3: (1/3)_n/(2/3)_n     p'=-1/9: (1/18)_n/(17/18)_n

Everything rational is a Pochhammer ratio, so
  b_mk (monic) = nu0 * q_mk,   h_m = mu0 * r_m,   q,r in Q
  normalized entry = gamma * q_mk / sqrt(r_m r_k),
  gamma = nu0/mu0,  nu0 = B(1/2, 1/6),  mu0 = B(1/2, 4/9).

All claims verified against 40-digit theta-quadrature of the DEFINING
integrals (split at pi/4, pi/2, 3pi/4 for the endpoint singularity).
"""

from fractions import Fraction as Fr
import sympy as sp
from mpmath import mp, mpf, quad, gamma as G, pi as PI, cos as mpcos, sin as mpsin, sqrt as mpsqrt

mp.dps = 40
NMAX = 6  # rows 0..NMAX

alpha = sp.Rational(-5, 9)
x = sp.symbols('x')

# ---------------------------------------------------------------- moments
def poch(a, n):
    """Pochhammer (a)_n as exact sympy Rational (a rational)."""
    r = sp.Integer(1)
    for j in range(n):
        r *= (a + j)
    return r

def moment_ratio(expnt, j):
    """int x^(2j) (1-x^2)^expnt dx / int (1-x^2)^expnt dx = (1/2)_j / (expnt+3/2)_j."""
    return poch(sp.Rational(1, 2), j) / poch(expnt + sp.Rational(3, 2), j)

MU = {2 * j: moment_ratio(alpha, j) for j in range(0, 4 * NMAX + 2)}          # W, mu0=1
NU = {2 * j: moment_ratio(sp.Rational(3, 2) * alpha, j) for j in range(0, 4 * NMAX + 2)}  # W^{3/2}, nu0=1

def functional(polycoeffs, M):
    """Apply moment functional M to a poly given as dict deg->Rational."""
    s = sp.Integer(0)
    for d, c in polycoeffs.items():
        if d % 2 == 0:
            s += c * M[d]
    return s

def coeffs(p):
    p = sp.Poly(sp.expand(p), x)
    return {m[0]: c for m, c in zip(p.monoms(), p.coeffs())}

def inner(p, q_, M):
    return functional(coeffs(sp.expand(p * q_)), M)

# ------------------------------------------------- monic OPs (Gram-Schmidt)
P = [sp.Integer(1)]
for n in range(1, NMAX + 1):
    pn = x ** n
    for j in range(n):
        pn -= inner(x ** n, P[j], MU) / inner(P[j], P[j], MU) * P[j]
    P.append(sp.expand(pn))

r = [inner(P[m], P[m], MU) for m in range(NMAX + 1)]           # norms / mu0
q = [[inner(P[m], P[k], NU) for k in range(NMAX + 1)] for m in range(NMAX + 1)]  # entries / nu0

# ------------------------------------------------------------ constants
mu0 = mp.beta(mpf(1) / 2, mpf(4) / 9)      # int (1-x^2)^alpha dx
nu0 = mp.beta(mpf(1) / 2, mpf(1) / 6)      # int (1-x^2)^{3alpha/2} dx
gamma_const = nu0 / mu0

# plan's claimed evaluation of I(-2/3,0) -- check it equals nu0
I_claim = PI * G(mpf(1) / 3) / (2 ** (-mpf(2) / 3) * G(mpf(2) / 3) ** 2)

print("=" * 72)
print("CONSTANT CHECKS")
print("=" * 72)
print(f"nu0 = B(1/2,1/6)              = {nu0}")
print(f"plan claim pi G(1/3) 2^(2/3)/G(2/3)^2 = {I_claim}")
print(f"  agree? diff = {nu0 - I_claim}")
print(f"mu0 = B(1/2,4/9)              = {mu0}")
print(f"gamma = nu0/mu0               = {gamma_const}")

# ------------------------------------------- Wallis-ratio verification
def I_quad(p, q_):
    f = lambda t: mpsin(t) ** p * mpcos(q_ * t)
    return quad(f, [0, PI / 4, PI / 2, 3 * PI / 4, PI])

print()
print("Wallis ratio I(p,2n)/I(p,0) vs (1-a)_n/(a)_n :")
for p_val, a_val, name in [(-mpf(2) / 3, sp.Rational(2, 3), "p=-2/3"),
                           (-mpf(1) / 9, sp.Rational(17, 18), "p=-1/9")]:
    I0 = I_quad(p_val, 0)
    for n in range(1, 4):
        num = I_quad(p_val, 2 * n) / I0
        pred = poch(1 - a_val, n) / poch(a_val, n)
        err = num - mpf(pred.p) / mpf(pred.q)
        print(f"  {name} n={n}: quad={num}  pred={pred}  err={err}")

# ------------------------------------------- entry verification (defining integral)
print()
print("=" * 72)
print("ENTRY VERIFICATION: monic b_mk = nu0*q_mk and h_m = mu0*r_m")
print("=" * 72)
Pnum = [sp.lambdify(x, P[m], 'mpmath') for m in range(NMAX + 1)]

def b_quad(m, k, p_val):
    # substitute t = u^3 near 0 and t = pi - u^3 near pi to kill the
    # sin(t)^p endpoint singularity (integrand becomes smooth: 3u^2 u^{3p})
    # sin(pi - u^3) = sin(u^3) exactly -- avoids cancellation of sin near pi
    left = quad(lambda u: Pnum[m](mpcos(u ** 3)) * Pnum[k](mpcos(u ** 3))
                * mpsin(u ** 3) ** p_val * 3 * u ** 2, [0, (PI / 2) ** (mpf(1) / 3)])
    right = quad(lambda u: Pnum[m](-mpcos(u ** 3)) * Pnum[k](-mpcos(u ** 3))
                 * mpsin(u ** 3) ** p_val * 3 * u ** 2, [0, (PI / 2) ** (mpf(1) / 3)])
    return left + right

for (m, k) in [(0, 0), (1, 1), (2, 0), (2, 2), (3, 1), (4, 2), (4, 4), (5, 5), (6, 6), (6, 0)]:
    lhs = b_quad(m, k, -mpf(2) / 3)
    rhs = nu0 * mpf(q[m][k].p) / mpf(q[m][k].q)
    print(f"  b[{m}][{k}]: quad={lhs}  exact={rhs}  err={lhs - rhs}")
for m in [1, 3, 5, 6]:
    lhs = b_quad(m, m, -mpf(1) / 9)
    rhs = mu0 * mpf(r[m].p) / mpf(r[m].q)
    print(f"  h[{m}]  : quad={lhs}  exact={rhs}  err={lhs - rhs}")

# ------------------------------------------------------------ exact data
print()
print("=" * 72)
print("EXACT DATA (all in Q; entry_mk = gamma * q_mk / sqrt(r_m r_k))")
print("=" * 72)
print("\nmonic norms r_m (units mu0), factored:")
for m in range(NMAX + 1):
    print(f"  r_{m} = {r[m]} = {sp.factorrat(sp.Rational(r[m]))}")

print("\nmonic entries q_mk (units nu0), factored (odd parity should vanish):")
for m in range(NMAX + 1):
    for k in range(m, NMAX + 1):
        v = q[m][k]
        if v == 0:
            continue
        print(f"  q[{m}][{k}] = {v} = {sp.factorrat(sp.Rational(v))}")

print("\nnormalized-entry rationals  t_mk = q_mk/sqrt(r_m r_k)  (entry = gamma*t):")
print("(printed as q_mk^2/(r_m r_k), factored -- square iff t rational)")
for m in range(NMAX + 1):
    for k in range(m, NMAX + 1):
        v = q[m][k]
        if v == 0:
            continue
        t2 = v ** 2 / (r[m] * r[k])
        rt = sp.sqrt(t2)
        tag = "RATIONAL" if rt.is_rational else "irrational"
        print(f"  t[{m}][{k}]^2 = {t2}   sqrt={rt}   [{tag}]")

# ratio structure along rows: does q_{0,2n} follow a sextic factorial?
print("\nrow-0 ratios q[0][2n+2]/q[0][2n]  and  r ratios r_{m+1}/r_m :")
for n in range(0, NMAX // 2):
    print(f"  q0 ratio n={n}: {q[0][2 * n + 2] / q[0][2 * n]}")
for m in range(NMAX):
    print(f"  r ratio  m={m}: {r[m + 1] / r[m]}")
