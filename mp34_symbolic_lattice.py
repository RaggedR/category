"""Exact symbolic proof of the MP lambda=3/4 forward matching lattice.

Forward matching integrals (mp34-matching.tex):

    J_mk = int_R P_m P_k Gamma(3/4+ix)^2 Gamma(3/4-ix) dx

with master generating function (mp34-matching.tex Thm gf)

    sum J_mk w^m u^k = (pi^2/2) U(3/2,1,zeta) / ((1-iw)^{3/2}(1-iu)^{3/2}),
    zeta = phi(w) phi(u),  phi(t) = (1+it)/(1-it).

Same Cayley factorization as lambda=1/4 (mp_symbolic_lattice.py):
zeta = (1+sigma)(1+tau), sigma = v/(1-v/2), v = 2iw.  The ONLY
family-dependent data:

  (1) prefactor exponent 3/2:  2^d d! [w^d](1-iw)^{-3/2} = i^d (2d+1)!!
  (2) kernel y(z) = U(3/2,1,z):  z y'' + (1-z) y' - (3/2) y = 0
      => at z=1, u'_n := 2^n y^(n)(1) obeys the INTEGER recurrence

          u'_{n+2} = -2n u'_{n+1} + (4n+6) u'_n,
          u'_0 = A' = U(3/2,1,1),  u'_1 = 2y'(1) = -3B',
          B' = U(5/2,2,1)  (since U'(a,b,z) = -a U(a+1,b+1,z)).

Everything else (Lah chains, Charlier matching sum) is verbatim from
the lambda=1/4 proof — this is claim 2 of the dictionary ("the
combinatorics is universal; the family only chooses u_n").

Checks:
  1. reproduce the 16 published table entries of mp34-matching.tex,
  2. integrality for all m,k <= N (default 8),
  3. quadrature of the defining integral (dps=50) for a spot set,
  4. the closed combinatorial formula vs the exact expansion (the proof).
"""

from fractions import Fraction
from math import factorial
import mpmath as mp

N = 8            # max index m,k
DPS = 50         # quadrature precision

# ----------------------------------------------------------------------
# Gaussian-rational truncated bivariate series: dict {(i,j): (re, im)}
# ----------------------------------------------------------------------

def gr(re=0, im=0):
    return (Fraction(re), Fraction(im))

def gr_add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def gr_mul(x, y):
    return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])

def gr_scale(x, r):
    return (x[0]*r, x[1]*r)

ZERO = gr()

def s_add(f, g):
    h = dict(f)
    for key, c in g.items():
        h[key] = gr_add(h.get(key, ZERO), c)
    return h

def s_mul(f, g):
    h = {}
    for (i1, j1), c1 in f.items():
        for (i2, j2), c2 in g.items():
            i, j = i1 + i2, j1 + j2
            if i > N or j > N:
                continue
            key = (i, j)
            h[key] = gr_add(h.get(key, ZERO), gr_mul(c1, c2))
    return h

def s_scale(f, r):
    return {k: gr_scale(c, r) for k, c in f.items()}

# ----------------------------------------------------------------------
# Building blocks.
# phi(w) = 1 + 2*sum_{l>=1} (iw)^l
# (1-iw)^{-3/2} = sum_n (3/2)_n/n! (iw)^n,  (3/2)_n = (2n+1)!!/2^n
# ----------------------------------------------------------------------

def binom(n, k):
    return factorial(n) // (factorial(k) * factorial(n - k))

def dfact_odd(n):
    """(2n+1)!! = 1*3*5*...*(2n+1)."""
    r = 1
    for x in range(1, 2*n + 2, 2):
        r *= x
    return r

def i_pow(l):
    return [gr(1, 0), gr(0, 1), gr(-1, 0), gr(0, -1)][l % 4]

def phi_series(var):
    f = {(0, 0): gr(1)}
    for l in range(1, N + 1):
        key = (l, 0) if var == 0 else (0, l)
        f[key] = gr_scale(i_pow(l), Fraction(2))
    return f

def pref_factor(var):
    f = {}
    for n in range(N + 1):
        key = (n, 0) if var == 0 else (0, n)
        f[key] = gr_scale(i_pow(n),
                          Fraction(dfact_odd(n), 2**n * factorial(n)))
    return f

# eps = phi(w) phi(u) - 1
eps = s_mul(phi_series(0), phi_series(1))
eps[(0, 0)] = gr_add(eps[(0, 0)], gr(-1))
pref = s_mul(pref_factor(0), pref_factor(1))

# ----------------------------------------------------------------------
# Taylor coefficients of y(z) = U(3/2,1,z) at z=1 in the (y, y')
# state algebra: y'' = ((z-1) y' + (3/2) y)/z.
# state (p,q) means p*y + q*y'; derivative rule below.
# ----------------------------------------------------------------------
import sympy as sp

zt = sp.symbols('z')
p, q = sp.Integer(1), sp.Integer(0)
pq_at_1 = []                 # (p_n(1), q_n(1)) as Fractions
for n in range(2*N + 1):
    p1 = sp.Rational(sp.together(p).subs(zt, 1))
    q1 = sp.Rational(sp.together(q).subs(zt, 1))
    pq_at_1.append((Fraction(int(p1.p), int(p1.q)),
                    Fraction(int(q1.p), int(q1.q))))
    p, q = (sp.cancel(sp.diff(p, zt) + q*sp.Rational(3, 2)/zt),
            sp.cancel(p + sp.diff(q, zt) + q*(zt - 1)/zt))

# ----------------------------------------------------------------------
# Assemble  GF/(pi^2/2) = pref * sum_n y^(n)(1)/n! * eps^n
# split into A'-part (coefficient of y(1)) and B'-part.
# y(1) = A',  y'(1) = -(3/2) B'.
# ----------------------------------------------------------------------

SA, SB = {}, {}
eps_pow = {(0, 0): gr(1)}
for n in range(2*N + 1):
    pn, qn = pq_at_1[n]
    inv_fact = Fraction(1, factorial(n))
    SA = s_add(SA, s_scale(eps_pow, pn * inv_fact))
    SB = s_add(SB, s_scale(eps_pow, qn * Fraction(-3, 2) * inv_fact))
    eps_pow = s_mul(eps_pow, eps)
SA, SB = s_mul(SA, pref), s_mul(SB, pref)

def T_exact(m, k):
    """Return (a, b) Gaussian rationals with U_mk = a*A' + b*B'."""
    scale = Fraction(2**m * factorial(m) * 2**k * factorial(k))
    a = gr_scale(SA.get((m, k), ZERO), scale)
    b = gr_scale(SB.get((m, k), ZERO), scale)
    return a, b

def strip_phase(m, k, a, b):
    if (m + k) % 2 == 0:
        assert a[1] == 0 and b[1] == 0, (m, k, a, b)
        return a[0], b[0]
    assert a[0] == 0 and b[0] == 0, (m, k, a, b)
    return a[1], b[1]

# ----------------------------------------------------------------------
# Check 1: the published table (mp34-matching.tex Table 1).
# ----------------------------------------------------------------------

paper = {
    (0, 0): (1, 0),           (1, 3): (7035, -11760),
    (0, 1): (3, -6),          (1, 4): (128643, -213630),
    (0, 2): (39, -60),        (1, 5): (2669529, -4440900),
    (0, 3): (513, -870),      (2, 2): (6609, -10920),
    (0, 4): (8673, -14280),   (2, 3): (112695, -187350),
    (0, 5): (165579, -276270),(2, 4): (2215047, -3679860),
    (1, 1): (33, -60),        (3, 3): (2078169, -3456420),
    (1, 2): (453, -750),      (3, 4): (43522881, -72346830),
}

print("=== Check 1: reproduce published table ===")
all_ok = True
for (m, k), (pa, pb) in sorted(paper.items()):
    a, b = strip_phase(m, k, *T_exact(m, k))
    ok = (a, b) == (pa, pb) or (a, b) == (-pa, -pb)
    all_ok &= ok
    print(f"  ({m},{k}): sym = {a}A' + {b}B'   paper = {pa}A' + {pb}B'   "
          f"{'OK' if ok else 'MISMATCH'}")
print("Table reproduced." if all_ok else "TABLE MISMATCH — STOP.")

# ----------------------------------------------------------------------
# Check 2: integrality for all m,k <= N.
# ----------------------------------------------------------------------

print(f"\n=== Check 2: integrality m,k <= {N} ===")
integral = True
for m in range(N + 1):
    for k in range(m, N + 1):
        a, b = strip_phase(m, k, *T_exact(m, k))
        for c in (a, b):
            if c.denominator != 1:
                integral = False
                print(f"  NON-INTEGER at ({m},{k}): {c}")
print("All entries in Z A' + Z B'." if integral else "INTEGRALITY FAILS.")

# ----------------------------------------------------------------------
# Check 3: quadrature of the defining integral (spot set).
# ----------------------------------------------------------------------

print(f"\n=== Check 3: quadrature spot checks (dps={DPS}) ===")
mp.mp.dps = DPS

def P_polys(nmax):
    """MP polynomials lambda=3/4, phi=pi/2:
    (n+1)P_{n+1} = 2x P_n - (n+1/2) P_{n-1}."""
    polys = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(2)]]
    for n in range(1, nmax):
        pn, pm = polys[n], polys[n-1]
        new = [mp.mpf(0)]*(n+2)
        for j, c in enumerate(pn):
            new[j+1] += 2*c
        for j, c in enumerate(pm):
            new[j] -= (n + mp.mpf('0.5'))*c
        polys.append([c/(n+1) for c in new])
    return polys

POL = P_polys(N)

def evalp(coeffs, x):
    v = mp.mpf(0)
    for c in reversed(coeffs):
        v = v*x + c
    return v

def J_quad(m, k):
    f = lambda x: (evalp(POL[m], x) * evalp(POL[k], x)
                   * mp.gamma(mp.mpf('0.75') + 1j*x)**2
                   * mp.gamma(mp.mpf('0.75') - 1j*x))
    return mp.quad(f, [-mp.inf, -5, 0, 5, mp.inf])

Ap = mp.hyperu(mp.mpf('1.5'), 1, 1)
Bp = mp.hyperu(mp.mpf('2.5'), 2, 1)
scale0 = mp.pi**2 / 2

spot = [(0, 0), (0, 1), (1, 1), (3, 5), (0, 7), (2, 6), (7, 7), (8, 8)]
quad_ok = True
for m, k in spot:
    a, b = T_exact(m, k)
    sym = (mp.mpf(a[0].numerator)/a[0].denominator
           + 1j*mp.mpf(a[1].numerator)/a[1].denominator) * Ap \
        + (mp.mpf(b[0].numerator)/b[0].denominator
           + 1j*mp.mpf(b[1].numerator)/b[1].denominator) * Bp
    num = 2**m * factorial(m) * 2**k * factorial(k) / scale0 * J_quad(m, k)
    err = abs(sym - num) / max(1, abs(num))
    ok = err < mp.mpf(10)**(-(DPS - 12))
    quad_ok &= ok
    print(f"  ({m},{k}): rel err = {mp.nstr(err, 3)}   {'OK' if ok else 'FAIL'}")
print("Quadrature checks pass." if quad_ok else "QUADRATURE FAILURE — STOP.")

# ----------------------------------------------------------------------
# Check 4: the closed combinatorial formula (the integrality PROOF).
#
#   U_mk = i^{m+k} sum_{d1,c1,d2,c2}
#            C(m,d1)(2d1+1)!! L(m-d1,c1) 2^{m-d1}
#          * C(k,d2)(2d2+1)!! L(k-d2,c2) 2^{k-d2}
#          * sum_j C(c1,j) C(c2,j) j! 2^j  u'_{c1+c2-j}
#
#   u'_{n+2} = -2n u'_{n+1} + (4n+6) u'_n,  u'_0 = A', u'_1 = -3B'.
# ----------------------------------------------------------------------

UREC = [(1, 0), (0, -3)]
for n in range(2*N):
    UREC.append((-2*n*UREC[n+1][0] + (4*n+6)*UREC[n][0],
                 -2*n*UREC[n+1][1] + (4*n+6)*UREC[n][1]))

# consistency: u'_n must equal 2^n y^(n)(1) from the derivative algebra
for n in range(2*N + 1):
    pn, qn = pq_at_1[n]
    assert (Fraction(UREC[n][0]), Fraction(UREC[n][1])) == \
        (2**n * pn, 2**n * qn * Fraction(-3, 2)), n

def lah(n, c):
    if n == c == 0:
        return 1
    if not (1 <= c <= n):
        return 0
    return binom(n - 1, c - 1) * factorial(n) // factorial(c)

def T_formula(m, k):
    sA = sB = 0
    for d1 in range(m + 1):
        for c1 in range(m - d1 + 1):
            wm = binom(m, d1) * dfact_odd(d1) * lah(m - d1, c1) * 2**(m - d1)
            if wm == 0:
                continue
            for d2 in range(k + 1):
                for c2 in range(k - d2 + 1):
                    wk = (binom(k, d2) * dfact_odd(d2)
                          * lah(k - d2, c2) * 2**(k - d2))
                    if wk == 0:
                        continue
                    for j in range(min(c1, c2) + 1):
                        wj = binom(c1, j) * binom(c2, j) * factorial(j) * 2**j
                        a, b = UREC[c1 + c2 - j]
                        sA += wm * wk * wj * a
                        sB += wm * wk * wj * b
    return sA, sB

print("\n=== Check 4: combinatorial formula (the proof) ===")
M4 = min(N, 6)
form_ok = True
for m in range(M4 + 1):
    for k in range(m, M4 + 1):
        a, b = T_exact(m, k)
        fa, fb = T_formula(m, k)
        ip = [(1, 0), (0, 1), (-1, 0), (0, -1)][(m + k) % 4]
        pred_a = (Fraction(fa*ip[0]), Fraction(fa*ip[1]))
        pred_b = (Fraction(fb*ip[0]), Fraction(fb*ip[1]))
        if (a, b) != (pred_a, pred_b):
            form_ok = False
            print(f"  MISMATCH at ({m},{k})")
print(f"Formula = exact table for all m,k <= {M4}:"
      f" {'PASS' if form_ok else 'FAIL'}")
print("u'_n =", UREC[:9])

print("\n=== Exact U_mk = a A' + b B', m <= k <= {} (phase stripped) ===".format(N))
for m in range(N + 1):
    for k in range(m, N + 1):
        a, b = strip_phase(m, k, *T_exact(m, k))
        print(f"  ({m},{k}): {a}A' + {b}B'")

print("\nDiagonal a_nn:", [strip_phase(n, n, *T_exact(n, n))[0]
                           for n in range(N + 1)])
print("Diagonal b_nn:", [strip_phase(n, n, *T_exact(n, n))[1]
                         for n in range(N + 1)])
