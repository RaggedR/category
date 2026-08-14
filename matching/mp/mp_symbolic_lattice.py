"""Exact symbolic computation of the forward dual-MP matching lattice.

Forward matching integrals (meixner-pollaczek-simple.tex, eq. J):

    J_mk = int_R P_m P_k Gamma(1/4+ix)^2 Gamma(1/4-ix) dx

with master generating function (eq. gf)

    sum J_mk w^m u^k = 2 pi^{3/2} e^{zeta/2} K0(zeta/2) / sqrt((1-iw)(1-iu)),
    zeta = phi(w) phi(u),  phi(t) = (1+it)/(1-it).

KEY FACTORIZATION (this script's foundation; verified below): with
v = 2iw, mu = 2iu, sigma = v/(1-v/2), tau = mu/(1-mu/2),

    zeta = (1+sigma)(1+tau),   i.e.   zeta - 1 = sigma + tau + sigma*tau,
    1/sqrt((1-iw)(1-iu)) = (sigma'(v) tau'(mu))^{1/4}.

We expand g(zeta) = e^{zeta/2} K0(zeta/2) around zeta = 1 through the
closed two-state derivative algebra on {KA, KB} = {e^{z/2}K0(z/2),
e^{z/2}K1(z/2)}:

    d/dz [p KA + q KB] = (p' + p/2 - q/2) KA + (q' - p/2 + q(1/2 - 1/z)) KB,

so every Taylor coefficient g^(n)(1) = p_n KA(1) + q_n KB(1) with
p_n, q_n in Z[1/2].  Combined with the Faa di Bruno expansion of
eps = zeta - 1 (Gaussian-rational Taylor coefficients, denominators
powers of 2 only), every

    T_mk = (2^m m! 2^k k! / h0) J_mk,   h0 = sqrt2 pi^{3/2}

is an EXACT element of Q A + Q B, A = sqrt(2e) K0(1/2), B = sqrt(2e) K1(1/2)
(note e^{1/2}K_nu(1/2) = A_or_B / sqrt2, and 2 pi^{3/2} / h0 = sqrt2).

Checks performed:
  1. reproduce all 17 PSLQ entries of the paper's table (m<=3, k<=5),
  2. produce the missing (3,5) entry,
  3. integrality a_mk, b_mk in Z for all m,k <= N (default 8),
  4. 40+ digit quadrature of the defining integral for a spot set.
"""

from fractions import Fraction
from math import factorial
import mpmath as mp

N = 8            # max index m,k
DPS = 50         # quadrature precision

# ----------------------------------------------------------------------
# Gaussian-rational truncated bivariate series: dict {(i,j): (re, im)}
# with re, im Fraction; truncation deg <= N in each variable.
# ----------------------------------------------------------------------

def gr(re=0, im=0):
    return (Fraction(re), Fraction(im))

def gr_add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def gr_mul(x, y):
    return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])

def gr_scale(x, r):          # r Fraction
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
# Building blocks.  i*w has series coefficient i at (1,0), etc.
# phi(w) = 1 + 2*sum_{l>=1} (iw)^l           (Cayley transform)
# (1-iw)^{-1/2} = sum binom(2n,n)/4^n (iw)^n  (prefactor factor)
# ----------------------------------------------------------------------

def binom(n, k):
    return factorial(n) // (factorial(k) * factorial(n - k))

def i_pow(l):
    """(i)^l as Gaussian rational."""
    return [gr(1, 0), gr(0, 1), gr(-1, 0), gr(0, -1)][l % 4]

def phi_series(var):         # var 0 -> w, 1 -> u
    f = {(0, 0): gr(1)}
    for l in range(1, N + 1):
        key = (l, 0) if var == 0 else (0, l)
        f[key] = gr_scale(i_pow(l), Fraction(2))
    return f

def pref_factor(var):
    f = {}
    for n in range(N + 1):
        key = (n, 0) if var == 0 else (0, n)
        f[key] = gr_scale(i_pow(n), Fraction(binom(2*n, n), 4**n))
    return f

# eps = phi(w) phi(u) - 1
eps = s_mul(phi_series(0), phi_series(1))
eps[(0, 0)] = gr_add(eps[(0, 0)], gr(-1))
pref = s_mul(pref_factor(0), pref_factor(1))

# ----------------------------------------------------------------------
# Taylor coefficients of g(z) = e^{z/2} K0(z/2) at z = 1 in the
# two-state algebra, via sympy rational functions of z.
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
    p, q = (sp.cancel(sp.diff(p, zt) + p/2 - q/2),
            sp.cancel(sp.diff(q, zt) - p/2 + q*(sp.Rational(1, 2) - 1/zt)))

# ----------------------------------------------------------------------
# Assemble  GF/(2 pi^{3/2}) = pref * sum_n g^(n)(1)/n! * eps^n
# split into KA-part and KB-part (KA = e^{1/2}K0(1/2) = A/sqrt2, etc.)
# ----------------------------------------------------------------------

SA, SB = {}, {}
eps_pow = {(0, 0): gr(1)}
for n in range(2*N + 1):
    pn, qn = pq_at_1[n]
    inv_fact = Fraction(1, factorial(n))
    SA = s_add(SA, s_scale(eps_pow, pn * inv_fact))
    SB = s_add(SB, s_scale(eps_pow, qn * inv_fact))
    eps_pow = s_mul(eps_pow, eps)
SA, SB = s_mul(SA, pref), s_mul(SB, pref)

# T_mk = 2^m m! 2^k k! * (a*A + b*B) with (a,b) the (m,k) coefficients
# (the sqrt2 from 2 pi^{3/2}/h0 cancels the 1/sqrt2 in KA = A/sqrt2).

def T_exact(m, k):
    """Return (a, b) Gaussian rationals with T_mk = a*A + b*B."""
    scale = Fraction(2**m * factorial(m) * 2**k * factorial(k))
    a = gr_scale(SA.get((m, k), ZERO), scale)
    b = gr_scale(SB.get((m, k), ZERO), scale)
    return a, b

def strip_phase(m, k, a, b):
    """Table convention: drop the factor i for m+k odd; entries then real."""
    if (m + k) % 2 == 0:
        assert a[1] == 0 and b[1] == 0, (m, k, a, b)
        return a[0], b[0]
    assert a[0] == 0 and b[0] == 0, (m, k, a, b)
    return a[1], b[1]

# ----------------------------------------------------------------------
# Check 1: the paper's PSLQ table (T_mk as coeff_A * A + coeff_B * B).
# ----------------------------------------------------------------------

paper = {
    (0, 0): (1, 0),        (1, 3): (2575, -1440),
    (0, 1): (3, -2),       (1, 4): (39411, -21978),
    (0, 2): (-23, 12),     (1, 5): (-707389, 394932),
    (0, 3): (-241, 138),   (2, 2): (2369, -1320),
    (0, 4): (3265, -1800), (2, 3): (33743, -18834),
    (0, 5): (52907, -29706), (2, 4): (-569719, 317892),
    (1, 1): (-21, 12),     (2, 5): (-11092229, 6190602),
    (1, 2): (-205, 114),   (3, 3): (-529549, 295572),
    (3, 4): (-9675553, 5399634),
}

print("=== Check 1: reproduce PSLQ table ===")
all_ok = True
for (m, k), (pa, pb) in sorted(paper.items()):
    a, b = strip_phase(m, k, *T_exact(m, k))
    # the paper's table is at a fixed phase; allow global sign per entry
    ok = (a, b) == (pa, pb) or (a, b) == (-pa, -pb)
    all_ok &= ok
    print(f"  ({m},{k}): sym = {a}A + {b}B   paper = {pa}A + {pb}B   "
          f"{'OK' if ok else 'MISMATCH'}")
print("Table reproduced." if all_ok else "TABLE MISMATCH — STOP.")

a35, b35 = strip_phase(3, 5, *T_exact(3, 5))
print(f"\n=== The missing entry ===\n  T_35 = {a35}A + {b35}B")

# ----------------------------------------------------------------------
# Check 2: integrality for all m,k <= N.
# ----------------------------------------------------------------------

print(f"\n=== Check 2: integrality m,k <= {N} ===")
worst = 1
integral = True
for m in range(N + 1):
    for k in range(m, N + 1):
        a, b = strip_phase(m, k, *T_exact(m, k))
        for c in (a, b):
            if c.denominator != 1:
                integral = False
                print(f"  NON-INTEGER at ({m},{k}): {c}")
            worst = max(worst, c.denominator)
print("All entries in Z A + Z B." if integral
      else f"INTEGRALITY FAILS (worst denominator {worst}).")

# ----------------------------------------------------------------------
# Check 3: quadrature of the defining integral (spot set), DPS digits.
# ----------------------------------------------------------------------

print(f"\n=== Check 3: quadrature spot checks (dps={DPS}) ===")
mp.mp.dps = DPS

def P_polys(nmax):
    """MP polynomials lambda=1/4, phi=pi/2 as mpmath-coeff lists."""
    polys = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(2)]]
    for n in range(1, nmax):
        # (n+1) P_{n+1} = 2x P_n - (n - 1/2) P_{n-1}
        pn, pm = polys[n], polys[n-1]
        new = [mp.mpf(0)]*(n+2)
        for j, c in enumerate(pn):
            new[j+1] += 2*c
        for j, c in enumerate(pm):
            new[j] -= (n - mp.mpf('0.5'))*c
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
                   * mp.gamma(mp.mpf('0.25') + 1j*x)**2
                   * mp.gamma(mp.mpf('0.25') - 1j*x))
    return mp.quad(f, [-mp.inf, -5, 0, 5, mp.inf])

A = mp.sqrt(2*mp.e) * mp.besselk(0, mp.mpf('0.5'))
B = mp.sqrt(2*mp.e) * mp.besselk(1, mp.mpf('0.5'))
h0 = mp.sqrt(2) * mp.pi**mp.mpf('1.5')

spot = [(0, 0), (0, 1), (1, 1), (3, 5), (0, 7), (2, 6), (7, 7), (8, 8)]
quad_ok = True
for m, k in spot:
    a, b = T_exact(m, k)
    sym = (mp.mpf(a[0].numerator)/a[0].denominator
           + 1j*mp.mpf(a[1].numerator)/a[1].denominator) * A \
        + (mp.mpf(b[0].numerator)/b[0].denominator
           + 1j*mp.mpf(b[1].numerator)/b[1].denominator) * B
    num = 2**m * factorial(m) * 2**k * factorial(k) / h0 * J_quad(m, k)
    err = abs(sym - num) / max(1, abs(num))
    ok = err < mp.mpf(10)**(-(DPS - 12))
    quad_ok &= ok
    print(f"  ({m},{k}): rel err = {mp.nstr(err, 3)}   {'OK' if ok else 'FAIL'}")
print("Quadrature checks pass." if quad_ok else "QUADRATURE FAILURE — STOP.")

# ----------------------------------------------------------------------
# Full exact table (A- and B-integer matrices), for the paper.
# ----------------------------------------------------------------------

print(f"\n=== Exact T_mk = a A + b B, m <= k <= {N} (phase stripped) ===")
for m in range(N + 1):
    for k in range(m, N + 1):
        a, b = strip_phase(m, k, *T_exact(m, k))
        print(f"  ({m},{k}): {a}A + {b}B")

print("\nDiagonal a_nn:", [strip_phase(n, n, *T_exact(n, n))[0]
                           for n in range(N + 1)])
print("Diagonal b_nn:", [strip_phase(n, n, *T_exact(n, n))[1]
                         for n in range(N + 1)])
print("Row 0 a_0k:", [strip_phase(0, k, *T_exact(0, k))[0]
                      for k in range(N + 1)])
print("Row 0 b_0k:", [strip_phase(0, k, *T_exact(0, k))[1]
                      for k in range(N + 1)])

# ----------------------------------------------------------------------
# Check 4: the closed combinatorial formula (the integrality PROOF).
#
# Kummer's equation zeta g'' + (1-zeta) g' - g/2 = 0 for
# g = sqrt(pi) U(1/2,1,zeta) gives, at zeta = 1, the INTEGER recurrence
#
#   u_{n+2} = -2n u_{n+1} + (4n+2) u_n,   u_0 = A,  u_1 = A - B,
#
# for u_n = 2^n sqrt2 g^(n)(1), so u_n in Z A + Z B.  Extracting
# coefficients through the factorization zeta - 1 = sigma + tau +
# sigma*tau (chains = Lah numbers, prefactor = defect double
# factorials via 2^d d! [w^d](1-iw)^{-1/2} = i^d (2d-1)!!) yields
#
#   T_mk = i^{m+k} sum_{d1,c1,d2,c2}
#            C(m,d1)(2d1-1)!! L(m-d1,c1) 2^{m-d1}
#          * C(k,d2)(2d2-1)!! L(k-d2,c2) 2^{k-d2}
#          * sum_j C(c1,j) C(c2,j) j! 2^j  u_{c1+c2-j}
#
# with L(n,c) = binom(n-1,c-1) n!/c! the Lah numbers.  Every weight is
# a positive integer, so T_mk in Z[i] A + Z[i] B, and the phase
# i^{m+k} produces the real/imaginary checkerboard.  We verify the
# formula against the exact expansion for all m,k <= min(N,6).
# ----------------------------------------------------------------------

UREC = [(1, 0), (1, -1)]
for n in range(2*N):
    UREC.append((-2*n*UREC[n+1][0] + (4*n+2)*UREC[n][0],
                 -2*n*UREC[n+1][1] + (4*n+2)*UREC[n][1]))

# consistency: u_n must equal 2^n sqrt2 g^(n)(1) from the derivative algebra
for n in range(2*N + 1):
    pn, qn = pq_at_1[n]
    assert (Fraction(UREC[n][0]), Fraction(UREC[n][1])) == (2**n*pn, 2**n*qn), n

def dfact(n):
    r = 1
    for x in range(1, 2*n, 2):
        r *= x
    return r

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
            wm = binom(m, d1) * dfact(d1) * lah(m - d1, c1) * 2**(m - d1)
            if wm == 0:
                continue
            for d2 in range(k + 1):
                for c2 in range(k - d2 + 1):
                    wk = (binom(k, d2) * dfact(d2)
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
print("u_n =", UREC[:9])
