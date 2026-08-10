"""Chebyshev-T as the lambda -> 0 degeneration of the connection-Gram theorem
(JACOBI-PLAN.md write-up step; companion to jacobi_sextic_connection.py).

Theorem (connection-Gram, verified there for alpha=-5/9): with
lam = alpha+1/2, mu = 3*alpha/2+1/2,

  <pi_m, pi_k>_{W^{3/2}} / I(3alpha+1, 0)
    = mf_m mf_k sum_j A_{(m-j)/2}(m) A_{(k-j)/2}(k) (h_j/h_0),
  A_kk(n) = (lam-mu)_kk (lam)_{n-kk} (n-2kk+mu) / (kk! (mu)_{n-kk+1}),
  h_j/h_0 = (2mu)_j mu / (j! (j+mu)),   mf_n = n!/(2^n (lam)_n).

Claim to verify here: the entries are RATIONAL functions of alpha, and at
alpha = -1/2 (Chebyshev-T, lam = 0, mu = -1/4) they equal the proved Wallis
closed form of jacobi-matching.tex Prop. chebT-closed:

  q^T_00      = 1
  q^T_{0,2n}  = 2^{1-2n} r(n),                      n >= 1
  q^T_mk      = 2^{2-m-k} [r(d) + r(s)] / 2,        m,k >= 1, m+k even
  r(n) = (1/4)_n / (3/4)_n = P(n)/Q(n),  s = (m+k)/2, d = |m-k|/2

(the 2-powers convert monic Chebyshev T~_n = 2^{1-n} T_n, n>=1; everything
in units of the mass I(-1/2,0) of W^{3/2}).

The point: lam = 0 kills C_n^lam but NOT the monic pi_n^lam -- mf_n has the
compensating 1/(lam)_n -- so the limit is evaluation of a rational function.
Equivalently this re-proves the Chebyshev closed form from the general
symmetric-Jacobi mechanism.

Run: python jacobi_lambda0_limit.py   (all-sympy, exact; ~seconds)
"""

import sympy as sp

a = sp.symbols('a')                       # alpha, symbolic
lam = a + sp.Rational(1, 2)
mu = sp.Rational(3, 2) * a + sp.Rational(1, 2)
NMAX = 8

def poch(z, n):
    r = sp.Integer(1)
    for j in range(n):
        r *= (z + j)
    return r

def A(kk, n):
    return (poch(lam - mu, kk) * poch(lam, n - kk) * (n - 2 * kk + mu)
            / (sp.factorial(kk) * poch(mu, n - kk + 1)))

def hratio(j):
    return poch(2 * mu, j) * mu / (sp.factorial(j) * (j + mu))

def mf(n):
    return sp.factorial(n) / (2 ** n * poch(lam, n))

def entry(m, k):
    s = sp.Integer(0)
    for j in range(min(m, k) % 2, min(m, k) + 1, 2):
        s += A((m - j) // 2, m) * A((k - j) // 2, k) * hratio(j)
    return mf(m) * mf(k) * s

# ---- Chebyshev targets (Wallis, proved) --------------------------------
def r(n):
    return poch(sp.Rational(1, 4), n) / poch(sp.Rational(3, 4), n)

def P(n):  # 1*5*...*(4n-3)
    return sp.prod([4 * j - 3 for j in range(1, n + 1)], start=sp.Integer(1))

def Q(n):  # 3*7*...*(4n-1)
    return sp.prod([4 * j - 1 for j in range(1, n + 1)], start=sp.Integer(1))

def cheb_target(m, k):
    m, k = min(m, k), max(m, k)
    if (m + k) % 2 == 1:
        return sp.Integer(0)
    s, d = (m + k) // 2, (k - m) // 2
    if m == 0 and k == 0:
        return sp.Integer(1)
    if m == 0:
        return sp.Integer(2) ** (1 - k) * r(s)
    return sp.Integer(2) ** (2 - m - k) * (r(d) + r(s)) / 2

# ---- verification ------------------------------------------------------
print("lambda->0 (alpha=-1/2) limit of connection-Gram vs Wallis closed form:")
allok = True
for m in range(NMAX + 1):
    for k in range(m, NMAX + 1):
        if (m + k) % 2 == 1:
            continue
        e = sp.cancel(sp.together(entry(m, k)))          # rational in a
        val = sp.cancel(e.subs(a, sp.Rational(-1, 2)))   # evaluate at alpha=-1/2
        tgt = sp.nsimplify(cheb_target(m, k))
        ok = sp.simplify(val - tgt) == 0
        allok = allok and ok
        print(f"  ({m},{k}): limit={val}  target={tgt}  {'OK' if ok else 'MISMATCH'}")
print(f"\nALL ENTRIES MATCH: {allok}")

# ---- integer form D_mk = (1/2)[P(d)Q(s)/Q(d) + P(s)] -------------------
print("\ninteger form check: Q(s) * 2^{m+k-2} * limit == (1/2)[P(d)Q(s)/Q(d)+P(s)], m,k>=1:")
allok2 = True
for m in range(1, NMAX + 1):
    for k in range(m, NMAX + 1):
        if (m + k) % 2 == 1:
            continue
        s, d = (m + k) // 2, (k - m) // 2
        e = sp.cancel(sp.together(entry(m, k))).subs(a, sp.Rational(-1, 2))
        D = sp.cancel(Q(s) * 2 ** (m + k - 2) * e)
        Dtgt = sp.Rational(1, 2) * (P(d) * Q(s) / Q(d) + P(s))
        ok = sp.simplify(D - Dtgt) == 0
        allok2 = allok2 and ok
        print(f"  ({m},{k}): D={D}  target={sp.nsimplify(Dtgt)}  {'OK' if ok else 'MISMATCH'}")
print(f"\nINTEGER FORM MATCHES: {allok2}")

# ---- the degenerate connection coefficients (for the paper's remark) ---
# a_j(n) := lim 2^{n-1} mf_n A_{(n-j)/2}(n) -- the coefficients of
# T_n = sum_j a_j(n) C_j^{-1/4}; check first few against direct expansion.
print("\ndegenerate connection coefficients a_j(n) at mu=-1/4 (T_n in C_j^mu):")
x = sp.symbols('x')
mu0 = sp.Rational(-1, 4)
C = [sp.Integer(1), 2 * mu0 * x]
for n in range(2, NMAX + 1):
    C.append(sp.expand((2 * (n - 1 + mu0) * x * C[n - 1]
                        - (n - 2 + 2 * mu0) * C[n - 2]) / n))
for n in range(1, 5):
    coeffs = {}
    for j in range(n % 2, n + 1, 2):
        aj = sp.cancel((2 ** (n - 1) * mf(n) * A((n - j) // 2, n))).subs(
            a, sp.Rational(-1, 2))
        coeffs[j] = sp.cancel(aj)
    rec = sp.expand(sum(c * C[j] for j, c in coeffs.items()))
    Tn = sp.chebyshevt(n, x)
    ok = sp.expand(rec - Tn) == 0
    print(f"  n={n}: {coeffs}  reconstructs T_n: {'OK' if ok else 'FAIL'}")
