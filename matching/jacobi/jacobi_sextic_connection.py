"""Connection-coefficient closed form for the sextic matching matrix
(JACOBI-PLAN.md Task 1, part 3 -- the structural theorem candidate).

Claim: with lambda = alpha + 1/2 (polys) and mu = 3*alpha/2 + 1/2 (weight
W^{3/2} = (1-x^2)^{mu - 1/2}), the Gegenbauer connection formula

  C_n^lam = sum_k A_k(n) C_{n-2k}^mu,
  A_k(n) = (lam-mu)_k (lam)_{n-k} (n-2k+mu) / (k! (mu)_{n-k+1})

turns each forward matching entry into a finite Pochhammer sum:

  <C_m^lam, C_k^lam>_{W^{3/2}} = sum_j A_{(m-j)/2}(m) A_{(k-j)/2}(k) h_j^mu

with h_j^mu / h_0^mu = the mu-Gegenbauer norm ratio
  h_j^mu = pi 2^{1-2mu} G(j+2mu) / (j! (j+mu) G(mu)^2)  =>
  h_j^mu/h_0^mu = (2mu)_j mu / (j! (j+mu)).

Verification: compare against the exact Gram-Schmidt values q_mk from
jacobi_sextic.py (monic frame), for all m,k <= 8, alpha = -5/9.
Then extract closed form for row 0 and the diagonal.
"""

import sympy as sp

x = sp.symbols('x')
alpha = sp.Rational(-5, 9)
lam = alpha + sp.Rational(1, 2)          # -1/18
mu = sp.Rational(3, 2) * alpha + sp.Rational(1, 2)   # -1/3
NMAX = 8

def poch(a, n):
    r = sp.Integer(1)
    for j in range(n):
        r *= (a + j)
    return r

# ---- independent exact values (Gram-Schmidt, as in jacobi_sextic.py) ----
def moment_ratio(expnt, j):
    return poch(sp.Rational(1, 2), j) / poch(expnt + sp.Rational(3, 2), j)

MU = {2 * j: moment_ratio(alpha, j) for j in range(0, 4 * NMAX + 2)}
NU = {2 * j: moment_ratio(sp.Rational(3, 2) * alpha, j) for j in range(0, 4 * NMAX + 2)}

def functional(p, M):
    p = sp.Poly(sp.expand(p), x)
    return sum(c * M[m[0]] for m, c in zip(p.monoms(), p.coeffs()) if m[0] % 2 == 0)

P = [sp.Integer(1)]
for n in range(1, NMAX + 1):
    pn = x ** n
    for j in range(n):
        pn -= functional(x ** n * P[j], MU) / functional(P[j] * P[j], MU) * P[j]
    P.append(sp.expand(pn))
q = [[functional(P[m] * P[k], NU) for k in range(NMAX + 1)] for m in range(NMAX + 1)]

# ---- connection-coefficient prediction --------------------------------
def A(kk, n, lam_, mu_):
    return (poch(lam_ - mu_, kk) * poch(lam_, n - kk) * (n - 2 * kk + mu_)
            / (sp.factorial(kk) * poch(mu_, n - kk + 1)))

def hratio(j, mu_):
    """h_j^mu / h_0^mu = (2mu)_j * mu / (j! (j+mu))."""
    return poch(2 * mu_, j) * mu_ / (sp.factorial(j) * (j + mu_))

def monic_factor(n, lam_):
    """pi_n = n!/(2^n (lam)_n) C_n^lam."""
    return sp.factorial(n) / (2 ** n * poch(lam_, n))

def entry_pred(m, k):
    """<pi_m, pi_k>_{W^{3/2}} / nu0 predicted by connection formula."""
    s = sp.Integer(0)
    for j in range(min(m, k) % 2, min(m, k) + 1, 2):
        s += A((m - j) // 2, m, lam, mu) * A((k - j) // 2, k, lam, mu) * hratio(j, mu)
    return monic_factor(m, lam) * monic_factor(k, lam) * s

print(f"lambda = {lam}, mu = {mu}, lambda-mu = {lam - mu}")
print("connection-formula prediction vs exact Gram-Schmidt q_mk:")
allok = True
for m in range(NMAX + 1):
    for k in range(m, NMAX + 1):
        if (m + k) % 2 == 1:
            continue
        pred = sp.together(entry_pred(m, k))
        ok = sp.simplify(pred - q[m][k]) == 0
        allok = allok and ok
        print(f"  ({m},{k}): {'OK' if ok else f'MISMATCH pred={pred} exact={q[m][k]}'}")
print(f"\nALL ENTRIES MATCH: {allok}")

# ---- closed forms -----------------------------------------------------
n_ = sp.symbols('n', positive=True, integer=True)
print("\nrow 0 closed form: q_{0,2n} = (2n)!/(2^{2n}(lam)_{2n}) *"
      " (lam-mu)_n (lam)_n mu / (n! (mu)_{n+1})")
for nn in range(0, 5):
    cf = (sp.factorial(2 * nn) / (2 ** (2 * nn) * poch(lam, 2 * nn))
          * poch(lam - mu, nn) * poch(lam, nn) * mu / (sp.factorial(nn) * poch(mu, nn + 1)))
    print(f"  n={nn}: closed={sp.cancel(cf)}  exact={q[0][2 * nn]}  "
          f"{'OK' if sp.simplify(cf - q[0][2 * nn]) == 0 else 'MISMATCH'}")

# generic alpha, symbolic: keep alpha symbolic and print the general row-0 term
a_s = sp.symbols('alpha')
lam_s = a_s + sp.Rational(1, 2)
mu_s = sp.Rational(3, 2) * a_s + sp.Rational(1, 2)
print(f"\ngeneral symmetric Jacobi: lam-mu = {sp.simplify(lam_s - mu_s)}  "
      f"(the '-alpha/2' Pochhammer -- vanishes iff alpha=0, Legendre)")
