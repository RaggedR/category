"""Weil explicit-formula quadratic form in the Meixner--Pollaczek basis.

RH <=> Weil positivity: Q(f) = sum_rho h(gamma_rho) >= 0 for all admissible
h = |fhat|^2.  We expand fhat on the critical line in the MP-1/4 dual basis
F_m = Gamma(1/4+ix) MPhat_m (x = r/2), polarize, and regularize:

    h_mn(r) = (r^2 + 1/4) * MPhat_m(r/2) MPhat_n(r/2) * |Gamma(1/4 + ir/2)|^2

The MP weight reappears as the Weil kernel.  The (r^2+1/4) factor cancels
the pole of Gamma(s/2) at s = 0,1 (x = -+ i/4) -- the archimedean factor's
pole sits exactly at the zeta pole.  Positive factor on R => still a
positivity family.  CAUTION (bug found in first run): the pole terms
h(+-i/2) are NOT killed -- zero x Gamma-pole leaves the finite residue
    h(+-i/2) = 2 Gamma(1/2) MPhat_m(+-i/4) MPhat_n(+-i/4),
an evaluation at the branch point x = i/4 (= s = 0, 1).

Explicit formula convention used (to be VALIDATED in stage 1 against a sum
over actual zeta zeros):  for h even, analytic in |Im r| <= 1/2+eps,
g(u) = (1/2pi) int h(r) e^{-iur} dr:

  sum_gamma h(gamma) = h(i/2) + h(-i/2) - g(0) log(pi)
                       + (1/2pi) int_R h(r) Re[psi(1/4 + ir/2)] dr
                       - 2 sum_{n>=2} Lambda(n) n^{-1/2} g(log n)

Parity note: MPhat_m has parity (-1)^m, so h_mn is even iff m = n mod 2;
zeros come in +-gamma pairs, so only the even part of h pairs with the
formula -> the form is block-diagonal in parity.  We compute both blocks.

Stages:
  1. Validate the formula with h(r) = exp(-(r/6)^2) against the first 60
     nontrivial zeros (all terms numerically visible at this width).
  2. Compute W_mn = A_mn - P_mn (arch minus primes) for the even block
     (m,n in 0,2,4) and odd block (1,3), indices m <= 4.  Check positivity.
  3. Recompute A_mn at dps=40 and PSLQ against the MP lattice constants
     A = sqrt(2e)K0(1/2), B = sqrt(2e)K1(1/2) times {1, log 2, pi, gamma}.
"""
import mpmath as mp
from sympy import primerange

mp.mp.dps = 30

MMAX = 4
PRIME_CUTOFF = 300
R_NODES = mp.linspace(0, 40, 120)   # oscillatory g-integrals, u <= log(300)


# ---------------------------------------------------------------- MP basis
def mp_poly(m, x):
    """MP^{1/4} polynomial, phi=pi/2: (j+1)P_{j+1} = 2xP_j - (j-1/2)P_{j-1}."""
    p0, p1 = mp.mpf(1), 2*x
    if m == 0:
        return p0
    for j in range(1, m):
        p0, p1 = p1, (2*x*p1 - (j - mp.mpf(1)/2)*p0) / (j + 1)
    return p1


def h_norm(m):
    """h_m = 2 pi Gamma(m+1/2) / (sqrt2 m!)."""
    return 2*mp.pi*mp.gamma(m + mp.mpf(1)/2) / (mp.sqrt(2)*mp.factorial(m))


def Wmp(x):
    g = mp.gamma(mp.mpf(1)/4 + 1j*x)
    return (g*mp.conj(g)).real


def h_mn(m, n, r):
    x = r/2
    return ((r**2 + mp.mpf(1)/4) * mp_poly(m, x) * mp_poly(n, x)
            / mp.sqrt(h_norm(m)*h_norm(n)) * Wmp(x))


# ------------------------------------------------- explicit-formula pieces
def arch_term(h):
    """(1/2pi) int_R h(r) [Re psi(1/4+ir/2) - log pi] dr  (h even)."""
    kern = lambda r: h(r) * (mp.re(mp.digamma(mp.mpf(1)/4 + 1j*r/2))
                             - mp.log(mp.pi))
    return (1/mp.pi) * mp.quad(kern, R_NODES)     # 2 * (1/2pi) * int_0^inf


def g_of(h, u):
    """g(u) = (1/2pi) int h(r) e^{-iur} dr = (1/pi) int_0^inf h cos(ur) dr."""
    return (1/mp.pi) * mp.quad(lambda r: h(r)*mp.cos(u*r), R_NODES)


def von_mangoldt_terms(cutoff):
    """[(n, Lambda(n))] for prime powers n <= cutoff."""
    out = []
    for p in primerange(2, cutoff + 1):
        lam = mp.log(p)
        q = p
        while q <= cutoff:
            out.append((q, lam))
            q *= p
    return sorted(out)


PP = von_mangoldt_terms(PRIME_CUTOFF)


def prime_term(h):
    """2 sum Lambda(n) n^{-1/2} g(log n)."""
    return 2*mp.fsum(lam / mp.sqrt(q) * g_of(h, mp.log(q)) for q, lam in PP)


def arithmetic_side(h, h_at_i2):
    return h_at_i2 * 2 + arch_term(h) - prime_term(h)


# ------------------------------------------------------------ stage 1
def validate():
    print("=" * 68)
    print("STAGE 1: validate explicit formula, h(r) = exp(-(r/6)^2)")
    print("=" * 68)
    h = lambda r: mp.exp(-(r/6)**2)
    h_i2 = mp.exp(-(mp.mpc(0, 0.5)/6)**2).real    # h(i/2) = e^{+1/144}
    rhs = arithmetic_side(h, h_i2)

    zero_sum = mp.mpf(0)
    NZ = 60
    for j in range(1, NZ + 1):
        gam = mp.im(mp.zetazero(j))
        zero_sum += 2*h(gam)          # zeros come in pairs +-gamma
    tail = 2*h(mp.im(mp.zetazero(NZ)))
    print(f"  sum over first {NZ} zero pairs : {mp.nstr(zero_sum, 15)}")
    print(f"  arithmetic side              : {mp.nstr(rhs, 15)}")
    print(f"  difference                   : {mp.nstr(zero_sum - rhs, 6)}"
          f"   (last-zero term ~ {mp.nstr(tail, 3)})")
    ok = abs(zero_sum - rhs) < max(tail*5, mp.mpf('1e-8'))
    print(f"  CONVENTION VALID: {ok}")
    return ok


# ------------------------------------------------------------ stage 2
def pole_term(m, n):
    """h(i/2) + h(-i/2) = 2 Gamma(1/2) [P^_mP^_n(i/4) + P^_mP^_n(-i/4)]."""
    nrm = mp.sqrt(h_norm(m)*h_norm(n))
    z = mp_poly(m, mp.mpc(0, 0.25)) * mp_poly(n, mp.mpc(0, 0.25)) / nrm
    zb = mp_poly(m, mp.mpc(0, -0.25)) * mp_poly(n, mp.mpc(0, -0.25)) / nrm
    return (2*mp.sqrt(mp.pi)*(z + zb)).real


def weil_matrix():
    print("\n" + "=" * 68)
    print(f"STAGE 2: W_mn in MP basis, m,n <= {MMAX}, prime cutoff {PRIME_CUTOFF}")
    print("=" * 68)
    A, P, Wq = {}, {}, {}
    for m in range(MMAX + 1):
        for n in range(m, MMAX + 1):
            if (m - n) % 2:
                continue                       # parity block-diagonal
            h = lambda r, m=m, n=n: h_mn(m, n, r)
            a = arch_term(h)
            p = prime_term(h)
            pole = pole_term(m, n)
            A[(m, n)] = A[(n, m)] = a
            P[(m, n)] = P[(n, m)] = p
            Wq[(m, n)] = Wq[(n, m)] = pole + a - p
            print(f"  ({m},{n}):  pole = {mp.nstr(pole, 10):>15s}   "
                  f"arch = {mp.nstr(a, 10):>15s}   "
                  f"prime = {mp.nstr(p, 10):>15s}   "
                  f"W = {mp.nstr(pole + a - p, 10):>15s}", flush=True)

    for block, name in [((0, 2, 4), "even"), ((1, 3), "odd")]:
        M = mp.matrix(len(block))
        for i, m in enumerate(block):
            for j, n in enumerate(block):
                M[i, j] = Wq[(m, n)]
        ev = mp.eigsy(M, eigvals_only=True)
        print(f"\n  {name} block eigenvalues: "
              + ", ".join(mp.nstr(e, 10) for e in ev))
        print(f"  positive definite: {all(e > 0 for e in ev)}")
    return A


# ------------------------------------------------------------ stage 3
def pslq_arch(A):
    print("\n" + "=" * 68)
    print("STAGE 3: PSLQ on archimedean entries (dps=40 recompute)")
    print("=" * 68)
    mp.mp.dps = 40
    K0, K1 = mp.besselk(0, mp.mpf(1)/2), mp.besselk(1, mp.mpf(1)/2)
    Ac = mp.sqrt(2*mp.e)*K0
    Bc = mp.sqrt(2*mp.e)*K1
    base = [('A', Ac), ('B', Bc)]
    mults = [('1', mp.mpf(1)), ('log2', mp.log(2)), ('pi', mp.pi),
             ('euler', mp.euler), ('logpi', mp.log(mp.pi))]
    cands = [(f"{bn}*{mn}", bv*mv) for bn, bv in base for mn, mv in mults]
    for (m, n) in sorted(k for k in A if k[0] <= k[1]):
        h = lambda r, m=m, n=n: h_mn(m, n, r)
        a = arch_term(h)
        vec = [a] + [cv for _, cv in cands]
        rel = mp.pslq(vec, maxcoeff=10**10, maxsteps=50000)
        if rel is not None and rel[0] != 0:
            resid = mp.fsum(r*v for r, v in zip(rel, vec))
            names = " + ".join(f"{c}*{nm}" for c, (nm, _) in zip(rel[1:], cands)
                               if c != 0)
            print(f"  A[{m},{n}] = {mp.nstr(a, 25)}\n"
                  f"      rel: {rel[0]}*v + {names} = 0   "
                  f"(resid {mp.nstr(abs(resid), 3)})")
        else:
            print(f"  A[{m},{n}] = {mp.nstr(a, 25)}: no relation")


if __name__ == "__main__":
    if validate():
        A = weil_matrix()
        pslq_arch(A)
    else:
        print("\nSTOP: convention mismatch -- fix before trusting stages 2-3.")
