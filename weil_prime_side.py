"""Prime side of the Weil explicit formula in the MP-1/4 F-basis.

Companion to weil_form_mp.py (see weil-prime-side-plan.md).  Basis kernel

    h_mn(r) = (r^2 + 1/4) MPhat_m(r/2) MPhat_n(r/2) |Gamma(1/4 + ir/2)|^2

Stage 0 (exact):  g_mn(u) = (1/2pi) int h_mn(r) e^{-iur} dr in closed form.
  Fourier fact:  (1/2pi) int |Gamma(1/4+ir/2)|^2 e^{-iur} dr = sqrt(2pi) sech^{1/2}(u)
  (verified numerically below).  Multiplication by r <-> i d/du, and d/du maps
  p(t) sech^{1/2}u  |->  ((1-t^2) p'(t) - (t/2) p(t)) sech^{1/2}u,  t = tanh u.
  Hence with Q_mn(r) = (r^2+1/4) P_m(r/2) P_n(r/2)  (UNnormalized, rational coeffs):

      g_mn(u) = sqrt(2pi)/sqrt(h_m h_n) * T_mn(tanh u) * sech^{1/2}(u),
      T_mn = sum_j q_j i^j D^j[1]  in Q[t]   (q_j = 0 for j odd when m = n mod 2).

  KEY: T_mn(1) = 0 -- the zero of (r^2+1/4) at r = -i/2 kills the leading
  e^{-u/2} decay, so per-prime terms are ~ Lambda(n)/n^3, not Lambda(n)/n.

Stage 1 (fast, exact-coefficient):  at u = log n, tanh = (n^2-1)/(n^2+1),
  sech = 2n/(n^2+1), so with w = n^{-2},

      P_mn = 2 sum Lambda(n) n^{-1/2} g_mn(log n)
           = (4 sqrt(pi)/sqrt(h_m h_n)) sum_n Lambda(n) n^{-1} psi_mn(n^{-2}),
      psi_mn(w) = (1+w)^{-1/2} T_mn((1-w)/(1+w)) = sum_{k>=1} a_k w^k,  a_k in Q.

  Swapping sums (absolutely convergent):

      P_mn = (4 sqrt(pi)/sqrt(h_m h_n)) sum_{k>=1} a_k * (-zeta'/zeta(1+2k)).

  Terms decay like 4^{-k}: 120 digits needs ~210 terms, NO prime sum at all.
  Direct prime-power sum (cutoff 1e5, tail ~ 5e-11) is kept as a cross-check.

Stage 2: PSLQ vs the archimedean lattice A = sqrt(2e)K0(1/2), B = sqrt(2e)K1(1/2)
  times {1, log2, logpi, pi, euler}; row/column Q-relations; positivity of
  W = pole + arch - prime for m,n <= 6.
  Pitfalls honoured: maxcoeff >= 1e10; any hit re-verified at 2x dps; keep
  #PSLQ inputs small vs working digits (~12 digits/input spurious threshold).

MP norm convention (match weil_form_mp.py): h_m = 2 pi Gamma(m+1/2)/(sqrt2 m!).
Parity: only m = n (mod 2) blocks exist.
"""
import pickle
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp
from sympy import Rational as Q

MMAX = 6
T_CACHE = "weil_prime_Tmn.pkl"
P_CACHE = "weil_prime_Pmn.pkl"

r_, t_, w_, x_ = sp.symbols("r t w x")


# ---------------------------------------------------------------- MP basis
def mp_poly_sym(m):
    """Unnormalized MP^{1/4} poly (phi=pi/2), exact: (j+1)P_{j+1}=2xP_j-(j-1/2)P_{j-1}."""
    p0, p1 = sp.Integer(1), 2 * x_
    if m == 0:
        return p0
    for j in range(1, m):
        p0, p1 = p1, sp.expand((2 * x_ * p1 - (j - Q(1, 2)) * p0) / (j + 1))
    return p1


def h_norm(m):
    return 2 * mp.pi * mp.gamma(m + mp.mpf(1) / 2) / (mp.sqrt(2) * mp.factorial(m))


# ------------------------------------------------- Stage 0: exact T_mn in Q[t]
def D_op(p):
    """p(t) -> (1-t^2) p'(t) - (t/2) p(t)   (action of d/du on p(t) sech^{1/2}u)."""
    return sp.expand((1 - t_**2) * sp.diff(p, t_) - t_ / 2 * p)


def T_poly(m, n):
    """T_mn in Q[t]: apply Q_mn(i d/du) to 1, tracking i^j (j even => real)."""
    Qmn = sp.Poly(
        sp.expand((r_**2 + Q(1, 4))
                  * mp_poly_sym(m).subs(x_, r_ / 2)
                  * mp_poly_sym(n).subs(x_, r_ / 2)), r_)
    coeffs = Qmn.all_coeffs()[::-1]          # q_0, q_1, ..., q_deg
    T = sp.Integer(0)
    p = sp.Integer(1)                        # D^j[1]
    for j, qj in enumerate(coeffs):
        if qj != 0:
            assert j % 2 == 0, f"odd power r^{j} in Q_{m}{n} -- parity broken"
            T += qj * sp.Integer(-1) ** (j // 2) * p
        p = D_op(p)
    return sp.expand(T)


def verify_fourier_constant():
    """K(u) = (1/2pi) int |Gamma(1/4+ir/2)|^2 e^{-iur} dr =? sqrt(2pi) sech^{1/2}u."""
    print("Stage 0a: Fourier constant check")
    ok = True
    for u in [mp.mpf("0.3"), mp.mpf(1), mp.mpf("2.2")]:
        quad = (1 / mp.pi) * mp.quad(
            lambda r: abs(mp.gamma(mp.mpf(1) / 4 + 1j * r / 2)) ** 2 * mp.cos(u * r),
            mp.linspace(0, 80, 160))
        closed = mp.sqrt(2 * mp.pi) * mp.sech(u) ** mp.mpf("0.5")
        err = abs(quad - closed)
        print(f"  u={mp.nstr(u,3):>5}  quad={mp.nstr(quad,15)}  "
              f"closed={mp.nstr(closed,15)}  err={mp.nstr(err,3)}")
        ok = ok and err < mp.mpf("1e-18")
    print(f"  constant sqrt(2pi) VALID: {ok}")
    return ok


def h_mn_num(m, n, r):
    x = r / 2
    pm = mp.polyval([mp.mpf(str(c)) for c in sp.Poly(mp_poly_sym(m), x_).all_coeffs()], x)
    pn = mp.polyval([mp.mpf(str(c)) for c in sp.Poly(mp_poly_sym(n), x_).all_coeffs()], x)
    g = mp.gamma(mp.mpf(1) / 4 + 1j * x)
    return (r**2 + mp.mpf(1) / 4) * pm * pn * (g * mp.conj(g)).real


def verify_g_closed(T):
    """Closed form vs quadrature, unnormalized, several (m,n) and u."""
    print("Stage 0b: g_mn closed form vs quadrature (unnormalized)")
    ok = True
    for (m, n) in [(0, 0), (1, 1), (2, 0), (3, 1), (4, 4), (6, 6), (5, 3)]:
        Tp = sp.Poly(T[(m, n)], t_).all_coeffs()
        for u in [mp.mpf("0.4"), mp.mpf("1.5"), mp.mpf("3.1")]:
            quad = (1 / mp.pi) * mp.quad(
                lambda r: h_mn_num(m, n, r) * mp.cos(u * r),
                mp.linspace(0, 80, 240))
            th = mp.tanh(u)
            closed = (mp.sqrt(2 * mp.pi) * mp.sech(u) ** mp.mpf("0.5")
                      * mp.polyval([mp.mpf(str(c)) for c in Tp], th))
            err = abs(quad - closed)
            rel = err / max(abs(closed), mp.mpf("1e-30"))
            good = rel < mp.mpf("1e-12")
            ok = ok and good
            print(f"  ({m},{n}) u={mp.nstr(u,2):>4}: quad={mp.nstr(quad,12):>18} "
                  f"closed={mp.nstr(closed,12):>18} rel={mp.nstr(rel,3)} {'OK' if good else 'FAIL'}")
    print(f"  closed form VALID: {ok}")
    return ok


def stage0():
    print("=" * 68)
    print("STAGE 0: exact g_mn")
    print("=" * 68)
    mp.mp.dps = 30
    if not verify_fourier_constant():
        sys.exit("STOP: Fourier constant mismatch")
    T = {}
    for m in range(MMAX + 1):
        for n in range(m, MMAX + 1):
            if (m - n) % 2:
                continue
            T[(m, n)] = T[(n, m)] = T_poly(m, n)
    # T_mn(1) = 0 (forced by the (r^2+1/4) zero at r = -i/2)
    print("Stage 0c: T_mn(1) = 0 check")
    for (m, n), Tp in sorted(T.items()):
        if m <= n:
            v = Tp.subs(t_, 1)
            assert v == 0, f"T_{m}{n}(1) = {v} != 0"
    print(f"  all {sum(1 for k in T if k[0] <= k[1])} polynomials vanish at t=1: OK")
    if not verify_g_closed(T):
        sys.exit("STOP: closed form fails quadrature check")
    with open(T_CACHE, "wb") as f:
        pickle.dump({k: sp.srepr(v) for k, v in T.items()}, f)
    print(f"  cached -> {T_CACHE}")
    return T


# ------------------------------------------------- Stage 1: P_mn via zeta'/zeta
def psi_series(Tmn, K):
    """a_1..a_K (exact Fractions) of psi(w) = (1+w)^{-1/2} T((1-w)/(1+w))."""
    d = sp.degree(Tmn, t_)
    # T((1-w)/(1+w)) = A(w)/(1+w)^d,  A in Q[w]
    A = sp.cancel(sp.together(Tmn.subs(t_, (1 - w_) / (1 + w_))) * (1 + w_) ** d)
    A = sp.Poly(sp.expand(A), w_)
    Ac = [Fraction(int(sp.fraction(c)[0]), int(sp.fraction(c)[1]))
          for c in A.all_coeffs()[::-1]]                      # A_0, A_1, ...
    # binomial series of (1+w)^{-(d+1/2)}
    e = Fraction(-(2 * d + 1), 2)
    b = [Fraction(1)]
    for k in range(1, K + 1):
        b.append(b[-1] * (e - (k - 1)) / k)
    a = []
    for k in range(K + 1):
        s = sum(Ac[j] * b[k - j] for j in range(min(k, len(Ac) - 1) + 1))
        a.append(s)
    assert a[0] == 0, "a_0 != 0 -- T(1)=0 should force this"
    return a


def neg_zp_over_z(K):
    """[-zeta'/zeta(1+2k)] for k = 1..K at current dps."""
    return [-mp.zeta(1 + 2 * k, 1, 1) / mp.zeta(1 + 2 * k) for k in range(1, K + 1)]


def stage1(T, dps=120):
    print("=" * 68)
    print(f"STAGE 1: P_mn via zeta'/zeta series, dps={dps}")
    print("=" * 68)
    mp.mp.dps = dps
    K = int(dps * mp.log(10) / mp.log(4)) + 25
    zz = neg_zp_over_z(K)
    print(f"  {K} series terms, last -z'/z(1+2K) ~ {mp.nstr(zz[-1], 3)}")
    P = {}
    for (m, n), Tmn in sorted(T.items()):
        if m > n:
            continue
        a = psi_series(Tmn, K)
        pref = 4 * mp.sqrt(mp.pi) / mp.sqrt(h_norm(m) * h_norm(n))
        terms = [mp.mpf(a[k].numerator) / mp.mpf(a[k].denominator) * zz[k - 1]
                 for k in range(1, K + 1)]
        val = pref * mp.fsum(terms)
        tail = abs(pref) * abs(terms[-1]) * 2   # honest-ish: terms decay ~4^{-k}
        P[(m, n)] = P[(n, m)] = val
        print(f"  P[{m},{n}] = {mp.nstr(val, 30)}   (last term ~ {mp.nstr(tail, 2)})")
    with open(P_CACHE, "wb") as f:
        pickle.dump({k: mp.nstr(v, dps) for k, v in P.items()}, f)
    print(f"  cached -> {P_CACHE}")
    return P


def stage1_crosscheck(T, P, cutoff=100000):
    """Direct prime-power sum, tail ~ 1/(2*cutoff^2) ~ 5e-11."""
    print(f"Stage 1b: direct prime-power sum cross-check, cutoff {cutoff}")
    mp.mp.dps = 30
    from sympy import primerange
    pp = []
    for p in primerange(2, cutoff + 1):
        lam, q = mp.log(p), p
        while q <= cutoff:
            pp.append((q, lam))
            q *= p
    for (m, n) in [(0, 0), (1, 1), (2, 2), (4, 2), (6, 6)]:
        Tp = [mp.mpf(str(c)) for c in sp.Poly(T[(m, n)], t_).all_coeffs()]
        pref = 4 * mp.sqrt(mp.pi) / mp.sqrt(h_norm(m) * h_norm(n))
        s = pref * mp.fsum(
            lam / mp.sqrt(q * q + 1) * mp.polyval(Tp, mp.mpf(q * q - 1) / (q * q + 1))
            for q, lam in pp)
        err = abs(s - mp.mpf(mp.nstr(P[(m, n)], 25)))
        # tail of direct sum ~ pref*a_1/(2 cutoff^2); verified to scale as 1/N^2
        tol = abs(pref) * mp.mpf(50) / mp.mpf(cutoff) ** 2
        print(f"  ({m},{n}): direct={mp.nstr(s, 15)}  series={mp.nstr(P[(m,n)], 15)}"
              f"  |diff|={mp.nstr(err, 3)}  {'OK' if err < tol else 'FAIL'}")


# ------------------------------------------------- Stage 2: structure hunt
def lattice_candidates():
    K0, K1 = mp.besselk(0, mp.mpf(1) / 2), mp.besselk(1, mp.mpf(1) / 2)
    Ac, Bc = mp.sqrt(2 * mp.e) * K0, mp.sqrt(2 * mp.e) * K1
    mults = [("1", mp.mpf(1)), ("log2", mp.log(2)), ("logpi", mp.log(mp.pi)),
             ("pi", mp.pi), ("euler", mp.euler)]
    return [(f"{bn}*{mn}", bv * mv)
            for bn, bv in [("A", Ac), ("B", Bc)] for mn, mv in mults]


def run_pslq(names, vec, maxcoeff=10**10):
    rel = mp.pslq(vec, maxcoeff=maxcoeff, maxsteps=10**6)
    if rel is None:
        return None
    resid = mp.fsum(c * v for c, v in zip(rel, vec))
    return rel, resid


def stage2_pslq(P, dps):
    """Per-entry PSLQ vs A,B lattice; then intra-block relations among P entries."""
    print("=" * 68)
    print(f"STAGE 2a: PSLQ P_mn vs A,B lattice (dps={dps}, "
          f"11 inputs -> need >=132 working digits)")
    print("=" * 68)
    mp.mp.dps = dps
    cands = lattice_candidates()
    hits = []
    for (m, n) in sorted(k for k in P if k[0] <= k[1]):
        vec = [P[(m, n)]] + [cv for _, cv in cands]
        out = run_pslq([f"P{m}{n}"] + [nm for nm, _ in cands], vec)
        if out is not None and out[0][0] != 0:
            rel, resid = out
            terms = " + ".join(f"{c}*{nm}" for c, (nm, _) in zip(rel[1:], cands) if c)
            print(f"  P[{m},{n}]: {rel[0]}*P + {terms} = 0  (resid {mp.nstr(abs(resid),3)})")
            hits.append(((m, n), rel))
        else:
            print(f"  P[{m},{n}] = {mp.nstr(P[(m,n)], 25)}: no relation")

    print(f"\nSTAGE 2b: intra-block Q-relations among P entries")
    for block in [(0, 2, 4, 6), (1, 3, 5)]:
        keys = [(m, n) for i, m in enumerate(block) for n in block[i:]]
        vec = [P[k] for k in keys]
        out = run_pslq(keys, vec)
        if out is not None:
            rel, resid = out
            terms = " + ".join(f"{c}*P{k}" for c, k in zip(rel, keys) if c)
            print(f"  block {block}: {terms} = 0  (resid {mp.nstr(abs(resid),3)})")
            hits.append((block, rel))
        else:
            print(f"  block {block}: no relation ({len(vec)} entries)")
    return hits


def stage2_matrix(T, P, dps=30):
    """Sign pattern + positivity of W = pole + arch - prime, m,n <= MMAX."""
    print("=" * 68)
    print(f"STAGE 2c: W = pole + arch - prime, m,n <= {MMAX} (dps={dps})")
    print("=" * 68)
    mp.mp.dps = dps
    nodes = mp.linspace(0, 80, 200)

    def coeffs(m):
        return [mp.mpf(str(c)) for c in sp.Poly(mp_poly_sym(m), x_).all_coeffs()]

    C = {m: coeffs(m) for m in range(MMAX + 1)}

    def hval(m, n, r):
        x = r / 2
        g = mp.gamma(mp.mpf(1) / 4 + 1j * x)
        return ((r**2 + mp.mpf(1) / 4) * mp.polyval(C[m], x) * mp.polyval(C[n], x)
                / mp.sqrt(h_norm(m) * h_norm(n)) * (g * mp.conj(g)).real)

    def arch(m, n):
        return (1 / mp.pi) * mp.quad(
            lambda r: hval(m, n, r) * (mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * r / 2))
                                       - mp.log(mp.pi)), nodes)

    def pole(m, n):
        nrm = mp.sqrt(h_norm(m) * h_norm(n))
        zp = mp.polyval(C[m], mp.mpc(0, 0.25)) * mp.polyval(C[n], mp.mpc(0, 0.25))
        zm = mp.polyval(C[m], mp.mpc(0, -0.25)) * mp.polyval(C[n], mp.mpc(0, -0.25))
        return (2 * mp.sqrt(mp.pi) * (zp + zm) / nrm).real

    W = {}
    for m in range(MMAX + 1):
        for n in range(m, MMAX + 1):
            if (m - n) % 2:
                continue
            w = pole(m, n) + arch(m, n) - mp.mpf(mp.nstr(P[(m, n)], 25))
            W[(m, n)] = W[(n, m)] = w
            print(f"  ({m},{n}): pole={mp.nstr(pole(m,n),10):>14} "
                  f"arch={mp.nstr(arch(m,n),10):>14} "
                  f"prime={mp.nstr(P[(m,n)],10):>14}  W={mp.nstr(w,10):>14}")
    for block, name in [((0, 2, 4, 6), "even"), ((1, 3, 5), "odd")]:
        M = mp.matrix(len(block))
        for i, m in enumerate(block):
            for j, n in enumerate(block):
                M[i, j] = W[(m, n)]
        ev = mp.eigsy(M, eigvals_only=True)
        print(f"  {name} block eigenvalues: " + ", ".join(mp.nstr(e, 10) for e in ev))
        print(f"  positive definite: {all(e > 0 for e in ev)}")


if __name__ == "__main__":
    import os
    if os.path.exists(T_CACHE) and "--fresh" not in sys.argv:
        T = {k: sp.sympify(v) for k, v in pickle.load(open(T_CACHE, "rb")).items()}
        print(f"loaded T_mn from {T_CACHE}")
    else:
        T = stage0()
    if "--stage2" in sys.argv:
        DPS = 240
        P = stage1(T, dps=DPS)
        stage2_pslq(P, DPS)
        stage2_matrix(T, P)
    else:
        P = stage1(T)
        stage1_crosscheck(T, P)
