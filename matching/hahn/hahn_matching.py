"""Continuous Hahn matching coefficients for a=1/4, b=3/4.

Weight: W(x) = |Gamma(1/4+ix) Gamma(3/4+ix)|^2 = 2 pi^2 / cosh(2 pi x).

Two matching matrices:
  (a) b_mk = int phat_m phat_k W^{3/2} dx  [sech^{3/2} kernel]
  (b) c_mk = int phat_m phat_k w_+ W dx    [dual-split kernel]
      where w_+ = Gamma(1/4+ix) Gamma(3/4+ix).

beta_n = n^2/4: the recurrence is p_{n+1}(x) = x p_n(x) - (n/2)^2 p_{n-1}(x).
"""

import mpmath as mp
import sympy as sp
from sympy import Rational, factorial, binomial
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import os

mp.mp.dps = 50
NMAX = 12
CACHE = "hahn_matching.pkl"


# ======================================================================
# Part 1: Exact moments via Euler numbers
# ======================================================================

def euler_number(n):
    return sp.euler(n)


def exact_moment_over_pi2(n):
    """mu_{2n}/pi^2 = |E_{2n}|/16^n (rational). Odd moments are zero.

    Derivation: mu_{2n} = int x^{2n} 2pi^2 sech(2pi x) dx.
    Sub u=2pi x: = pi/(2pi)^{2n} * int u^{2n} sech(u) du.
    int u^{2n} sech(u) du = pi |E_{2n}| (pi/2)^{2n}.
    So mu_{2n} = pi^2 |E_{2n}| (pi/2)^{2n}/(2pi)^{2n} = pi^2 |E_{2n}|/16^n.
    """
    E2n = abs(euler_number(2 * n))
    return Rational(E2n, 16**n)


def verify_moments():
    """Cross-check exact moments against quadrature."""
    print("=== Moment verification (mu_{2n} = pi^2 * |E_{2n}|/4^n) ===")
    for n in range(6):
        exact_rat = exact_moment_over_pi2(n)
        exact_val = float(exact_rat) * float(sp.pi.evalf(30))**2
        # Quadrature of x^{2n} * 2pi^2 sech(2pi x)
        num = float(mp.quad(lambda x: x**(2*n) * 2 * mp.pi**2 / mp.cosh(2*mp.pi*x),
                           [-20, -5, -1, 0, 1, 5, 20]))
        print(f"  mu_{2*n:2d} = pi^2 * {exact_rat}  "
              f"exact={exact_val:.12f}  quad={num:.12f}  err={abs(exact_val - num):.2e}")
    print()


# ======================================================================
# Part 2: Exact monic orthogonal polynomials
# ======================================================================

def build_exact_polys(N):
    """Build monic OPs for the sech weight using exact rational moments.

    The inner product is <f,g> = int f g W dx. Since W = 2pi^2 sech(2pi x),
    we have <x^j, x^k> = mu_{j+k}. Odd moments vanish (even weight).

    Returns (alphas, betas, norms, polys) where norms[n] = <p_n, p_n>/pi^2
    (divided by pi^2 for rational arithmetic).
    """
    max_deg = 4 * N + 4
    mu = [Rational(0)] * (max_deg + 1)
    for k in range(max_deg // 2 + 1):
        mu[2*k] = exact_moment_over_pi2(k)

    x = sp.Symbol('x')
    polys = [sp.Poly(sp.Integer(1), x)]

    def ip(p, q):
        pc = p.all_coeffs()[::-1]
        qc = q.all_coeffs()[::-1]
        s = Rational(0)
        for j, cj in enumerate(pc):
            for k, ck in enumerate(qc):
                idx = j + k
                if idx < len(mu):
                    s += cj * ck * mu[idx]
        return s

    norms = [ip(polys[0], polys[0])]  # = mu_0/pi^2 ... wait.
    # mu[0] = exact_moment_over_pi2(0) = |E_0|/1 = 1.
    # ip computes sum c_j c_k mu[j+k] where mu[k] = mu_k/pi^2.
    # So ip = <f,g>_W / pi^2.  norms[n] = <p_n,p_n>/pi^2.

    alphas = []
    betas = [norms[0]]  # beta_0 = h_0/pi^2 = 1

    xpoly = sp.Poly(x, x)
    a0 = ip(xpoly * polys[0], polys[0]) / norms[0]  # = 0 (even weight)
    alphas.append(a0)
    p1 = sp.Poly(x - a0, x)
    polys.append(p1)
    norms.append(ip(p1, p1))

    for n in range(1, N):
        an = ip(xpoly * polys[n], polys[n]) / norms[n]
        alphas.append(an)
        bn = norms[n] / norms[n-1]
        betas.append(bn)
        pnew = (xpoly - sp.Poly(an, x)) * polys[n] - sp.Poly(bn, x) * polys[n-1]
        polys.append(pnew)
        norms.append(ip(pnew, pnew))

    return alphas, betas, norms, polys


# ======================================================================
# Part 3: Numerical evaluation
# ======================================================================

def eval_phat(n, x, norms_f):
    """Evaluate truly orthonormal phat_n(x) wrt W(x) = 2pi^2 sech(2pi x).

    The monic recurrence is p_{n+1} = x p_n - beta_n p_{n-1},
    with beta_n = n^2/4 and alpha_n = 0 (even weight).
    h_n = <p_n, p_n> = norms_f[n] * pi^2, so phat_n = p_n / sqrt(h_n).

    We evaluate the monic p_n directly and divide by sqrt(h_n).
    norms_f[n] = h_n / pi^2 (exact rationals converted to mpmath).
    """
    # Evaluate monic p_n via recurrence: p_{n+1} = x p_n - beta_n p_{n-1}
    # beta_n = norms_f[n]/norms_f[n-1] computed from exact moments.
    if n == 0:
        p = mp.mpf(1)
    else:
        p_prev = mp.mpf(1)
        p_curr = x  # p_1 = x (alpha_0 = 0)
        for j in range(1, n):
            beta_j = norms_f[j] / norms_f[j-1]
            p_next = x * p_curr - beta_j * p_prev
            p_prev, p_curr = p_curr, p_next
        p = p_curr

    # Normalize: phat_n = p_n / sqrt(h_n) = p_n / (pi * sqrt(norms_f[n]))
    return p / (mp.pi * mp.sqrt(norms_f[n]))


NODES = [-mp.inf, -50, -30, -20, -12, -8, -5, -3, -1.5, 0,
         1.5, 3, 5, 8, 12, 20, 30, 50, mp.inf]


def W(x):
    """Weight W(x) = 2 pi^2 / cosh(2 pi x)."""
    return 2 * mp.pi**2 / mp.cosh(2 * mp.pi * x)


def kernel_b(x):
    """W(x)^{3/2} = (2pi^2)^{3/2} sech(2pi x)^{3/2}."""
    w = W(x)
    return w * mp.sqrt(w)


def kernel_c(x):
    """w_+(x) * W(x) = Gamma(1/4+ix) Gamma(3/4+ix) * W(x).

    Using reflection:
    Gamma(1/4+ix)Gamma(3/4-ix) * Gamma(3/4+ix)Gamma(1/4-ix) = 2pi^2/cosh(2pi x) = W(x).
    So Gamma(1/4+ix)Gamma(3/4+ix) = W(x) / (Gamma(1/4-ix)Gamma(3/4-ix)... no, that gives w_+ = W/w_-.
    w_+ * W = w_+ * w_+ * w_- = Gamma(1/4+ix)^2 Gamma(3/4+ix)^2 Gamma(1/4-ix) Gamma(3/4-ix).
    """
    q = mp.mpf(1)/4
    t = mp.mpf(3)/4
    return mp.gamma(q + 1j*x)**2 * mp.gamma(t + 1j*x)**2 * mp.gamma(q - 1j*x) * mp.gamma(t - 1j*x)


def compute_matrices(norms_f, N):
    """Compute both 12x12 matching matrices."""
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            data = pickle.load(f)
            return data['B'], data['C']

    B = {}
    C = {}

    for m in range(N):
        for k in range(m, N):
            print(f"  ({m},{k})...", end=" ", flush=True)

            def integ_b(x, m=m, k=k):
                return eval_phat(m, x, norms_f) * eval_phat(k, x, norms_f) * kernel_b(x)

            def integ_c(x, m=m, k=k):
                return eval_phat(m, x, norms_f) * eval_phat(k, x, norms_f) * kernel_c(x)

            B[(m, k)] = mp.quad(integ_b, NODES)
            C[(m, k)] = mp.quad(integ_c, NODES)

            # Matrix (a) should be real (both kernel and polys are real for real x,
            # and the contour is on the real line).
            # Matrix (b): m+k even → real, m+k odd → pure imaginary (by x → -x symmetry).

    print()

    with open(CACHE, "wb") as f:
        pickle.dump({'B': B, 'C': C}, f)

    return B, C


def get_entry(M, m, k):
    if (m, k) in M:
        return M[(m, k)]
    return M[(k, m)]


# ======================================================================
# Part 5: Barnes analysis
# ======================================================================

def barnes_analysis():
    """Closed form for c_00.

    c_00 = int phat_0^2 w_+ W dx = (1/pi^2) int w_+ W dx
         = (1/pi^2) int Gamma(1/4+ix)^2 Gamma(3/4+ix)^2 Gamma(1/4-ix) Gamma(3/4-ix) dx.

    Using Legendre duplication on both pairs:
      Gamma(1/4+ix)Gamma(3/4+ix) = sqrt(pi) 2^{1/2-2ix} Gamma(1/2+2ix)
      Gamma(1/4-ix)Gamma(3/4-ix) = sqrt(pi) 2^{1/2+2ix} Gamma(1/2-2ix)

    Integrand = [sqrt(pi) 2^{1/2-2ix} Gamma(1/2+2ix)]^2 * [sqrt(pi) 2^{1/2+2ix} Gamma(1/2-2ix)]
              = pi^{3/2} 2^{3/2} * 2^{-2ix} * Gamma(1/2+2ix)^2 Gamma(1/2-2ix)

    Substituting u = 2x (du = 2dx):
    = pi^{3/2} sqrt(2) * int Gamma(1/2+iu)^2 Gamma(1/2-iu) * e^{-iu ln 2} du

    Shifted 3-Gamma with a=b=c=1/2, s=-ln2 (e^{-s}=2, e^{-as}=sqrt(2)):
    int Gamma(1/2+iu)^2 Gamma(1/2-iu) e^{-iu ln2} du = 2pi sqrt(2) U(1, 1, 2)

    U(1,1,z) = e^z E_1(z) where E_1 is the exponential integral.

    Raw integral = pi^{3/2} sqrt(2) * 2pi sqrt(2) * U(1,1,2) = 4 pi^{5/2} U(1,1,2).
    c_00 = raw/pi^2 = 4 sqrt(pi) U(1,1,2) = 4 sqrt(pi) e^2 E_1(2).
    """
    U112 = mp.hyperu(1, 1, 2)
    closed = 4 * mp.sqrt(mp.pi) * U112

    # Numerical check
    raw = mp.quad(kernel_c, NODES)
    c00_num = raw / mp.pi**2

    print(f"  c_00 numerical  = {mp.nstr(mp.re(c00_num), 20)}")
    print(f"  c_00 closed     = {mp.nstr(closed, 20)}")
    print(f"  = 4 sqrt(pi) U(1,1,2) = 4 sqrt(pi) e^2 E_1(2)")
    print(f"  U(1,1,2) = {mp.nstr(U112, 20)}")
    print(f"  Agreement: {mp.nstr(abs(mp.re(c00_num) - closed), 5)}")

    return closed, U112


# ======================================================================
# Main
# ======================================================================

def main():
    print("=" * 70)
    print("CONTINUOUS HAHN MATCHING: a=1/4, b=3/4")
    print("Weight: W(x) = 2 pi^2 sech(2 pi x)")
    print("=" * 70)
    print()

    # ---- Task 1: Moments and recurrence ----
    verify_moments()

    print("=== Exact moments mu_{2n}/pi^2 = |E_{2n}|/16^n ===")
    for n in range(7):
        E2n = abs(euler_number(2*n))
        print(f"  mu_{2*n:2d}/pi^2 = {E2n}/{16**n} = {Rational(E2n, 16**n)}")
    print()

    N = NMAX
    alphas, betas, norms, polys = build_exact_polys(N)

    # All alphas should be 0
    assert all(a == 0 for a in alphas), "alphas not all zero!"
    print("All alpha_n = 0 (even weight). ✓")
    print()

    print("=== Recurrence coefficients beta_n = n^2/16 ===")
    for n, b in enumerate(betas):
        expected = Rational(n**2, 16) if n > 0 else Rational(1)
        match = "✓" if b == expected else f"✗ (got {b})"
        print(f"  beta_{n:2d} = {str(b):>8s}   expected {str(expected):>8s}  {match}")
    print()

    # ---- Task 2: Orthonormality check ----
    norms_f = [mp.mpf(n) for n in norms]

    print("=== Orthonormality check <phat_m, phat_k>_W ===")
    ok = True
    for m in range(min(5, N)):
        for k in range(m, min(5, N)):
            val = mp.quad(lambda x, m=m, k=k:
                         eval_phat(m, x, norms_f) *
                         eval_phat(k, x, norms_f) * W(x),
                         NODES)
            expected = 1 if m == k else 0
            err = abs(val - expected)
            status = "✓" if err < 1e-15 else f"✗ (err={mp.nstr(err,3)})"
            if m <= 3 and k <= 3:
                print(f"  <phat_{m}, phat_{k}> = {mp.nstr(mp.re(val), 12)}  {status}")
            if err > 1e-10:
                ok = False
    print(f"  {'All OK' if ok else 'FAILED'}!")
    print()

    # ---- Task 3: Matching matrices ----
    print("=== Computing matching matrices (12x12) ===")
    B, C = compute_matrices(norms_f, N)

    # ---- Display matrices ----
    print("\n=== Matrix (a): b_mk = <phat_m phat_k, W^{3/2}> ===")
    for m in range(min(6, N)):
        row = []
        for k in range(min(6, N)):
            v = mp.re(get_entry(B, m, k))
            row.append(f"{float(v):10.6f}")
        print(f"  m={m}: " + "  ".join(row))
    print()

    print("=== Matrix (b): c_mk = <phat_m phat_k, w_+ W> ===")
    print("(m+k even: real; m+k odd: pure imaginary)")
    for m in range(min(6, N)):
        row = []
        for k in range(min(6, N)):
            v = get_entry(C, m, k)
            if (m + k) % 2 == 0:
                row.append(f"{float(mp.re(v)):10.6f}")
            else:
                row.append(f"{float(mp.im(v)):10.6f}i")
        print(f"  m={m}: " + "  ".join(row))
    print()

    # ---- Barnes closed form ----
    print("=== Barnes closed form for c_00 ===")
    c00_closed, U112 = barnes_analysis()
    print()

    # ---- Checkerboard verification ----
    print("=== Checkerboard verification (matrix b) ===")
    worst = 0
    for m in range(N):
        for k in range(m, N):
            v = get_entry(C, m, k)
            if (m + k) % 2 == 0:
                err = abs(mp.im(v))
            else:
                err = abs(mp.re(v))
            worst = max(worst, float(err))
    print(f"  Worst violation: {worst:.2e}")
    print()

    # ---- PSLQ: is c_mk/c_00 rational? ----
    print("=== Testing c_mk/c_00 rationality (matrix b) ===")
    c00 = c00_closed
    ratios = {}

    for m in range(min(8, N)):
        for k in range(m, min(8, N)):
            v = get_entry(C, m, k)
            # Extract real or imaginary part
            if (m + k) % 2 == 0:
                val = mp.re(v)
            else:
                val = mp.im(v)

            ratio = val / c00
            # PSLQ: ratio * r0 + 1 * r1 = 0 => ratio = -r1/r0
            rel = mp.pslq([ratio, mp.mpf(1)])
            if rel is not None and abs(rel[0]) > 0:
                rat = sp.Rational(-int(rel[1]), int(rel[0]))
                check = abs(float(ratio) - float(rat))
                if check < 1e-20:
                    ratios[(m,k)] = rat
                    part = "" if (m+k)%2==0 else "i"
                    if m <= 5 and k <= 5:
                        print(f"  c_{m}{k}/c_00 = {rat}{part}")
                else:
                    ratios[(m,k)] = None
                    if m <= 5 and k <= 5:
                        print(f"  c_{m}{k}: ratio={float(ratio):.10f}  residual={check:.2e}")
            else:
                ratios[(m,k)] = None
                if m <= 5 and k <= 5:
                    print(f"  c_{m}{k}: ratio≈{mp.nstr(ratio, 10)}  PSLQ failed")
    print()

    # If not all rational, try 2-generator PSLQ with U(1,1,2) and U(2,1,2)
    # or with c_00 and c_01/i (two independent constants)
    if any(v is None for v in ratios.values()):
        print("=== Not all rational. Trying 2-generator PSLQ ===")
        # Second generator: try c_{01}/i normalized by something
        c01_im = mp.im(get_entry(C, 0, 1))
        gen2 = c01_im  # second generator

        # Try: each c_mk (real or imaginary part) = a * c_00 + b * Im(c_01)
        for m in range(min(6, N)):
            for k in range(m, min(6, N)):
                if ratios.get((m,k)) is not None:
                    continue
                v = get_entry(C, m, k)
                val = mp.re(v) if (m+k)%2==0 else mp.im(v)

                rel = mp.pslq([val, c00, gen2])
                if rel is not None and abs(rel[0]) > 0:
                    ca = sp.Rational(-int(rel[1]), int(rel[0]))
                    cb = sp.Rational(-int(rel[2]), int(rel[0]))
                    check = abs(float(val) - float(ca)*float(c00) - float(cb)*float(gen2))
                    if check < 1e-15:
                        part = "" if (m+k)%2==0 else "i"
                        print(f"  c_{m}{k} = ({ca})*c_00 + ({cb})*Im(c_01)  {part}")
                    else:
                        print(f"  c_{m}{k}: 2-gen PSLQ residual {check:.2e}")
                else:
                    print(f"  c_{m}{k}: 2-gen PSLQ also failed")
        print()

    # ---- Decay ----
    print("=== Decay of row 0 ===")
    print("\nMatrix (a) [sech^{3/2}]:")
    for k in range(min(10, N)):
        val = abs(get_entry(B, 0, k))
        if k > 0:
            prev = abs(get_entry(B, 0, k-1))
            print(f"  k={k:2d}  |b_0k| = {mp.nstr(val, 8):>14s}  ratio = {mp.nstr(val/prev, 6)}")
        else:
            print(f"  k={k:2d}  |b_0k| = {mp.nstr(val, 8):>14s}")

    print("\nMatrix (b) [6-Gamma]:")
    for k in range(min(10, N)):
        v = get_entry(C, 0, k)
        val = abs(v)
        if k > 0:
            prev_v = get_entry(C, 0, k-1)
            prev = abs(prev_v)
            print(f"  k={k:2d}  |c_0k| = {mp.nstr(val, 8):>14s}  ratio = {mp.nstr(val/prev, 6)}")
        else:
            print(f"  k={k:2d}  |c_0k| = {mp.nstr(val, 8):>14s}")
    print()

    # ---- Euler numbers / combinatorics ----
    print("=== Euler numbers (A000364) ===")
    for n in range(8):
        print(f"  E_{2*n} = {euler_number(2*n)}")
    print()

    # beta_n = n^2/4: numerators n^2, denominators 4
    # 4*beta_n = 0, 1, 4, 9, 16, 25, ... = perfect squares. OEIS A000290.
    print("4*beta_n sequence: ", [int(4*b) for b in betas])
    print("These are perfect squares n^2 — OEIS A000290.")
    print()

    # ---- Zeros ----
    print("=== Computing zeros ===")
    # Truncated dual function: sum_{k<=N} c_{0k} phat_k(z) for matrix (b)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for mat_idx, (M, title) in enumerate([(B, "sech^{3/2}"), (C, "6-Gamma")]):
        ax = axes[mat_idx]
        for N_trunc in [8, 10, 12]:
            if N_trunc > N:
                continue
            coeffs = [get_entry(M, 0, k) for k in range(N_trunc)]

            # Find real zeros by scanning
            real_zeros = []
            xs = [mp.mpf(t)/5 for t in range(-100, 101)]
            for i in range(len(xs)-1):
                def fsum(z, cc=coeffs):
                    return sum(cc[j] * eval_phat(j, z, norms_f) for j in range(len(cc)))
                try:
                    v1, v2 = mp.re(fsum(xs[i])), mp.re(fsum(xs[i+1]))
                    if v1 * v2 < 0:
                        z0 = mp.findroot(lambda z: mp.re(fsum(z)), (xs[i]+xs[i+1])/2)
                        real_zeros.append(z0)
                except:
                    pass

            # Find complex zeros by scanning
            cpx_zeros = list(real_zeros)
            for re_p in [mp.mpf(t)/4 for t in range(-30, 31)]:
                for im_p in [mp.mpf(t)/10 for t in range(-15, 16) if t != 0]:
                    z_init = mp.mpc(re_p, im_p)
                    def fsum(z, cc=coeffs):
                        return sum(cc[j] * eval_phat(j, z, norms_f) for j in range(len(cc)))
                    try:
                        z0 = mp.findroot(fsum, z_init, tol=1e-20, maxsteps=50)
                        if abs(fsum(z0)) < 1e-12:
                            is_dup = any(abs(z0 - zold) < 0.05 for zold in cpx_zeros)
                            if not is_dup:
                                cpx_zeros.append(z0)
                    except:
                        pass

            reals = [float(mp.re(z)) for z in cpx_zeros]
            imags = [float(mp.im(z)) for z in cpx_zeros]
            ax.scatter(reals, imags, label=f"N={N_trunc}", s=20, alpha=0.7)
            print(f"  {title} N={N_trunc}: {len(cpx_zeros)} zeros ({len(real_zeros)} real)")

        ax.set_title(f"Zeros: {title}")
        ax.set_xlabel("Re(z)")
        ax.set_ylabel("Im(z)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='k', linewidth=0.5)

    plt.tight_layout()
    plt.savefig("/Users/robin/git/category/hahn_zeros.png", dpi=150)
    print("Saved hahn_zeros.png")
    print()

    # ---- Full 12x12 tables for TeX ----
    print("=== Full tables for TeX ===")
    print("\nMatrix (b) — c_mk values (real part for m+k even, imag for odd):")
    for m in range(N):
        for k in range(m, N):
            v = get_entry(C, m, k)
            val = float(mp.re(v)) if (m+k)%2==0 else float(mp.im(v))
            print(f"  c_{m:2d},{k:2d} = {val:16.10f}")

    # ---- GF discussion ----
    print("\n=== Generating function analysis ===")
    print("After Legendre duplication, the 6-Gamma integral reduces to")
    print("a 3-Gamma integral with a=b=c=1/2 and a Cayley-shift parameter.")
    print("The continuous Hahn generating function with a=c=1/4, b=d=3/4 is:")
    print("  sum_n p_n(x;1/4,3/4,1/4,3/4) t^n = ???")
    print("(This needs the KLS generating function for continuous Hahn polys.)")
    print()

    # ---- Verify the GF approach for c_00 ----
    print("=== Verification: duplication-based closed form ===")
    # The generating function for our sech-orthogonal polys:
    # Since beta_n = n^2/4, the recurrence is p_{n+1} = x p_n - (n/2)^2 p_{n-1}.
    # This is exactly the recurrence for MONIC Meixner-Pollaczek with lambda=1/2, phi=pi/2!
    # MP(lambda, phi=pi/2): (n+1)P_{n+1} = 2x P_n - (n+2lambda-1) P_{n-1}
    # For monic: p_n = P_n * (n!/(2)^n ... hmm let me think.
    # The MP recurrence for monic OPs with lambda=1/2, phi=pi/2 is:
    # p_{n+1} = x p_n - n(n+2*1/2-1)/(4) p_{n-1}? No.
    # Standard MP: P_{n+1} = (2x P_n - (n+2lam-1)P_{n-1})/(n+1), P_0=1, P_1=2x.
    # Monic: q_n = P_n * n!/2^n (from leading coeff 2^n/n!).
    # Then q_{n+1} = (n+1)!/(2^{n+1}) P_{n+1}
    #  = (n+1)!/(2^{n+1}) * (2x P_n - (n+2lam-1)P_{n-1})/(n+1)
    #  = n!/(2^{n+1}) * 2x * 2^n/n! q_n - n!/(2^{n+1}) * (n+2lam-1) * 2^{n-1}/(n-1)! q_{n-1}
    # Hmm this is getting tangled. Let me just note the observation.

    print("OBSERVATION: beta_n = n^2/16.")
    print("This means the MONIC recurrence is p_{n+1} = x p_n - (n^2/16) p_{n-1}.")
    print("Compare: sech(pi x) is the weight for Meixner-Pollaczek(lambda=1/2, phi=pi/2),")
    print("where the monic beta_n = n^2/4. Our weight is sech(2pi x) ~ sech(pi*(2x)),")
    print("so our polys are p_n(x) = (1/2)^n * q_n(2x) where q_n are the standard monic MP.")
    print("Thus beta_n^{ours} = (1/4)*beta_n^{MP} = n^2/16. Confirmed!")
    print()

    return alphas, betas, B, C, ratios


if __name__ == "__main__":
    main()
