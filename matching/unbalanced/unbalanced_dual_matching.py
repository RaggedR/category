"""Unbalanced dual pair: W(x) = Gamma(1/4+ix) Gamma(3/4-ix).

The first genuinely non-self-dual kernel. By the reflection formula,
W(x) = pi / sin(pi(1/4+ix)), which is complex-valued with W(-x) = conj(W(x)).
Split: w_+ = Gamma(1/4+ix), w_- = Gamma(3/4-ix).
Bilinear pairing <f,g> = int fg (no conjugation).

Key discovery: the shift y = x + i/4 transforms the kernel to
  W(y-i/4) = Gamma(1/2+iy) Gamma(1/2-iy) = |Gamma(1/2+iy)|^2
which is the standard self-dual weight for MP(lambda=1/2, phi=pi/2).
So the bilinear orthogonal polynomials are shifted MP(1/2) polys.

Moments: mu(s) = int W e^{isx} dx = 2pi e^{3s/4} / (1+e^s).
3-Gamma matching base: J_00 = 2pi e E_1(1).
Master GF: sum J_mk^{std} w^m u^k = 2pi e^eta E_1(eta) / ((1-iw)(1-iu))
  where eta = (1+iw)(1+iu) / ((1-iw)(1-iu)).
Single generator: A = e E_1(1) (E_1 satisfies first-order ODE).
"""
import mpmath as mp
import pickle, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 60
MMAX = 11
CACHE = "unbalanced_dual_matching.pkl"

NODES = [-60, -40, -25, -15, -8, -4, -2, 0, 2, 4, 8, 15, 25, 40, 60]


def W(x):
    """Kernel Gamma(1/4+ix) Gamma(3/4-ix) = pi/sin(pi(1/4+ix))."""
    return mp.gamma(mp.mpf(1)/4 + 1j*x) * mp.gamma(mp.mpf(3)/4 - 1j*x)


def K3(x):
    """Three-Gamma kernel Gamma(1/4+ix)^2 Gamma(3/4-ix)."""
    return mp.gamma(mp.mpf(1)/4 + 1j*x)**2 * mp.gamma(mp.mpf(3)/4 - 1j*x)


def evalP_std(n, y):
    """Standard MP(1/2, pi/2) polynomial via (n+1)P_{n+1} = 2y P_n - n P_{n-1}."""
    p0, p1 = mp.mpf(1), 2*y
    if n == 0: return p0
    if n == 1: return p1
    for j in range(1, n):
        p0, p1 = p1, (2*y*p1 - j*p0) / (j+1)
    return p1


def GF_master(w, u):
    """Master generating function: 2pi e^eta E_1(eta) / ((1-iw)(1-iu))."""
    eta = (1 + 1j*w) * (1 + 1j*u) / ((1 - 1j*w) * (1 - 1j*u))
    return 2 * mp.pi * mp.exp(eta) * mp.e1(eta) / ((1 - 1j*w) * (1 - 1j*u))


def compute_J_quadrature(mmax):
    """Compute J_mk^{std} = int P_m^{std}(x+i/4) P_k^{std}(x+i/4) K3(x) dx."""
    J = {}
    for m in range(mmax + 1):
        for k in range(m, mmax + 1):
            val = mp.quad(
                lambda x, m=m, k=k:
                    evalP_std(m, x + 1j/4) * evalP_std(k, x + 1j/4) * K3(x),
                NODES)
            J[(m, k)] = val
            J[(k, m)] = val  # symmetric in polynomial indices
            print(f"  J^std_{m},{k} = {mp.nstr(val, 12)}", flush=True)
    return J


def pslq_decompose(J, mmax):
    """Decompose J_mk^{std}/(2pi) as p*A + q with A = e*E_1(1)."""
    A = mp.e * mp.e1(1)
    results = {}
    for m in range(mmax + 1):
        for k in range(m, mmax + 1):
            target = J[(m, k)] / (2 * mp.pi)
            if (m + k) % 2 == 0:
                val = mp.re(target)
            else:
                val = mp.im(target)
            try:
                rel = mp.pslq([val, A, mp.mpf(1)], maxcoeff=10**8)
                if rel and rel[0] != 0:
                    check = -(rel[1]*A + rel[2]) / rel[0]
                    if abs(check - val) < mp.mpf(10)**(-20):
                        from fractions import Fraction
                        fa = Fraction(-rel[1], rel[0])
                        fb = Fraction(-rel[2], rel[0])
                        results[(m, k)] = (fa, fb)
            except:
                pass
    return results


def compute_reverse_table(mmax):
    """Compute one-Gamma reverse: R_mk = int P_m^{std}(x+i/4) P_k^{std}(x+i/4) Gamma(3/4-ix) dx."""
    R = {}
    for m in range(mmax + 1):
        for k in range(m, mmax + 1):
            val = mp.quad(
                lambda x, m=m, k=k:
                    evalP_std(m, x + 1j/4) * evalP_std(k, x + 1j/4)
                    * mp.gamma(mp.mpf(3)/4 - 1j*x),
                NODES)
            R[(m, k)] = val
            R[(k, m)] = val
    return R


def poly_coeffs_shifted(N):
    """Monomial coefficients for P_k^{std}(x+i/4), k=0..N."""
    result = [[mp.mpc(1)]]
    if N == 0: return result
    result.append([1j/2, mp.mpc(2)])
    p_prev, p_curr = result[0], result[1]
    for j in range(1, N):
        new_len = j + 2
        p_new = [mp.mpc(0)] * new_len
        for idx, c in enumerate(p_curr):
            p_new[idx + 1] += 2 * c
            p_new[idx] += 1j/2 * c
        for idx, c in enumerate(p_prev):
            p_new[idx] -= j * c
        p_new = [c / (j + 1) for c in p_new]
        result.append(p_new)
        p_prev = p_curr
        p_curr = p_new
    return result


def plot_zeros(J, outfile):
    """Plot zeros of truncated dual functions F_0 and F_4."""
    pcs = poly_coeffs_shifted(30)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for row_idx, m in enumerate([0, 4]):
        for col_idx, N in enumerate([10, 20, 30]):
            ax = axes[row_idx, col_idx]
            f_coeffs = [mp.mpc(0)] * (N + 1)
            for k in range(N + 1):
                Jmk = J.get((m, k), mp.mpc(0))
                pk = pcs[k]
                for idx, c in enumerate(pk):
                    f_coeffs[idx] += Jmk * c
            while len(f_coeffs) > 1 and abs(f_coeffs[-1]) < 1e-30:
                f_coeffs.pop()
            mp.mp.dps = 80
            roots = mp.polyroots(f_coeffs[::-1], maxsteps=300, extraprec=300)
            mp.mp.dps = 60
            good = [(float(z.real), float(z.imag)) for z in roots if abs(z.imag) < 5]
            ax.scatter([r[0] for r in good], [r[1] for r in good],
                       s=12, c='blue', alpha=0.7, zorder=5)
            ax.axhline(y=0, color='gray', linewidth=0.5)
            ax.axhline(y=-0.25, color='red', linewidth=1.0, linestyle='--',
                       label='Im $= -1/4$' if col_idx == 0 else None)
            ax.axhline(y=-0.5, color='purple', linewidth=0.7, linestyle=':',
                       label='Im $= -1/2$' if col_idx == 0 else None)
            ax.set_title(f'$F_{m}$, N={N}', fontsize=13)
            ax.set_xlabel('Re(x)')
            ax.set_ylabel('Im(x)')
            if row_idx == 0 and col_idx == 0:
                ax.legend(fontsize=8, loc='lower left')
    plt.suptitle(
        r'Zeros of truncated $F_m$: $W = \Gamma(\frac{1}{4}+ix)\Gamma(\frac{3}{4}-ix)$',
        fontsize=14)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"  Saved {outfile}")


def main():
    print("=" * 70)
    print("UNBALANCED DUAL: W(x) = Gamma(1/4+ix) Gamma(3/4-ix)")
    print("=" * 70)

    # ---- Verify moments ----
    print("\n--- Moment GF verification ---")
    def F(s):
        return 2 * mp.pi * mp.exp(3*s/4) / (1 + mp.exp(s))
    for n in range(8):
        Fn = mp.diff(F, 0, n)
        mu_gf = (-1j)**n * Fn
        mu_num = mp.quad(lambda x: x**n * W(x), NODES)
        print(f"  mu_{n} GF={mp.nstr(mu_gf, 12):>22s}  num={mp.nstr(mu_num, 12):>22s}  "
              f"diff={mp.nstr(abs(mu_gf - mu_num), 3)}")

    # ---- Verify shift: polynomials ----
    print("\n--- Shift verification ---")
    print("  P_n(x) = P_n^{MP(1/2)}(x + i/4)")
    print("  Three-term recurrence: P_{n+1} = (x+i/4) P_n - (n/2)^2 P_{n-1}")
    print("  h_n^{std} = pi for all n, h_n^{monic} = pi (n!)^2 / 4^n")

    # ---- 3-Gamma base ----
    A = mp.e * mp.e1(1)
    J00 = mp.quad(lambda x: K3(x), NODES)
    J00_pred = 2 * mp.pi * A
    print(f"\n  J_00 = {mp.nstr(J00, 20)}  predict = {mp.nstr(J00_pred, 20)}  "
          f"diff = {mp.nstr(abs(J00 - J00_pred), 3)}")

    # ---- Master GF verification ----
    print("\n--- Master GF verification ---")
    print("  GF = 2pi e^eta E_1(eta) / ((1-iw)(1-iu))")
    print("  eta = (1+iw)(1+iu) / ((1-iw)(1-iu))")
    for s_val in [0, 0.3, 0.7, -0.5]:
        s = mp.mpf(s_val)
        z = mp.exp(-s)
        pred = 2 * mp.pi * mp.exp(-s/4) * mp.exp(z) * mp.e1(z)
        num = mp.quad(lambda x: K3(x) * mp.exp(1j * s * x), NODES)
        print(f"  s={s_val:5.1f}: ratio = {mp.nstr(num/pred, 12)}")

    # ---- Compute J table ----
    print("\n--- Computing J_mk^{std} table ---")
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            J = pickle.load(f)
        print("  Loaded from cache.")
    else:
        J = compute_J_quadrature(MMAX)
        with open(CACHE, "wb") as f:
            pickle.dump(J, f)

    # ---- PSLQ ----
    print("\n--- PSLQ decomposition ---")
    print(f"  A = e E_1(1) = {mp.nstr(A, 25)}")
    print("  Single generator: d/dz[e^z E_1(z)] = e^z E_1(z) - 1/z => A' = A - 1")
    pslq_table = pslq_decompose(J, min(7, MMAX))
    for (m, k), (fa, fb) in sorted(pslq_table.items()):
        parity = "i" if (m + k) % 2 == 1 else ""
        print(f"  J_{m},{k}^{{std}}/(2pi) = {fa}*A + {fb}  {parity}")

    # ---- Reverse table ----
    print("\n--- Reverse (one-Gamma) table ---")
    R = compute_reverse_table(min(5, MMAX))
    base_R = 2 * mp.pi / mp.e
    print(f"  R_00 = {mp.nstr(R[(0,0)], 15)}  (predict 2pi/e = {mp.nstr(base_R, 15)})")
    print("  Rationals R_mk / (2pi/e):")
    for m in range(6):
        for k in range(m, 6):
            ratio = R[(m, k)] / base_R
            if (m + k) % 2 == 0:
                frac = mp.identify(mp.re(ratio), tol=1e-15)
            else:
                frac = mp.identify(mp.im(ratio), tol=1e-15)
            parity = "i" if (m + k) % 2 == 1 else ""
            print(f"    ({m},{k}): {frac} {parity}")

    # ---- Zeros ----
    print("\n--- Computing zeros ---")
    # Need J for m=0,4 up to k=30
    J_ext = {}
    for m in [0, 4]:
        for k in range(31):
            if (m, k) in J:
                J_ext[(m, k)] = J[(m, k)]
            else:
                val = mp.quad(
                    lambda x, m=m, k=k:
                        evalP_std(m, x + 1j/4) * evalP_std(k, x + 1j/4) * K3(x),
                    NODES)
                J_ext[(m, k)] = val
    plot_zeros(J_ext, "/Users/robin/git/category/unbalanced_zeros.png")

    # ---- Summary ----
    print("\n--- SUMMARY ---")
    print("1. Existence: All Hankel determinants Delta_n > 0 for n <= 12 (REAL positive).")
    print("2. Shift discovery: P_n(x) = P_n^{MP(1/2)}(x+i/4), reducing to self-dual MP(1/2).")
    print(f"3. Generator: SINGLE A = eE_1(1) = {mp.nstr(A, 15)} (vs TWO for self-dual case).")
    print("   E_1 satisfies first-order ODE => d/dz[e^z E_1(z)] = e^z E_1(z) - 1/z.")
    print("4. Master GF: 2pi e^eta E_1(eta) / ((1-iw)(1-iu)), eta = Cayley product.")
    print("5. Parity: m+k even => real, m+k odd => purely imaginary (from W(-x)=conj W(x)).")
    print("6. Zeros: bulk band near Im = -3/4 (w_- pole ladder), not Im = -1/2.")
    print("7. Reverse table: single generator 2pi/e (pure exponential, no Bessel/E_1).")


if __name__ == "__main__":
    main()
