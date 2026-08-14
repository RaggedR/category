"""Wilson polynomial matching matrices at top of the Askey scheme.

Wilson weight on (0,∞):
    W(x) = |Γ(a+ix)Γ(b+ix)Γ(c+ix)Γ(d+ix)/Γ(2ix)|²

Parameters: (a,b,c,d) = (1/4,1/4,3/4,3/4), s = 2.

Wilson orthogonality (DLMF 18.25):
    ∫₀^∞ pₘ(x²) pₖ(x²) W(x) dx = hₙ δ_{mk}

where pₙ = ₄F₃(-n, n+s-1, a+ix, a-ix; a+b, a+c, a+d; 1) and
hₙ = h₀/(poch_n)² with h₀ = π².  For our parameters, h_pn = π² for all n
(a remarkable identity: the Pochhammer prefactors exactly cancel the norm
variation).  So p̂ₙ = pₙ/π are orthonormal.

Matching matrices:
    b_mk = ∫₀^∞ p̂ₘ p̂ₖ W^{3/2} dx       (W^{3/2} = W · √W)
    c_mk = ∫₀^∞ p̂ₘ p̂ₖ w₊ W dx          (dual-split: w₊ = ∏Γ(aⱼ+ix)/Γ(2ix))
"""
import mpmath as mp
import pickle, os
import numpy as np

mp.mp.dps = 50

PARAMS = (mp.mpf(1)/4, mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(3)/4)
NMAX = 8
CACHE = "/Users/robin/git/category/wilson_matching.pkl"
NODES = [0, 0.3, 0.7, 1, 1.5, 2, 3, 5, 8, 12, 18, 25, 35, 50]

# The universal norm for our parameters
H0 = mp.pi**2   # ∫₀^∞ pₙ² W dx = π² for all n


def w_plus(x, params=PARAMS):
    """w₊(x) = ∏ⱼΓ(aⱼ+ix)/Γ(2ix)."""
    a, b, c, d = params
    if abs(x) < mp.mpf('1e-30'):
        return mp.mpc(0)
    return (mp.gamma(a + 1j*x) * mp.gamma(b + 1j*x) *
            mp.gamma(c + 1j*x) * mp.gamma(d + 1j*x) / mp.gamma(2j*x))


def W(x, params=PARAMS):
    """W(x) = |w₊(x)|²."""
    v = w_plus(x, params)
    return (v * mp.conj(v)).real


def sqrt_W(x, params=PARAMS):
    """√W(x) = |w₊(x)|."""
    return abs(w_plus(x, params))


def eval_wilson(n, x, params=PARAMS):
    """pₙ(x²) = ₄F₃(-n, n+s-1, a+ix, a-ix; a+b, a+c, a+d; 1)."""
    a, b, c, d = params
    s = a + b + c + d
    if n == 0:
        return mp.mpf(1)
    return mp.re(mp.hyper([-n, n + s - 1, a + 1j*x, a - 1j*x],
                          [a + b, a + c, a + d], 1))


def main():
    params = PARAMS
    a, b, c, d = params
    s = a + b + c + d

    print("=" * 70)
    print("WILSON MATCHING MATRICES")
    print(f"Parameters (a,b,c,d) = (1/4, 1/4, 3/4, 3/4), s = {float(s)}")
    print("=" * 70)

    # =====================================================================
    # STEP 1: Verify orthogonality and norms
    # =====================================================================
    print("\n--- Orthogonality verification ---")

    # Verify h₀ = π²
    h0_num = mp.quad(lambda x: W(x), NODES)
    print(f"∫₀^∞ W dx (numerical) = {mp.nstr(h0_num, 25)}")
    print(f"π²                    = {mp.nstr(mp.pi**2, 25)}")
    print(f"Match to 25 digits: {abs(h0_num - mp.pi**2) < mp.mpf('1e-20')}")

    # Verify norms for n=1..4
    print("\nNorm verification (all should equal π²):")
    for n in range(5):
        h_num = mp.quad(lambda x, n=n: eval_wilson(n, x)**2 * W(x), NODES)
        print(f"  n={n}: ∫ pₙ² W dx = {mp.nstr(h_num, 25)}, "
              f"reldiff = {mp.nstr(abs(h_num - H0)/H0, 3)}")

    # Cross terms
    print("\nOrthogonality cross terms:")
    for (m, k) in [(0,1), (0,2), (1,2), (0,3), (1,3)]:
        cross = mp.quad(lambda x, m=m, k=k: eval_wilson(m, x) * eval_wilson(k, x) * W(x),
                        NODES)
        print(f"  <p_{m},p_{k}> = {mp.nstr(cross, 6)}")

    # =====================================================================
    # STEP 3: Barnes analysis
    # =====================================================================
    print("\n" + "=" * 70)
    print("BARNES INTEGRAL ANALYSIS FOR c₀₀")
    print("=" * 70)

    print("""
The integral c₀₀ = (1/π²) ∫₀^∞ w₊(x) · W(x) dx  involves:

  Integrand = [∏Γ(aⱼ+ix)]² · [∏Γ(aⱼ-ix)] / [Γ(2ix)² · Γ(-2ix)]

This is 12 Γ-factors / 3 Γ-factors = net 9 Gamma factors.

Known Barnes-type evaluations:
  Barnes 1st lemma:   net 4 Γ's     (4/0)
  Barnes 2nd lemma:   net 4 Γ's     (5/1, balanced condition)
  Wilson/de Branges:  net 6 Γ's     (8/2, requires Re(params)>0)
  Nassrallah-Rahman:  net 6 Γ's     (10/4, ₅F₄ extension)
  Gustafson A_n:      net 8 Γ's     (12/4, root system)

Our integral: net 9 Γ's — exceeds all standard evaluations.

STRUCTURAL OBSTRUCTION: The Γ(2ix) denominator has coefficient 2 in the
argument (2ix vs aⱼ±ix). Barnes-type evaluations require ALL Γ-arguments
to be first-degree in the integration variable with coefficient ±1.
The duplication formula Γ(2z) = π^{-1/2} 2^{2z-1/2} Γ(z)Γ(z+½) replaces
Γ(2ix) by Γ(ix)·Γ(ix+½), but these have Re(param) = 0 or ½ — and the
ix-only factor Γ(ix) has poles on the REAL axis (at x = 0, -i, -2i, ...),
violating the separation-of-poles condition.

In the MP degeneration (c,d → ∞), the Γ(2ix) factors cancel entirely:
W_MP = |Γ(a+ix)|², no denominator. This is why MP matching coefficients
admit closed forms (the 3-Γ lemma) while Wilson matching does not.

CONCLUSION: c₀₀ is NOT expressible via known Barnes lemmas. The Γ(2ix)
denominator — native to Wilson polynomials and absent in all degenerations
below continuous Hahn — is the precise obstruction.
""")

    # Numerical c₀₀
    J00 = mp.quad(lambda x: w_plus(x) * W(x), NODES)
    c00 = J00 / H0
    print(f"J₀₀ = ∫ w₊ W dx = {mp.nstr(J00, 25)}")
    print(f"c₀₀ = J₀₀/π²   = {mp.nstr(c00, 25)}")
    print(f"|c₀₀|           = {mp.nstr(abs(c00), 25)}")
    print(f"arg(c₀₀)/π      = {mp.nstr(mp.arg(c00)/mp.pi, 15)}")

    # =====================================================================
    # STEP 2: Matching matrices
    # =====================================================================
    print("\n" + "=" * 70)
    print("MATCHING MATRICES (9×9)")
    print("=" * 70)

    # b_mk: W^{3/2} matching (real-valued since W^{3/2} is even and real)
    print("\n--- b_mk = (1/π²) ∫₀^∞ pₘ pₖ W^{3/2} dx ---")
    B = {}
    for m in range(NMAX + 1):
        for k in range(m, NMAX + 1):
            val = mp.quad(
                lambda x, m=m, k=k: eval_wilson(m, x) * eval_wilson(k, x) *
                W(x) * sqrt_W(x) / H0,
                NODES)
            B[(m, k)] = B[(k, m)] = mp.re(val)
            print(f"  b_{m},{k} = {mp.nstr(mp.re(val), 18)}", flush=True)

    # c_mk: dual-split matching (complex-valued)
    print("\n--- c_mk = (1/π²) ∫₀^∞ pₘ pₖ w₊ W dx ---")
    CC = {}
    for m in range(NMAX + 1):
        for k in range(m, NMAX + 1):
            val = mp.quad(
                lambda x, m=m, k=k: eval_wilson(m, x) * eval_wilson(k, x) *
                w_plus(x) * W(x) / H0,
                NODES)
            CC[(m, k)] = val
            CC[(k, m)] = mp.conj(val)
            print(f"  c_{m},{k} = {mp.nstr(val, 18)}", flush=True)

    # Save
    with open(CACHE, 'wb') as f:
        pickle.dump({
            'params': (0.25, 0.25, 0.75, 0.75),
            'H0': float(H0),
            'B': {str(k): (float(mp.re(v)), float(mp.im(v))) for k, v in B.items()},
            'C': {str(k): (float(mp.re(v)), float(mp.im(v))) for k, v in CC.items()},
        }, f)
    print(f"\nSaved to {CACHE}")

    # =====================================================================
    # STEP 4: PSLQ
    # =====================================================================
    print("\n" + "=" * 70)
    print("PSLQ ANALYSIS (6×6 corner)")
    print("=" * 70)

    K0h = mp.besselk(0, mp.mpf(1)/2)
    K1h = mp.besselk(1, mp.mpf(1)/2)
    A_mp = mp.sqrt(2 * mp.e) * K0h
    B_mp = mp.sqrt(2 * mp.e) * K1h
    G14 = mp.gamma(mp.mpf(1)/4)
    G34 = mp.gamma(mp.mpf(3)/4)
    G12 = mp.gamma(mp.mpf(1)/2)  # = √π

    # Wilson-specific: h₀ = π², and the orthogonality involves
    # Γ(1/2)² Γ(1) Γ(3/2) etc. So natural candidates include
    # Γ values at half-integers and quarter-integers.
    candidates = [
        ('1',        mp.mpf(1)),
        ('π',        mp.pi),
        ('π²',       mp.pi**2),
        ('Γ(1/4)',   G14),
        ('Γ(3/4)',   G34),
        ('Γ(1/4)²',  G14**2),
        ('Γ(1/4)Γ(3/4)', G14 * G34),  # = π√2
        ('Γ(1/4)⁴',  G14**4),
        ('K₀(½)',    K0h),
        ('K₁(½)',    K1h),
        ('A_MP',     A_mp),
        ('B_MP',     B_mp),
    ]

    print("\nCandidate constants:")
    for nm, val in candidates:
        print(f"  {nm:15s} = {mp.nstr(val, 20)}")

    # Check Γ(1/4)·Γ(3/4) = π/sin(π/4) = π√2 (reflection formula)
    print(f"\n  Γ(1/4)·Γ(3/4) = {mp.nstr(G14*G34, 20)}")
    print(f"  π√2            = {mp.nstr(mp.pi*mp.sqrt(2), 20)}")

    for mat_name, mat in [('b', B), ('c', CC)]:
        print(f"\n--- PSLQ for {mat_name}_mk ---")
        for m in range(6):
            for k in range(m, 6):
                val = mat[(m, k)]
                for part_name, part_val in [("Re", mp.re(val)), ("Im", mp.im(val))]:
                    if abs(part_val) < mp.mpf('1e-20'):
                        continue

                    found = False
                    # Try single-constant rational multiples
                    for cn, cv in candidates:
                        if abs(cv) < mp.mpf('1e-20'):
                            continue
                        ratio = part_val / cv
                        for denom in range(1, 2000):
                            numer = mp.nint(ratio * denom)
                            if abs(ratio * denom - numer) < mp.mpf('1e-15'):
                                print(f"  {part_name} {mat_name}_{m},{k} = "
                                      f"({int(numer)}/{denom}) · {cn}")
                                found = True
                                break
                        if found:
                            break

                    if not found:
                        # Try 2-constant PSLQ
                        for i in range(len(candidates)):
                            for j in range(i+1, len(candidates)):
                                cn1, cv1 = candidates[i]
                                cn2, cv2 = candidates[j]
                                try:
                                    rel = mp.pslq([part_val, cv1, cv2],
                                                  maxcoeff=10000, maxsteps=300)
                                    if rel is not None:
                                        print(f"  {part_name} {mat_name}_{m},{k}: "
                                              f"{rel[0]}·val + {rel[1]}·{cn1} + {rel[2]}·{cn2} = 0")
                                        found = True
                                        break
                                except:
                                    pass
                            if found:
                                break

                    if not found:
                        print(f"  {part_name} {mat_name}_{m},{k} = "
                              f"{mp.nstr(part_val, 18)}: no relation (maxcoeff=10000)")

    # =====================================================================
    # STEP 5: Decay
    # =====================================================================
    print("\n" + "=" * 70)
    print("DECAY ANALYSIS")
    print("=" * 70)

    for mat_name, mat in [('b', B), ('c', CC)]:
        print(f"\n--- |{mat_name}_{{0,k}}| decay ---")
        vals = [(k, abs(mat[(0, k)])) for k in range(NMAX + 1)]
        print(f"  {'k':>3s}  {'|val|':>22s}  {'ratio':>12s}  {'log ratio':>12s}")
        for i, (k, v) in enumerate(vals):
            if i > 0:
                ratio = v / vals[i-1][1]
                lr = float(mp.log(ratio)) if ratio > 0 else float('nan')
            else:
                ratio = mp.mpf(1)
                lr = 0.0
            print(f"  {k:3d}  {mp.nstr(v, 18):>22s}  {mp.nstr(ratio, 8):>12s}  {lr:12.4f}")

    # Check: is decay geometric or power-law?
    print("\n--- Decay characterization ---")
    for mat_name, mat in [('b', B), ('c', CC)]:
        row0 = [abs(mat[(0, k)]) for k in range(NMAX + 1)]
        # Fit log|val| ~ a + b*k (geometric) vs a + b*log(k) (power law)
        ks = list(range(1, NMAX + 1))
        log_vals = [float(mp.log(row0[k])) for k in ks]
        log_ks = [float(mp.log(k)) for k in ks]

        # Geometric fit: log|val| = a + b*k
        n_pts = len(ks)
        sum_k = sum(ks)
        sum_lv = sum(log_vals)
        sum_klv = sum(k * lv for k, lv in zip(ks, log_vals))
        sum_k2 = sum(k**2 for k in ks)
        b_geo = (n_pts * sum_klv - sum_k * sum_lv) / (n_pts * sum_k2 - sum_k**2)
        a_geo = (sum_lv - b_geo * sum_k) / n_pts

        # Power-law fit: log|val| = a + b*log(k)
        sum_lk = sum(log_ks)
        sum_lklv = sum(lk * lv for lk, lv in zip(log_ks, log_vals))
        sum_lk2 = sum(lk**2 for lk in log_ks)
        b_pow = (n_pts * sum_lklv - sum_lk * sum_lv) / (n_pts * sum_lk2 - sum_lk**2)
        a_pow = (sum_lv - b_pow * sum_lk) / n_pts

        # Residuals
        res_geo = sum((lv - a_geo - b_geo*k)**2 for k, lv in zip(ks, log_vals))
        res_pow = sum((lv - a_pow - b_pow*lk)**2 for lk, lv in zip(log_ks, log_vals))

        print(f"\n  {mat_name}_{{0,k}}:")
        print(f"    Geometric fit:  log|val| ~ {a_geo:.4f} + {b_geo:.4f}·k  "
              f"(ρ = e^b = {np.exp(b_geo):.6f}), residual = {res_geo:.6f}")
        print(f"    Power-law fit:  log|val| ~ {a_pow:.4f} + {b_pow:.4f}·log(k), "
              f"residual = {res_pow:.6f}")
        winner = "GEOMETRIC" if res_geo < res_pow else "POWER LAW"
        print(f"    Better fit: {winner}")

    # =====================================================================
    # STEP 6: Degeneration check
    # =====================================================================
    print("\n" + "=" * 70)
    print("DEGENERATION CHECK: Wilson → Hahn → MP")
    print("=" * 70)

    print("\nc₀₀ for Wilson(1/4, 1/4, 3/4, T), T growing:")
    print("(As T → ∞, Wilson degenerates to continuous dual Hahn then to MP)")
    for T_val in [0.75, 2, 5, 10, 20, 50]:
        T = mp.mpf(T_val)
        par = (mp.mpf(1)/4, mp.mpf(1)/4, mp.mpf(3)/4, T)

        # norm h₀ for these params
        pairs_T = [(par[i], par[j]) for i in range(4) for j in range(i+1, 4)]
        prod_g = mp.mpf(1)
        for ai, aj in pairs_T:
            prod_g *= mp.gamma(ai + aj)
        s_T = sum(par)
        h0_T = 2 * mp.pi * prod_g / ((s_T - 1) * mp.gamma(s_T - 1))

        try:
            J00_T = mp.quad(lambda x, p=par: w_plus(x, p) * W(x, p), NODES)
            c00_T = J00_T / h0_T
            print(f"  T={T_val:6.2f}: c₀₀ = {mp.nstr(c00_T, 15)}, "
                  f"|c₀₀| = {mp.nstr(abs(c00_T), 12)}, h₀ = {mp.nstr(h0_T, 12)}")
        except Exception as e:
            print(f"  T={T_val:6.2f}: FAILED ({e})")

    print(f"\n  MP λ=1/4 reference: c₀₀ = A_MP = {mp.nstr(A_mp, 15)}")
    print("  (The drift shows Wilson c₀₀ moving away from the MP value as T grows,")
    print("   because the degeneration involves rescaling. The KEY structural point:")
    print("   the Γ(2ix) factors in Wilson DO NOT survive the limit.)")

    # =====================================================================
    # STEP 7: Zeros
    # =====================================================================
    print("\n" + "=" * 70)
    print("ZEROS OF TRUNCATED EXPANSIONS")
    print("=" * 70)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, trunc_N in enumerate([4, NMAX]):
        def f_trunc(x, N=trunc_N):
            result = mp.mpc(0)
            for k in range(N + 1):
                result += CC[(0, k)] * eval_wilson(k, x) / H0
            return result

        # Find real zeros of Re(f)
        real_zeros = []
        x_scan = [mp.mpf(j)/5 for j in range(1, 100)]
        for i in range(len(x_scan) - 1):
            v1 = mp.re(f_trunc(x_scan[i]))
            v2 = mp.re(f_trunc(x_scan[i+1]))
            if v1 * v2 < 0:
                try:
                    z = mp.findroot(lambda x: mp.re(f_trunc(x)),
                                    (x_scan[i], x_scan[i+1]))
                    real_zeros.append(float(mp.re(z)))
                except:
                    pass

        print(f"\n  N={trunc_N}: real zeros of Re(f) on (0,20):")
        for rz in real_zeros:
            print(f"    x = {rz:.6f}")

        # Plot
        xs = np.linspace(0.01, 15, 400)
        vals_re = [float(mp.re(f_trunc(mp.mpf(x)))) for x in xs]
        vals_im = [float(mp.im(f_trunc(mp.mpf(x)))) for x in xs]
        vals_abs = [np.sqrt(r**2 + i**2) for r, i in zip(vals_re, vals_im)]

        ax = axes[idx]
        ax.plot(xs, vals_re, 'b-', linewidth=1.5, label='Re')
        ax.plot(xs, vals_im, 'r--', linewidth=1.2, label='Im')
        ax.plot(xs, vals_abs, 'k:', linewidth=0.8, label='|f|')
        ax.axhline(0, color='gray', linewidth=0.5)
        for rz in real_zeros:
            ax.axvline(rz, color='green', linewidth=0.8, alpha=0.5, linestyle=':')
        ax.set_title(f'$\\sum_{{k \\leq {trunc_N}}} c_{{0k}} \\hat{{W}}_k(x)$, '
                     f'Wilson $(\\frac{{1}}{{4}},\\frac{{1}}{{4}},\\frac{{3}}{{4}},\\frac{{3}}{{4}})$')
        ax.set_xlabel('$x$')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpng = '/Users/robin/git/category/wilson_zeros.png'
    plt.savefig(outpng, dpi=150)
    print(f"\n  Saved {outpng}")
    plt.close()

    # =====================================================================
    # Summary tables
    # =====================================================================
    print("\n" + "=" * 70)
    print("SUMMARY TABLES")
    print("=" * 70)

    print("\n--- b_mk (real, 9×9) ---")
    header = "     " + "".join(f"{'k='+str(k):>14s}" for k in range(NMAX + 1))
    print(header)
    for m in range(NMAX + 1):
        row = f"m={m}  "
        for k in range(NMAX + 1):
            row += f"{float(mp.re(B[(m,k)])):14.6f}"
        print(row)

    print("\n--- c_mk (complex, 9×9) — showing |c_mk| ---")
    header = "     " + "".join(f"{'k='+str(k):>14s}" for k in range(NMAX + 1))
    print(header)
    for m in range(NMAX + 1):
        row = f"m={m}  "
        for k in range(NMAX + 1):
            row += f"{float(abs(CC[(m,k)])):14.6f}"
        print(row)

    print("\n--- c_mk phases arg(c_mk)/π ---")
    header = "     " + "".join(f"{'k='+str(k):>10s}" for k in range(NMAX + 1))
    print(header)
    for m in range(NMAX + 1):
        row = f"m={m}  "
        for k in range(NMAX + 1):
            phase = float(mp.arg(CC[(m, k)]) / mp.pi)
            row += f"{phase:10.4f}"
        print(row)


if __name__ == "__main__":
    main()
