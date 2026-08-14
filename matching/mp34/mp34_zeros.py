"""Zeros of truncated dual MP approximants for lambda=3/4.

F_m^(N) = sum_{k<=N} c_mk Phat_k(x).

For lambda=1/4: zeros converged to Im x = -1/2 (the lambda=3/4 line).
For lambda=3/4: expect Im x = +/- something related to next lambda, or the poles of Gamma(3/4+ix).
Gamma(3/4+ix) has poles at x = i(3/4+n) for n=0,1,2,...
So F_0 = Gamma(3/4+ix) has poles at x = 3i/4, 7i/4, ...

The zeros of the approximant should converge to some line Im x = const.
For lambda=1/4 the line was Im x = -1/2 = -(3/4 - 1/4).
Conjecture for lambda=3/4: Im x = -something.

Note: Gamma(3/4-ix) has poles at x = -i(3/4+n) (below real axis).
The "dual" structure suggests: zeros of F_m^(N) for lambda=3/4 might go to
Im x = -1/4 (the lambda=1/4 partner) or Im x = -1.
"""
import mpmath as mp
import pickle
import os

mp.mp.dps = 60

# Load the quadrature cache
with open("mp34_dual_matching.pkl", "rb") as f:
    nrm, C = pickle.load(f)

h0 = nrm[0]
pref = mp.pi**2 / 2

def h34(n):
    return h0 * mp.gamma(n + mp.mpf(3)/2) / (mp.gamma(mp.mpf(3)/2) * mp.fac(n))

def c(m, k):
    return C[(min(m,k), max(m,k))]

def GF_34(w, u):
    """CORRECTED GF: pi^2/2 * U(3/2,1,q) / ((1+iw)(1+iu))^{3/2}."""
    q = (1 - 1j*w) * (1 - 1j*u) / ((1 + 1j*w) * (1 + 1j*u))
    exp = mp.mpf(3) / 2
    return pref * mp.hyperu(mp.mpf(3)/2, 1, q) / ((1 + 1j*w) * (1 + 1j*u))**exp

# Verify the fix
print(f"GF(0,0) = {mp.nstr(GF_34(0,0), 12)}")
print(f"J_00 = {mp.nstr(c(0,0)*h34(0), 12)}")

# =====================================================================
# Extract row coefficients via corrected GF
# =====================================================================
KMAX = 40
CACHE_ROW = "mp34_row_coeffs_v2.pkl"

if os.path.exists(CACHE_ROW):
    with open(CACHE_ROW, "rb") as f:
        C_rows = pickle.load(f)
    print(f"\nLoaded {CACHE_ROW}")
else:
    print("\nComputing row coefficients via corrected GF...")
    Mw, Mu = 64, 128
    r = mp.mpf(1) / 4
    ws = [r * mp.expjpi(2 * mp.mpf(a) / Mw) for a in range(Mw)]
    us = [r * mp.expjpi(2 * mp.mpf(b) / Mu) for b in range(Mu)]
    vals = [[GF_34(ws[a], us[b]) for b in range(Mu)] for a in range(Mw)]

    C_rows = {}
    for m in (0, 4):
        inner = [sum(vals[a][b] * mp.expjpi(-2*mp.mpf(m*a)/Mw) for a in range(Mw)) / (Mw * r**m)
                 for b in range(Mu)]
        for k in range(KMAX+1):
            J_gf = sum(inner[b] * mp.expjpi(-2*mp.mpf(k*b)/Mu) for b in range(Mu)) / (Mu * r**k)
            C_rows[(m, k)] = J_gf / mp.sqrt(h34(m) * h34(k))

    with open(CACHE_ROW, "wb") as f:
        pickle.dump(C_rows, f)
    print(f"Saved {CACHE_ROW}")

# Verify against quadrature
print("\nVerification vs quadrature (rows 0,4, k<=11):")
for m in (0, 4):
    for k in range(min(12, m), 12):
        c_gf = C_rows.get((m, k), 0)
        c_quad = c(m, k) if m <= k else c(k, m)
        err = abs(c_gf - c_quad)
        if err > 1e-10:
            print(f"  c_{m},{k}: GF={mp.nstr(abs(c_gf),6)}, quad={mp.nstr(abs(c_quad),6)}, err={mp.nstr(err,3)}")

max_err = max(abs(C_rows.get((m,k),0) - (c(m,k) if m<=k else c(k,m)))
              for m in (0,4) for k in range(12))
print(f"Max error: {mp.nstr(max_err, 5)}")

print("\nDecay of |c_0k|:")
for k in range(0, KMAX+1, 4):
    print(f"  k={k:2d}  |c_0k| = {mp.nstr(abs(C_rows.get((0,k),0)), 8)}")

# =====================================================================
# Polynomial coefficients for lambda=3/4
# =====================================================================
from fractions import Fraction

def mp34_poly_coeffs(nmax):
    """Monomial coefficients of P_n for lambda=3/4."""
    polys = [[Fraction(1)], [Fraction(0), Fraction(2)]]
    for n in range(1, nmax):
        pn, pm = polys[n], polys[n-1]
        nxt = [Fraction(0)] * (n + 2)
        for j, c_coeff in enumerate(pn):
            nxt[j+1] += 2 * c_coeff
        for j, c_coeff in enumerate(pm):
            nxt[j] -= (Fraction(n) + Fraction(1, 2)) * c_coeff
        polys.append([x / (n + 1) for x in nxt])
    return polys

POLYS34 = mp34_poly_coeffs(KMAX + 2)

def approximant_coeffs(m, N, row_cache=C_rows):
    """Coefficients of F_m^(N) = sum_{k<=N} c_mk P_k / sqrt(h_k)."""
    coeffs = [mp.mpc(0)] * (N + 1)
    for k in range(N + 1):
        wk = row_cache.get((m, k), 0) / mp.sqrt(h34(k))
        for j, co in enumerate(POLYS34[k]):
            coeffs[j] += wk * mp.mpf(co.numerator) / mp.mpf(co.denominator)
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def zeros_of_approx(m, N):
    cs = approximant_coeffs(m, N)
    return mp.polyroots(cs[::-1], maxsteps=500, extraprec=600)

# =====================================================================
# Compute zeros
# =====================================================================
print("\n=== Computing zeros of approximants ===")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

poles_34 = [0.75 + n for n in range(3)]  # poles of Gamma(3/4+ix) at x = i*(3/4+n)

# Left: m=0 for N = 10, 20, 40
colors = plt.cm.viridis([0.15, 0.45, 0.7, 0.95])
all_zeros_0 = {}
for col, N in zip(colors, [10, 20, 30, 40]):
    print(f"Computing zeros for m=0, N={N}...")
    try:
        zs = zeros_of_approx(0, N)
        all_zeros_0[N] = zs
        axes[0].scatter([float(z.real) for z in zs],
                        [float(z.imag) for z in zs],
                        s=18, color=col, label=f"N={N}")
        print(f"  {len(zs)} zeros found")
    except Exception as e:
        print(f"  Error: {e}")

axes[0].scatter([0]*len(poles_34), poles_34, marker="x", s=70, color="crimson",
                label=r"poles of $\Gamma(\frac{3}{4}+ix)$")
axes[0].axhline(0, color="gray", lw=0.5)

# Find where zeros accumulate
print("\nAnalyzing Im part of bulk zeros:")
for N, zs in all_zeros_0.items():
    # Filter to bulk (away from real axis and not near poles)
    bulk = [z for z in zs if abs(z.imag) > 0.1 and abs(z.imag) < 5]
    if bulk:
        imag_vals = [float(z.imag) for z in bulk]
        median_im = sorted(imag_vals)[len(imag_vals)//2]
        print(f"  N={N:2d}: {len(bulk)} bulk zeros, median Im = {median_im:.4f}")
        # Check if they cluster below real axis at some Im = -c
        neg_bulk = [float(z.imag) for z in bulk if z.imag < 0]
        if neg_bulk:
            print(f"         negative Im: min={min(neg_bulk):.4f}, median={sorted(neg_bulk)[len(neg_bulk)//2]:.4f}")

# Determine the accumulation line
# For lambda=1/4, accumulation was at Im x = -1/2.
# For lambda=3/4, the dual pole structure suggests Im x = -1/4 (lambda=1/4 line).
# Or: Im x = -(1-3/4) = -1/4? Or Im x = -3/4? Or Im x = -(lambda_next)?
# From the lambda=1/4 result: "zeros accumulate on Im x = -1/2 where Gamma(1/4+ix) = Gamma(3/4+it)
# -- the lambda=3/4 dual Bump-Ng line". So for lambda=3/4, the accumulation should be at
# Im x = -3/4? (the lambda=3/4+1/2 = 5/4 dual?) Or at -1/4 (the other Bump-Ng line)?

# Let's look at where zeros actually are
if all_zeros_0:
    N_max = max(all_zeros_0.keys())
    zs_max = all_zeros_0[N_max]
    bulk_neg = sorted([z for z in zs_max if -5 < z.imag < -0.05],
                      key=lambda z: abs(z.imag))
    if bulk_neg:
        print(f"\nZeros with -5 < Im z < -0.05 (N={N_max}):")
        imag_vals = [float(z.imag) for z in bulk_neg]
        # Find the mode
        from collections import Counter
        rounded = [round(v, 1) for v in imag_vals]
        print(f"  Rounded Im values: {Counter(rounded).most_common(5)}")
        print(f"  First few: {imag_vals[:10]}")

# Try annotating possible lines
for y_line, label in [(-0.25, "Im=-1/4"), (-0.75, "Im=-3/4"), (-1.0, "Im=-1"), (-0.5, "Im=-1/2")]:
    axes[0].axhline(y_line, color="gray", lw=0.4, ls=":", alpha=0.5)
    axes[0].text(8, y_line + 0.05, label, fontsize=6, color="gray")

axes[0].set_title(r"zeros of $F_0^{(N)}$ ($\lambda=\frac{3}{4}$)")
axes[0].set_xlabel("Re x")
axes[0].set_ylabel("Im x")
axes[0].set_ylim(-3, 3)
axes[0].legend(fontsize=8)

# Right: m=4, N=40
print("\nComputing zeros for m=4, N=40...")
try:
    z4 = zeros_of_approx(4, 40)
    axes[1].scatter([float(z.real) for z in z4],
                    [float(z.imag) for z in z4],
                    s=18, color="steelblue", label=r"zeros of $F_4^{(40)}$")

    # Genuine zeros of P_4 (lambda=3/4)
    p4_coeffs = [mp.mpf(co.numerator)/mp.mpf(co.denominator) for co in POLYS34[4][::-1]]
    r4 = mp.polyroots(p4_coeffs)
    axes[1].scatter([float(z.real) for z in r4], [0]*len(r4),
                    s=90, facecolors="none", edgecolors="k",
                    label=r"zeros of $P_4$ (genuine)")
except Exception as e:
    print(f"  Error: {e}")

axes[1].scatter([0]*len(poles_34), poles_34, marker="x", s=70, color="crimson")
axes[1].axhline(0, color="gray", lw=0.5)
for y_line in [-0.25, -0.75, -1.0]:
    axes[1].axhline(y_line, color="gray", lw=0.4, ls=":", alpha=0.5)
axes[1].set_title(r"$m=4$, $N=40$ ($\lambda=\frac{3}{4}$)")
axes[1].set_xlabel("Re x")
axes[1].set_ylim(-3, 3)
axes[1].legend(fontsize=8)

plt.tight_layout()
plt.savefig("mp34_approximant_zeros.png", dpi=150)
print("\nSaved mp34_approximant_zeros.png")

# Precise analysis of accumulation line
print("\n=== Zeros analysis for m=0 ===")
for N, zs in sorted(all_zeros_0.items()):
    # Find the minimum |Im z| for negative imaginary zeros
    neg = [z for z in zs if z.imag < -0.01]
    if neg:
        closest = min(neg, key=lambda z: abs(z.imag))
        # The bulk: zoom into the region where zeros cluster
        bulk = sorted([z for z in neg if abs(z.real) < N/3 and z.imag > -3],
                      key=lambda z: abs(z.imag))
        if bulk:
            imag_bulk = [float(z.imag) for z in bulk]
            print(f"  N={N:2d}: {len(bulk)} bulk zeros in lower half")
            print(f"    Im range: [{min(imag_bulk):.4f}, {max(imag_bulk):.4f}]")
            print(f"    Median Im: {sorted(imag_bulk)[len(imag_bulk)//2]:.4f}")
