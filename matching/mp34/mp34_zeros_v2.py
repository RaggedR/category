"""Zeros of truncated dual MP approximants for lambda=3/4.
CORRECTED GF: sum J_mk w^m u^k = pi^2/2 * U(3/2, 1, zeta) / [(1-iw)(1-iu)]^{3/2}
where zeta = (1+iw)(1+iu) / ((1-iw)(1-iu)) = phi(w) * phi(u).

Note: for lambda=1/4 the GF was:
  sum J_mk w^m u^k = 2pi^{3/2} e^{zeta/2} K_0(zeta/2) / sqrt((1-iw)(1-iu))
Both GFs have the SAME (1-iw)(1-iu) denominator structure.
This is because both F_m functions have poles ABOVE the real axis (not below).
"""
import mpmath as mp
import pickle
import os
from fractions import Fraction

mp.mp.dps = 60

with open("mp34_dual_matching.pkl", "rb") as f:
    nrm, C = pickle.load(f)

h0 = nrm[0]
pref = mp.pi**2 / 2

def h34(n):
    return h0 * mp.gamma(n + mp.mpf(3)/2) / (mp.gamma(mp.mpf(3)/2) * mp.fac(n))

def c(m, k):
    return C[(min(m,k), max(m,k))]

def GF_34_correct(w, u):
    """CORRECT GF: pi^2/2 * U(3/2,1,zeta) / [(1-iw)(1-iu)]^{3/2}
    where zeta = (1+iw)(1+iu)/((1-iw)(1-iu)).
    """
    zeta = (1 + 1j*w) * (1 + 1j*u) / ((1 - 1j*w) * (1 - 1j*u))
    exp = mp.mpf(3) / 2
    return pref * mp.hyperu(mp.mpf(3)/2, 1, zeta) / ((1 - 1j*w) * (1 - 1j*u))**exp

# Verify
print(f"GF(0,0) = {mp.nstr(GF_34_correct(0,0), 12)}")
print(f"J_00 = {mp.nstr(c(0,0)*h34(0), 12)}")

# Check J_01:
J01_gf_deriv = mp.diff(lambda u: GF_34_correct(0, u), 0)
J01_quad = c(0,1) * mp.sqrt(h34(0)*h34(1))
print(f"J_01 (GF derivative): {mp.nstr(J01_gf_deriv, 12)}")
print(f"J_01 (quadrature): {mp.nstr(J01_quad, 12)}")

# Extract row coefficients
KMAX = 40
CACHE_ROW = "mp34_row_coeffs_v3.pkl"

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
    vals = [[GF_34_correct(ws[a], us[b]) for b in range(Mu)] for a in range(Mw)]

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

# Verify
max_err = max(abs(C_rows.get((m,k),0) - c(m,k))
              for m in (0,) for k in range(12))
print(f"\nMax error row 0 (k<=11): {mp.nstr(max_err, 5)}")

print("\nDecay of |c_0k|:")
for k in range(0, KMAX+1, 4):
    print(f"  k={k:2d}  |c_0k| = {mp.nstr(abs(C_rows.get((0,k),0)), 8)}")

# =====================================================================
# Polynomial coefficients
# =====================================================================
def mp34_poly_coeffs(nmax):
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

def approximant_coeffs(m, N):
    coeffs = [mp.mpc(0)] * (N + 1)
    for k in range(N + 1):
        wk = C_rows.get((m, k), 0) / mp.sqrt(h34(k))
        for j, co in enumerate(POLYS34[k]):
            coeffs[j] += wk * mp.mpf(co.numerator) / mp.mpf(co.denominator)
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def zeros_of_approx(m, N):
    cs = approximant_coeffs(m, N)
    return mp.polyroots(cs[::-1], maxsteps=500, extraprec=600)

# =====================================================================
# Compute zeros and plot
# =====================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

poles_34 = [0.75 + n for n in range(3)]  # Im part of poles i*(3/4+n)

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

# Analyze accumulation
print("\n=== Zeros analysis ===")
for N, zs in sorted(all_zeros_0.items()):
    neg_z = sorted([z for z in zs if z.imag < -0.05], key=lambda z: z.imag)
    pos_z = sorted([z for z in zs if z.imag > 0.05], key=lambda z: -z.imag)
    print(f"N={N:2d}:")
    if neg_z:
        neg_imag = [float(z.imag) for z in neg_z]
        print(f"  Neg Im: {len(neg_imag)} zeros, range [{min(neg_imag):.4f}, {max(neg_imag):.4f}]")
    if pos_z:
        pos_imag = [float(z.imag) for z in pos_z]
        print(f"  Pos Im: {len(pos_imag)} zeros, range [{min(pos_imag):.4f}, {max(pos_imag):.4f}]")

# Determine the accumulation line
# F_0 = Gamma(3/4+ix) has poles at x = i*(3/4+n) in UPPER half plane.
# The approximant zeros might approach a line below the real axis.
# For lambda=1/4: poles at x = i*(1/4+n), accumulation at Im = -1/2.
# The -1/2 was related to -lambda (or -(lambda+1/2))?
# For lambda=3/4: poles at x = i*(3/4+n), accumulation might be at Im = -3/4?
# Or Im = -1/4 (the other branch)?

# Check: for N=40, where do bulk zeros cluster?
if 40 in all_zeros_0:
    zs40 = all_zeros_0[40]
    # Get all imaginary parts
    all_im = sorted([float(z.imag) for z in zs40])
    print(f"\nN=40 Im parts: {[f'{v:.3f}' for v in sorted(set([round(v,2) for v in all_im]))]}")

# The poles of Gamma(3/4+ix) are at x = i(3/4+n), not marked on Im axis:
# Im x = 3/4, 7/4, 11/4, ...

# Plot reference lines
for y_line, label, style in [
    (-0.25, r"Im$\,x=-\frac{1}{4}$", "crimson"),
    (-0.75, r"Im$\,x=-\frac{3}{4}$", "blue"),
    (-0.5, r"Im$\,x=-\frac{1}{2}$", "green"),
    (0.5, r"Im$\,x=+\frac{1}{2}$", "orange"),
]:
    axes[0].axhline(y_line, color=style, lw=0.8, ls="--", alpha=0.6, label=label)

axes[0].scatter([0]*len(poles_34), poles_34, marker="x", s=70, color="crimson",
                label=r"poles of $\Gamma(\frac{3}{4}+ix)$")
axes[0].axhline(0, color="gray", lw=0.5)
axes[0].set_title(r"zeros of $F_0^{(N)}$ ($\lambda=\frac{3}{4}$, $F_0$ has no zeros)")
axes[0].set_xlabel("Re x")
axes[0].set_ylabel("Im x")
axes[0].set_ylim(-2, 2)
axes[0].legend(fontsize=7, ncol=2)

# Right panel: m=4
print("\nComputing zeros for m=4, N=40...")
try:
    z4 = zeros_of_approx(4, 40)
    axes[1].scatter([float(z.real) for z in z4],
                    [float(z.imag) for z in z4],
                    s=18, color="steelblue", label=r"zeros of $F_4^{(40)}$")
    # Genuine zeros of P_4 (lambda=3/4): real zeros
    p4_coeffs = [mp.mpf(co.numerator)/mp.mpf(co.denominator) for co in POLYS34[4][::-1]]
    r4 = mp.polyroots(p4_coeffs)
    axes[1].scatter([float(z.real) for z in r4], [0]*len(r4),
                    s=90, facecolors="none", edgecolors="k",
                    label=r"zeros of $P_4$ (genuine)")
    print(f"  {len(z4)} zeros found")
    print(f"  Genuine P_4 zeros: {[mp.nstr(z.real, 5) for z in r4 if z.real > 0]}")
except Exception as e:
    print(f"  Error: {e}")

axes[1].scatter([0]*len(poles_34), poles_34, marker="x", s=70, color="crimson")
for y_line, style in [(-0.25, "crimson"), (-0.75, "blue"), (0.5, "orange")]:
    axes[1].axhline(y_line, color=style, lw=0.8, ls="--", alpha=0.6)
axes[1].axhline(0, color="gray", lw=0.5)
axes[1].set_title(r"$m=4$, $N=40$ ($\lambda=\frac{3}{4}$)")
axes[1].set_xlabel("Re x")
axes[1].set_ylim(-2, 2)
axes[1].legend(fontsize=8)

plt.tight_layout()
plt.savefig("mp34_approximant_zeros.png", dpi=150)
print("\nSaved mp34_approximant_zeros.png")

# Precise numerical check
print("\n=== Precise zeros analysis (N=40) ===")
if 40 in all_zeros_0:
    zs40 = all_zeros_0[40]
    bulk_neg = sorted([z for z in zs40 if z.imag < -0.05 and abs(z.real) < 20],
                      key=lambda z: z.imag, reverse=True)
    bulk_pos = sorted([z for z in zs40 if z.imag > 0.05 and abs(z.real) < 20],
                      key=lambda z: -z.imag)
    print(f"Zeros with Im < -0.05 (bulk):")
    for z in bulk_neg[:10]:
        print(f"  {float(z.real):+.5f} + {float(z.imag):.5f}i")
    print(f"\nZeros with Im > 0.05 (bulk):")
    for z in bulk_pos[:10]:
        print(f"  {float(z.real):+.5f} + {float(z.imag):.5f}i")
