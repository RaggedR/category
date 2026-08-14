"""Single-panel: complex zeros of F_4^(40) for MP lambda=3/4,
with dashed accumulation line Im x = -1/2."""
import mpmath as mp
import pickle
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))

mp.mp.dps = 60

with open("mp34_dual_matching.pkl", "rb") as f:
    nrm, C = pickle.load(f)

h0 = nrm[0]
pref = mp.pi**2 / 2


def h34(n):
    return h0 * mp.gamma(n + mp.mpf(3)/2) / (mp.gamma(mp.mpf(3)/2) * mp.fac(n))


def GF_34_correct(w, u):
    zeta = (1 + 1j*w) * (1 + 1j*u) / ((1 - 1j*w) * (1 - 1j*u))
    return pref * mp.hyperu(mp.mpf(3)/2, 1, zeta) / ((1 - 1j*w) * (1 - 1j*u))**(mp.mpf(3)/2)


KMAX = 40
CACHE_ROW = "mp34_row_coeffs_v3.pkl"

with open(CACHE_ROW, "rb") as f:
    C_rows = pickle.load(f)
print(f"Loaded {CACHE_ROW}")

from fractions import Fraction


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


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m, N = 4, 40
print(f"Computing zeros of F_{m}^({N}) for MP lambda=3/4...")
z4 = zeros_of_approx(m, N)
print(f"  {len(z4)} zeros found")

# Genuine zeros of P_4 (lambda=3/4)
p4_coeffs = [mp.mpf(co.numerator)/mp.mpf(co.denominator) for co in POLYS34[4][::-1]]
r4 = mp.polyroots(p4_coeffs)

fig, ax = plt.subplots(figsize=(7, 5))

ax.scatter([float(z.real) for z in z4],
           [float(z.imag) for z in z4],
           s=18, color="steelblue", label=r"zeros of $F_4^{(40)}$", zorder=3)

ax.scatter([float(z.real) for z in r4],
           [float(z.imag) for z in r4],
           s=90, facecolors="none", edgecolors="k", linewidths=1.5,
           label=r"zeros of $P_4$ (genuine)", zorder=4)

# Pole markers
poles_34 = [0.75, 1.75, 2.75]
ax.scatter([0]*len(poles_34), poles_34, marker="x", s=70, color="crimson",
           label=r"poles of $\Gamma(\frac{3}{4}+ix)$", zorder=4)

# Accumulation line Im x = -1/2
ax.axhline(-0.5, color="navy", lw=1.2, ls="--",
           label=r"Im$\,x = -\frac{1}{2}$", alpha=0.8)
ax.axhline(0, color="gray", lw=0.5)

ax.set_title(r"Zeros of $F_4^{(40)}$ (MP $\lambda=\frac{3}{4}$, $m=4$, $N=40$)", fontsize=12)
ax.set_xlabel("Re(x)")
ax.set_ylabel("Im(x)")
ax.set_ylim(-2, 2)
ax.legend(fontsize=9, loc="upper right")
plt.tight_layout()
plt.savefig("mp34_zeros_single.png", dpi=150)
print("Saved mp34_zeros_single.png")
