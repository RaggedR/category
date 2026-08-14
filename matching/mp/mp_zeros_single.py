"""Single-panel: zeros of Phi_4^(40) (MP lambda=1/4) in the complex plane,
with genuine zeros of P_4 marked as open circles (rug on real axis)."""
import pickle
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import mpmath as mp
mp.mp.dps = 90

# Load matching coefficients and polynomials
from mp_matching import POLYS

with open("mp_matching_coeffs.pkl", "rb") as f:
    NRM, B = pickle.load(f)


def approximant_coeffs(m, N):
    coeffs = [mp.mpf(0)] * (N + 1)
    for k in range(m % 2, N + 1, 2):
        wk = B[(m, k)] / mp.sqrt(NRM[k])
        for j, c in enumerate(POLYS[k]):
            coeffs[j] += wk * mp.mpf(c.numerator) / mp.mpf(c.denominator)
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def zeros(m, N):
    cs = approximant_coeffs(m, N)
    return mp.polyroots(cs[::-1], maxsteps=300, extraprec=400)


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

m, N = 4, 40
print(f"Computing zeros of Phi_{m}^({N})...")
z4 = zeros(m, N)

# All zeros are real for MP lambda=1/4
real_zs = sorted(float(z.real) for z in z4 if abs(z.imag) < 1e-6)
cmplx_zs = [z for z in z4 if abs(z.imag) >= 1e-6]
print(f"  {len(real_zs)} real zeros, {len(cmplx_zs)} complex zeros")

# Genuine zeros of P_4 (lambda=1/4)
p4_coeffs = [mp.mpf(c.numerator) / mp.mpf(c.denominator) for c in POLYS[4][::-1]]
r4 = sorted(float(z.real) for z in mp.polyroots(p4_coeffs) if abs(z.imag) < 1e-6)
print(f"  P_4 genuine zeros: {r4}")

fig, ax = plt.subplots(figsize=(7, 5))

# Scatter all zeros in complex plane
ax.scatter([float(z.real) for z in z4],
           [float(z.imag) for z in z4],
           s=18, color="crimson", label=r"zeros of $\Phi_4^{(40)}$", zorder=3)

# Mark genuine P_4 zeros at Im=0 with open circles
ax.scatter(r4, [0] * len(r4),
           s=90, facecolors="none", edgecolors="k", linewidths=1.5,
           label=r"zeros of $P_4$ (genuine)", zorder=4)

ax.axhline(0, color="gray", lw=0.5)
ax.axvline(0, color="gray", lw=0.5)
ax.set_title(r"Zeros of $\Phi_4^{(40)}$ (MP $\lambda=\frac{1}{4}$, $m=4$, $N=40$)", fontsize=12)
ax.set_xlabel("Re(x)")
ax.set_ylabel("Im(x)")
ax.legend(fontsize=10)
ax.set_xlim(-13, 13)
plt.tight_layout()
plt.savefig("mp_zeros_single.png", dpi=150)
print("Saved mp_zeros_single.png")
