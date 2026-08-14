"""Single-panel: zeros of Laguerre f_{4,40} (alpha=0, row=4, N=40)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import comb, factorial
from fractions import Fraction
import mpmath


mpmath.mp.dps = 40


def gf_coefficients_alpha0_row(row, N):
    """b_{row,k} for alpha=0 via GF 2/(3-w-u-wu), k=0..N."""
    coeffs = []
    m = row
    for k in range(N + 1):
        val = Fraction(0)
        for c in range(min(m, k) + 1):
            n = m + k - c
            coeff = Fraction(factorial(n), factorial(m-c) * factorial(k-c) * factorial(c))
            val += coeff * Fraction(3**c, 3**(m+k))
        val *= Fraction(2, 3)
        coeffs.append(val)
    return coeffs


def compute_zeros(row, N_trunc=40, alpha=0):
    bmk_coeffs = gf_coefficients_alpha0_row(row, N_trunc)
    poly = np.zeros(N_trunc + 1)
    for k in range(N_trunc + 1):
        bk = float(bmk_coeffs[k])
        norm_sq = float(mpmath.gamma(k + alpha + 1) / mpmath.fac(k))
        norm_k = np.sqrt(norm_sq)
        for j in range(k + 1):
            coeff = comb(k, j) * ((-1)**j) / float(factorial(j))
            poly[j] += bk * coeff / norm_k
    poly_coeffs = poly[::-1]
    roots = np.roots(poly_coeffs)
    return roots


row, N = 4, 40
print(f"Computing zeros for row={row}, N={N}...")
roots = compute_zeros(row, N_trunc=N)

real_mask = np.abs(roots.imag) < 0.1
colors = np.where(real_mask, "blue", "red")

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(roots[real_mask].real, roots[real_mask].imag,
           c="blue", s=18, label="real zeros", zorder=3)
ax.scatter(roots[~real_mask].real, roots[~real_mask].imag,
           c="red", s=18, label="complex zeros", zorder=3)
ax.axhline(0, color="gray", lw=0.5)
ax.axvline(0, color="gray", lw=0.5)
ax.set_title(r"Zeros of Laguerre $f_{4,N}$, $N=40$ ($\alpha=0$)", fontsize=13)
ax.set_xlabel("Re(z)")
ax.set_ylabel("Im(z)")
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig("laguerre_zeros_single.png", dpi=150)
print("Saved laguerre_zeros_single.png")
