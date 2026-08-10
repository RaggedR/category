"""Single-panel: zeros of psi_4^(40) with H_4 genuine zeros marked."""
from fractions import Fraction
from math import comb, factorial
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

mp.mp.dps = 120


def dfact(n):
    if n <= 0:
        return 1
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def D(m, k):
    if (m + k) % 2:
        return 0
    tot = 0
    for c in range(0, min(m, k) + 1):
        if (m - c) % 2 or (k - c) % 2:
            continue
        tot += (comb(m, c) * comb(k, c) * factorial(c)
                * dfact(m - c - 1) * dfact(k - c - 1)
                * 4 ** c * (-2) ** ((m + k - 2 * c) // 2))
    return tot


def hermite_coeffs(n):
    a, b = [1], [0, 2]
    if n == 0:
        return a
    for i in range(1, n):
        nxt = [0] * (i + 2)
        for j, cj in enumerate(b):
            nxt[j + 1] += 2 * cj
        for j, cj in enumerate(a):
            nxt[j] -= 2 * i * cj
        a, b = b, nxt
    return b


def approximant_coeffs(m, N):
    poly = [Fraction(0)] * (N + 1)
    for k in range(N + 1):
        d = D(m, k)
        if d == 0:
            continue
        amk = Fraction(d, 3 ** ((m + k) // 2) * 2 ** k * factorial(k))
        for j, cj in enumerate(hermite_coeffs(k)):
            poly[j] += amk * cj
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def zeros(m, N):
    coeffs = approximant_coeffs(m, N)
    cs = [mp.mpf(c.numerator) / mp.mpf(c.denominator) for c in coeffs[::-1]]
    return mp.polyroots(cs, maxsteps=200, extraprec=300)


m, N = 4, 40
print(f"Computing zeros of psi_{m}^({N})...")
zs = zeros(m, N)
xs = [float(z.real) for z in zs]
ys = [float(z.imag) for z in zs]

# Genuine zeros of H_4 = 16x^4 - 48x^2 + 12
h4 = np.roots([16, 0, -48, 0, 12])

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(xs, ys, s=18, color="crimson", label=r"zeros of $\psi_4^{(40)}$", zorder=3)
ax.scatter(h4.real, h4.imag, s=90, facecolors="none",
           edgecolors="k", linewidths=1.5, label=r"zeros of $H_4$ (genuine)", zorder=4)
ax.axhline(0, color="gray", lw=0.5)
ax.axvline(0, color="gray", lw=0.5)
ax.set_title(r"Zeros of $\psi_4^{(40)}$ (Hermite, $m=4$, $N=40$)", fontsize=13)
ax.set_xlabel("Re(z)")
ax.set_ylabel("Im(z)")
ax.legend(fontsize=10)
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("hermite_zeros_single.png", dpi=150)
print("Saved hermite_zeros_single.png")
