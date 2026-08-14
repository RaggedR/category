"""Zeros of the polynomial approximants psi_m^(N) = sum_{k<=N} a_{mk} H_k.

a_{mk} propto D_{mk} / (3^{(m+k)/2} 2^k k!)   (m-fixed overall constants dropped;
they do not affect zeros).  D_{mk} from the matching formula (6) of the paper:

D_{mk} = sum_{c: m-c, k-c even >= 0} C(m,c) C(k,c) c! (m-c-1)!! (k-c-1)!! (-2)^{(m+k-2c)/2}

Coefficients are exact rationals -> build monomial coefficients exactly,
then root-find with mpmath at high precision.
"""
from fractions import Fraction
from math import comb, factorial
import mpmath as mp

mp.mp.dps = 120


def dfact(n):  # double factorial with (-1)!! = 1
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
    """Monomial coefficients of physicist's H_n, as list c[j]*x^j."""
    # H_{n+1} = 2x H_n - 2n H_{n-1}
    a, b = [1], [0, 2]  # H_0, H_1
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
    """Exact monomial coefficients (Fractions) of p_N = sum_k a_{mk} H_k."""
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


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Left: m=0 (Gaussian), several N, zeros scaled by 1/sqrt(2N)
    colors = plt.cm.viridis([0.15, 0.45, 0.7, 0.95])
    for col, N in zip(colors, [20, 30, 40, 56]):
        zs = zeros(0, N)
        sc = 1 / mp.sqrt(2 * N)
        xs = [float(z.real * sc) for z in zs]
        ys = [float(z.imag * sc) for z in zs]
        axes[0].scatter(xs, ys, s=18, color=col, label=f"N={N}")
    axes[0].set_title(r"zeros of $\psi_0^{(N)}$, scaled by $1/\sqrt{2N}$")
    axes[0].axhline(0, color="gray", lw=0.5)
    axes[0].axvline(0, color="gray", lw=0.5)
    axes[0].legend()
    axes[0].set_aspect("equal")

    # Right: m=4, N=40 — genuine zeros of H_4 vs spurious zeros
    N = 40
    zs = zeros(4, N)
    xs = [float(z.real) for z in zs]
    ys = [float(z.imag) for z in zs]
    axes[1].scatter(xs, ys, s=18, color="crimson", label=f"zeros of $\\psi_4^{{(40)}}$")
    import numpy as np
    h4 = np.roots([16, 0, -48, 0, 12])  # H_4 = 16x^4 - 48x^2 + 12
    axes[1].scatter(h4.real, h4.imag, s=90, facecolors="none",
                    edgecolors="k", label="zeros of $H_4$")
    axes[1].set_title(r"$m=4$, $N=40$ (unscaled)")
    axes[1].axhline(0, color="gray", lw=0.5)
    axes[1].legend()
    axes[1].set_aspect("equal")

    plt.tight_layout()
    plt.savefig("approximant_zeros.png", dpi=150)
    print("saved approximant_zeros.png")

    # quantitative check: how close are the 4 genuine zeros?
    genuine = sorted([z for z in zeros(4, N) if abs(z.imag) < 0.1], key=lambda z: z.real)
    print("near-real zeros of psi_4^(40):", [mp.nstr(z, 8) for z in genuine])
    print("zeros of H_4: +-0.52464866, +-1.65068012")
