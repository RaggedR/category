"""Zeros of the dual MP approximants  F_m^(N) = sum_{k<=N} c_mk Phat_k.

Coefficients from mp_dual_row_coeffs.pkl (rows m=0,4 to k=40, generated
from the master generating function by mp_dual_row_coeffs.py).  Unlike the
sqrt(W) case (mp_approximant_zeros.py: all zeros real, equioscillation),
here the coefficients are COMPLEX (phase law c_mk = (-i)^{|m-k|} |c_mk|),
so the approximant has no conjugate symmetry and its zeros live in the
complex plane.  The target F_0 = Gamma(1/4+ix) has NO zeros at all --
Gamma never vanishes -- but poles at x = i(1/4+n) in the upper half plane.
Observed at N=11: the zeros accumulate on the line Im x = -1/2, where
Gamma(1/4+ix) = Gamma(3/4+it) -- the lambda=3/4 dual Bump-Ng line.
"""
import pickle
import mpmath as mp
from mp_matching import POLYS

mp.mp.dps = 90

with open("mp_dual_row_coeffs.pkl", "rb") as f:
    C = pickle.load(f)


def h(n):
    return mp.sqrt(2) * mp.pi**mp.mpf(1.5) * mp.binomial(2*n, n) / 4**n


def approximant_coeffs(m, N):
    """Coefficients (in x) of F_m^(N) = sum_{k<=N} c_mk P_k / sqrt(h_k)."""
    coeffs = [mp.mpc(0)] * (N + 1)
    for k in range(N + 1):
        wk = C[(m, k)] / mp.sqrt(h(k))
        for j, co in enumerate(POLYS[k]):
            coeffs[j] += wk * mp.mpf(co.numerator) / mp.mpf(co.denominator)
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def zeros(m, N):
    cs = approximant_coeffs(m, N)
    return mp.polyroots(cs[::-1], maxsteps=500, extraprec=600)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    poles = [0.25 + n for n in range(3)]

    # Left: m=0 -- zeros in the complex plane for increasing N.
    colors = plt.cm.viridis([0.15, 0.45, 0.7, 0.95])
    for col, N in zip(colors, [16, 24, 32, 40]):
        zs = zeros(0, N)
        axes[0].scatter([float(z.real) for z in zs],
                        [float(z.imag) for z in zs],
                        s=18, color=col, label=f"N={N}")
    axes[0].scatter([0] * len(poles), poles, marker="x", s=70, color="crimson",
                    label=r"poles of $\Gamma(\frac{1}{4}+ix)$")
    axes[0].axhline(0, color="gray", lw=0.5)
    axes[0].axhline(-0.5, color="crimson", lw=0.7, ls="--",
                    label=r"Im$\,x=-\frac{1}{2}$  ($\lambda=\frac{3}{4}$ line)")
    axes[0].set_title(r"zeros of $F_0^{(N)}$  ($F_0$ itself has no zeros)")
    axes[0].set_xlabel("Re x")
    axes[0].set_ylabel("Im x")
    axes[0].legend(fontsize=8)

    # Right: m=4, N=40 -- spurious zeros vs the genuine zeros of P_4.
    z4 = zeros(4, 40)
    axes[1].scatter([float(z.real) for z in z4],
                    [float(z.imag) for z in z4],
                    s=18, color="steelblue", label=r"zeros of $F_4^{(40)}$")
    r4 = mp.polyroots([mp.mpf(co.numerator) / mp.mpf(co.denominator)
                       for co in POLYS[4][::-1]])
    axes[1].scatter([float(z.real) for z in r4], [0] * len(r4),
                    s=90, facecolors="none", edgecolors="k",
                    label=r"zeros of $P_4$ (genuine)")
    axes[1].scatter([0] * len(poles), poles, marker="x", s=70, color="crimson")
    axes[1].axhline(0, color="gray", lw=0.5)
    axes[1].axhline(-0.5, color="crimson", lw=0.7, ls="--")
    axes[1].set_title(r"$m=4$, $N=40$")
    axes[1].set_xlabel("Re x")
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("mp_dual_approximant_zeros.png", dpi=150)
    print("saved mp_dual_approximant_zeros.png")

    # how close to Im = -1/2 does the bulk get?
    for N in [16, 24, 32, 40]:
        zs = sorted(zeros(0, N), key=lambda z: float(z.real))
        mid = [z for z in zs if 1 < abs(float(z.real)) < N / 4]
        if mid:
            devs = [abs(float(z.imag) + 0.5) for z in mid]
            print(f"N={N:2d}: bulk zeros {len(mid)}, "
                  f"max |Im+1/2| = {max(devs):.4f}, "
                  f"median = {sorted(devs)[len(devs)//2]:.4f}")

    print("\nzeros of F_0^(40):")
    for z in sorted(zeros(0, 40), key=lambda z: float(z.real)):
        print(f"  {float(z.real): .6f} {float(z.imag):+.6f}i")

    print("\nzeros of F_4^(40) nearest the genuine P_4 zeros "
          "(genuine: +-0.297, +-1.632):")
    for z in sorted(z4, key=lambda z: abs(z)):
        if abs(z) < 2.5:
            print(f"  {float(z.real): .6f} {float(z.imag):+.6f}i")
