"""Zeros of the MP polynomial approximants  Phi_m^(N) = sum_{k<=N} b_{mk} Phat_k.

Uses the cached matching coefficients from mp_matching.py (rows m=0,4 to k=40).
MP analogue of approximant_zeros.py.  With the exponentially decaying weight
|Gamma(1/4+ix)|^2 ~ e^{-pi|x|}, the balance law predicts a spurious-zero onset
LINEAR in N (Hermite case: sqrt(N)), so the left panel scales zeros by 1/N.
"""
import pickle
import mpmath as mp
from mp_matching import POLYS

mp.mp.dps = 90

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


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: m=0 — zero positions vs N.  Nearly all zeros are REAL
    # (contrast Hermite: complex spurious arcs).  Complex ones: open markers.
    colors = plt.cm.viridis([0.15, 0.45, 0.7, 0.95])
    for col, N in zip(colors, [16, 24, 32, 40]):
        zs = zeros(0, N)
        rr = [float(z.real) for z in zs if abs(z.imag) < 1e-6]
        cc = [z for z in zs if z.imag > 1e-6]
        axes[0].scatter(rr, [N] * len(rr), s=14, color=col)
        if cc:
            axes[0].scatter([float(z.real) for z in cc], [N] * len(cc),
                            s=50, facecolors="none", edgecolors="crimson",
                            label=f"complex (N={N}), |Im| up to "
                                  f"{max(float(z.imag) for z in cc):.1f}")
    axes[0].set_title(r"zeros of $\Phi_0^{(N)}$ (nearly all real; linear spread)")
    axes[0].set_xlabel("x (Re for complex)")
    axes[0].set_ylabel("N")
    axes[0].legend(loc="upper left", fontsize=8)

    # Right: m=4, N=40 — three rug levels: approximant, P_40, genuine P_4
    N = 40
    z4 = zeros(4, N)
    a40 = sorted(float(z.real) for z in z4 if abs(z.imag) < 1e-6)
    p40 = mp.polyroots([mp.mpf(c.numerator) / mp.mpf(c.denominator)
                        for c in POLYS[40][::-1]], maxsteps=300, extraprec=400)
    p40 = sorted(float(z.real) for z in p40)
    r4 = mp.polyroots([mp.mpf(c.numerator) / mp.mpf(c.denominator)
                       for c in POLYS[4][::-1]])
    r4 = sorted(float(z.real) for z in r4)
    axes[1].scatter(a40, [2] * len(a40), s=18, color="crimson",
                    label=r"zeros of $\Phi_4^{(40)}$")
    axes[1].scatter(p40, [1] * len(p40), s=18, color="steelblue",
                    label=r"zeros of $P_{40}$")
    axes[1].scatter(r4, [0] * len(r4), s=90, facecolors="none",
                    edgecolors="k", label=r"zeros of $P_4$ (genuine)")
    axes[1].set_yticks([0, 1, 2])
    axes[1].set_yticklabels([r"$P_4$", r"$P_{40}$", r"$\Phi_4^{(40)}$"])
    axes[1].set_ylim(-0.8, 2.8)
    axes[1].set_xlim(-12, 12)
    axes[1].set_title(r"$m=4$, $N=40$: all real; genuine zeros NOT captured")
    axes[1].set_xlabel("x")
    axes[1].legend(loc="upper right", fontsize=8)

    plt.tight_layout()
    plt.savefig("mp_approximant_zeros.png", dpi=150)
    print("saved mp_approximant_zeros.png")

    print("zeros of Phi_4^(40):", ["%.6f" % z for z in a40],
          " (+ %d complex)" % (len(z4) - len(a40)))
    print("zeros of P_4:", ["%.6f" % z for z in r4])
