"""Single-panel: zeros of Chebyshev-T A_40, with Bernstein ellipses and
genuine/spurious distinction."""
import mpmath as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 50


def I_sin_cos(p, q):
    p = mp.mpf(p)
    q = int(abs(q))
    num = mp.pi * mp.gamma(p + 1)
    den = mp.power(2, p) * mp.gamma((p + q + 2)/2) * mp.gamma((p - q + 2)/2)
    return num / den


def b_T_0k(k):
    p = mp.mpf(-1)/2
    integral = I_sin_cos(p, k)
    nm = 1/mp.sqrt(mp.pi)
    nk = 1/mp.sqrt(mp.pi) if k == 0 else mp.sqrt(2/mp.pi)
    return nm * nk * integral


def cheb_T_zeros(N):
    c = []
    for k in range(N+1):
        b = float(b_T_0k(k))
        if k == 0:
            c.append(b / np.sqrt(np.pi))
        else:
            c.append(b * np.sqrt(2/np.pi))
    poly = np.polynomial.chebyshev.Chebyshev(c)
    return poly.roots()


def bernstein_ellipse(rho, npts=500):
    a = (rho + 1/rho) / 2
    b = (rho - 1/rho) / 2
    theta = np.linspace(0, 2*np.pi, npts)
    return a * np.cos(theta), b * np.sin(theta)


N = 40
print(f"Computing Chebyshev-T zeros, N={N}...")
zs = cheb_T_zeros(N)
print(f"  {len(zs)} zeros computed")

real_mask = (np.abs(zs.imag) < 1e-10) & (np.abs(zs.real) <= 1.01)
spur_mask = ~real_mask

fig, ax = plt.subplots(figsize=(7, 5))

ax.scatter(zs[real_mask].real, zs[real_mask].imag,
           s=18, color="tab:blue", label=f"genuine (in $[-1,1]$)", zorder=4)
ax.scatter(zs[spur_mask].real, zs[spur_mask].imag,
           s=18, color="tab:blue", facecolors="none",
           marker="o", linewidths=0.8,
           label=f"spurious", zorder=4)

for rho, ls, lw in [(1.1, ':', 1.0), (1.5, '--', 1.0), (2.0, '-', 1.0)]:
    xe, ye = bernstein_ellipse(rho)
    ax.plot(xe, ye, ls, color='gray', lw=lw, alpha=0.6,
            label=f'Bernstein $\\rho={rho}$')

ax.axhline(0, color='k', linewidth=0.5)
ax.axvline(-1, color='k', linestyle=':', linewidth=0.8, alpha=0.4)
ax.axvline(1, color='k', linestyle=':', linewidth=0.8, alpha=0.4)
ax.set_xlim(-3.5, 3.5)
ax.set_ylim(-2.5, 2.5)
ax.set_aspect("equal")
ax.set_xlabel("Re(z)", fontsize=11)
ax.set_ylabel("Im(z)", fontsize=11)
ax.set_title(f"Zeros of Chebyshev-T $A_{{40}}$ ($\\alpha=\\beta=-\\frac{{1}}{{2}}$, $N={N}$)", fontsize=12)
ax.legend(fontsize=9, ncol=2)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("cheb_zeros_single.png", dpi=150)
print("Saved cheb_zeros_single.png")
