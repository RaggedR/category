"""
Zeros of truncated Chebyshev-T matching expansions.
A_N(x) = ∑_{k=0}^{N} b_{0k} P̂_k^T(x)
       = ∑_{k=0}^{N} c_k T_k(x)  in Chebyshev basis

where c_0 = b_{00}/√π, c_k = b_{0k}√(2/π) for k≥1.

We use the closed-form b_{0k} for arbitrarily large N (no truncation at 12).
Then plot zeros in ℂ and compare to Bernstein ellipse with parameter ρ.
"""

import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import pickle

mp.mp.dps = 50

# ─── Reconstruct b_{0k} from closed form for large N ───────────────────────
def I_sin_cos(p, q):
    """I(p,q) = ∫₀^π sin^p θ cos(qθ) dθ
             = cos(πq/2) · π Γ(p+1) / (2^p Γ((p+q+2)/2) Γ((p−q+2)/2)).
    The phase factor cos(πq/2) is essential: it kills all odd q and
    contributes (−1)^{q/2} for even q (which cancels the sign of the
    reflected Gamma, making e.g. all I(−1/2, 2n) positive)."""
    p = mp.mpf(p)
    q = int(abs(q))
    if q % 2 == 1:
        return mp.mpf(0)
    phase = (-1) ** (q // 2)
    num = phase * mp.pi * mp.gamma(p + 1)
    den = mp.power(2, p) * mp.gamma((p + q + 2)/2) * mp.gamma((p - q + 2)/2)
    return num / den

def b_T_0k(k):
    """b_{0k} for Cheb-T, m=0."""
    p = mp.mpf(-1)/2
    integral = I_sin_cos(p, k)  # since m=0, m+k=k and |m-k|=k
    nm = 1/mp.sqrt(mp.pi)
    nk = 1/mp.sqrt(mp.pi) if k == 0 else mp.sqrt(2/mp.pi)
    return nm * nk * integral

def b_U_0k(k):
    """b_{0k} for Cheb-U, m=0."""
    hU = mp.pi / 2
    p = mp.mpf(3)/2
    I_diff = I_sin_cos(p, k) - I_sin_cos(p, k+2)
    return I_diff / (2 * hU)

def cheb_T_zeros(N):
    """Zeros of A_N^T(x) = ∑_{k=0}^N c_k T_k(x)."""
    # Chebyshev coefficients c_k
    c = []
    for k in range(N+1):
        b = float(b_T_0k(k))
        if k == 0:
            c.append(b / np.sqrt(np.pi))
        else:
            c.append(b * np.sqrt(2/np.pi))
    poly = np.polynomial.chebyshev.Chebyshev(c)
    roots = poly.roots()
    return roots

def cheb_U_zeros(N):
    """Zeros of A_N^U(x) = ∑_{k=0}^N b_{0k} P̂_k^U(x)."""
    # P̂_k^U = U_k / sqrt(π/2)
    # A_N^U(x) = ∑ b_{0k}/√(π/2) · U_k(x)
    # U_k(cosθ) = sin((k+1)θ)/sinθ
    # We need polynomial coefficients for ∑ d_k U_k(x).
    # Convert U-coefficients to Chebyshev via U_k = 2T_... or direct.
    # Easiest: convert to monomial via numpy
    d = []
    for k in range(N+1):
        d.append(float(b_U_0k(k)) / np.sqrt(np.pi/2))
    # U_k coefficients → chebychev u polynomial object
    poly = np.polynomial.chebyshev.Chebyshev(
        np.polynomial.chebyshev.cheb2poly(
            np.polynomial.polynomial.Polynomial(
                np.zeros(N+1)  # placeholder
            ).coef
        )
    )
    # Better: use np.polynomial for Chebyshev-U directly
    # U_k is not in numpy standard, but we know: U_0=1,U_1=2x,U_n=2x U_{n-1}-U_{n-2}
    # Convert U-basis to monomial
    u_mono = [[1.0], [2.0, 0.0]]  # U_0=1 as poly [const], U_1=2x
    # Wait: poly [a,b] means a + b*x in numpy Polynomial
    # U_0 = 1 → [1]
    # U_1 = 2x → [0, 2]
    u_mono = {0: np.array([1.0]), 1: np.array([0.0, 2.0])}
    for k in range(2, N+1):
        # U_k = 2x U_{k-1} - U_{k-2}
        um1 = u_mono[k-1]
        um2 = u_mono[k-2]
        # 2x * um1: shift by 1 and multiply by 2
        xu = np.zeros(len(um1)+1)
        xu[1:] = 2*um1
        # subtract um2
        um2_pad = np.zeros(len(xu))
        um2_pad[:len(um2)] = um2
        u_mono[k] = xu - um2_pad
    # Build total polynomial in monomial basis
    max_len = max(len(u_mono[k]) for k in range(N+1))
    total = np.zeros(max_len)
    for k in range(N+1):
        coeff = d[k]
        mono = u_mono[k]
        total[:len(mono)] += coeff * mono
    # Find roots
    roots = np.roots(total[::-1])  # np.roots expects highest-degree first
    return roots

# ─── Bernstein ellipse ──────────────────────────────────────────────────────
def bernstein_ellipse(rho, npts=500):
    """Bernstein ellipse: z = (ρ e^{iθ} + ρ^{-1}e^{-iθ})/2 = cosh(σ+iθ)
    with σ = log(ρ), so semi-axes a=(ρ+1/ρ)/2, b=(ρ-1/ρ)/2."""
    a = (rho + 1/rho) / 2
    b = (rho - 1/rho) / 2
    theta = np.linspace(0, 2*np.pi, npts)
    x = a * np.cos(theta)
    y = b * np.sin(theta)
    return x, y

# ─── MAIN PLOT ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    Ns = [10, 20, 40]

    # Collect zeros for Cheb-T
    all_zeros_T = {}
    for N in Ns:
        z = cheb_T_zeros(N)
        all_zeros_T[N] = z
        print(f"N={N}: {len(z)} zeros computed")

    # Collect zeros for Cheb-U
    all_zeros_U = {}
    for N in Ns:
        z = cheb_U_zeros(N)
        all_zeros_U[N] = z
        print(f"U N={N}: {len(z)} zeros computed")

    # ─── Figure ─────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    colors = ['tab:blue', 'tab:orange', 'tab:green']
    markers = ['o', 's', '^']

    for ax, zeros_dict, title, label in [
        (axes[0], all_zeros_T, "Chebyshev-T  (α=β=−½)", "T"),
        (axes[1], all_zeros_U, "Chebyshev-U  (α=β=+½)", "U"),
    ]:
        for (N, col, mk) in zip(Ns, colors, markers):
            zs = zeros_dict[N]
            # real zeros in [−1,1] (genuine)
            real_mask = (np.abs(zs.imag) < 1e-10) & (np.abs(zs.real) <= 1.01)
            # spurious: outside or complex
            spur_mask = ~real_mask
            ax.scatter(zs[real_mask].real, zs[real_mask].imag,
                       s=15, color=col, marker=mk, label=f'N={N} genuine', zorder=4)
            ax.scatter(zs[spur_mask].real, zs[spur_mask].imag,
                       s=15, color=col, marker=mk, facecolors='none',
                       label=f'N={N} spurious', zorder=4)

        # Draw Bernstein ellipses for rho = 1.1, 1.5, 2.0
        for rho, ls in [(1.1, ':'), (1.5, '--'), (2.0, '-')]:
            xe, ye = bernstein_ellipse(rho)
            ax.plot(xe, ye, ls, color='gray', alpha=0.5,
                    label=f'Bernstein ρ={rho}' if label == "T" else '_nolegend_')

        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_xlim(-3.5, 3.5)
        ax.set_ylim(-2.5, 2.5)
        ax.set_xlabel('Re(z)', fontsize=12)
        ax.set_ylabel('Im(z)', fontsize=12)
        ax.set_title(title, fontsize=13)
        ax.legend(fontsize=8, ncol=2)
        ax.set_aspect('equal')
        # Mark [−1,1] interval
        ax.axvline(-1, color='k', linestyle=':', linewidth=0.8, alpha=0.5)
        ax.axvline(1, color='k', linestyle=':', linewidth=0.8, alpha=0.5)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Zeros of truncated matching expansions  $A_N(x) = \\sum_{k\\leq N} b_{0k}\\hat{P}_k(x)$\n"
                 "Dashed ellipses: Bernstein ellipses; hollow symbols = spurious zeros",
                 fontsize=12)
    plt.tight_layout()
    plt.savefig('/Users/robin/git/category/jacobi_zeros.png', dpi=150, bbox_inches='tight')
    print("Saved jacobi_zeros.png")
    plt.close()

    # ─── Print zero analysis ─────────────────────────────────────────────────
    print("\n=== ZERO ANALYSIS ===")
    for N in Ns:
        zs = all_zeros_T[N]
        real_in = zs[(np.abs(zs.imag)<1e-10) & (np.abs(zs.real)<=1.01)]
        spur = zs[~((np.abs(zs.imag)<1e-10) & (np.abs(zs.real)<=1.01))]
        print(f"\nCheb-T N={N}:")
        print(f"  {len(real_in)} genuine zeros in [−1,1]")
        print(f"  {len(spur)} spurious zeros")
        if len(spur) > 0:
            mods = np.abs(spur)
            # Bernstein parameter: ρ = |z| + sqrt(z²-1)?
            # For z = (ρe^{iθ} + ρ^{-1}e^{-iθ})/2, |z| ≤ (ρ+1/ρ)/2 = cosh(log ρ)
            # So ρ = |z| + sqrt(|z|²-1) for real z outside [−1,1]
            rhos = []
            for z in spur:
                if abs(z.real) > 1 and abs(z.imag) < 1e-6:
                    rho = abs(z.real) + np.sqrt(abs(z.real)**2 - 1)
                    rhos.append(rho)
            if rhos:
                print(f"  Real spurious Bernstein ρ: {[f'{r:.3f}' for r in rhos]}")
            # Complex spurious: max |z|
            cmx = spur[np.abs(spur.imag) > 1e-6]
            if len(cmx) > 0:
                print(f"  Complex spurious |z| range: [{np.min(np.abs(cmx)):.3f}, {np.max(np.abs(cmx)):.3f}]")

    print("\n=== Expected Bernstein ellipse behavior ===")
    print("For Chebyshev approximation with decay b_{0k} ~ k^{-1/2}:")
    print("The 'balance' approximation: b_{0N} ~ ρ^{-N} for a point on Bernstein ellipse ρ.")
    print("But b_{0k} decays as k^{-1/2}, NOT geometrically → ρ→1⁺ as N→∞.")
    print("This means zeros accumulate on [−1,1] itself (the Bernstein ellipse ρ=1).")
    print("Compare to Hermite: geometric decay → ρ=√3 fixed Bernstein-like radius.")
    print("Compare to MP: power-law decay → zeros also crowd toward real axis.")
