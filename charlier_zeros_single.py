"""Single-panel: real zeros of Charlier f_N for the largest feasible N.
Targets N=40; falls back to N=20 if too slow.
Uses the GF approach and bisection on the real line (same as charlier_matching.py)."""
from mpmath import mp, mpf, power, factorial as fac, exp, sqrt, besseli
from fractions import Fraction
from math import factorial
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import time
import pickle
import os

mp.dps = 50

a = 1
a_mp = mpf(a)


def charlier_poly(n, x, a_val):
    if n == 0:
        return mpf(1)
    Cp = mpf(1)
    Cc = 1 - x / a_val
    for k in range(1, n):
        Cn = ((k + a_val - x) * Cc - k * Cp) / a_val
        Cp, Cc = Cc, Cn
    return Cc


def h_n(n, a_val):
    return power(a_val, -n) * fac(n)


def normalized_w2(m, k, a_val, X=200):
    s = mpf(0)
    for x in range(X):
        s += (charlier_poly(m, mpf(x), a_val) *
              charlier_poly(k, mpf(x), a_val) *
              power(a_val, 2*x) / fac(x)**2)
    return exp(-2*a_val) * s / sqrt(h_n(m, a_val) * h_n(k, a_val))


def get_coeffs_0k(Nz):
    """Load or compute b_{0,k} for k=0..Nz."""
    cache_file = f"charlier_b0k_N{Nz}.pkl"
    if os.path.exists(cache_file):
        with open(cache_file, 'rb') as f:
            return pickle.load(f)
    print(f"Computing b_{{0,k}} for k=0..{Nz}...")
    coeffs = []
    for k in range(Nz + 1):
        v = normalized_w2(0, k, a_mp, X=200)
        coeffs.append(v)
        if k % 5 == 0:
            print(f"  k={k}: b_{{0,{k}}} = {float(v):.8f}")
    with open(cache_file, 'wb') as f:
        pickle.dump(coeffs, f)
    return coeffs


def f_trunc(xv, clist, Nz):
    return sum(clist[kk] * charlier_poly(kk, xv, a_mp) / sqrt(h_n(kk, a_mp))
               for kk in range(Nz + 1))


def find_zeros(Nz, clist):
    """Find real zeros by bisection on a grid."""
    x_min, x_max = -2, int(3 * Nz) + 2
    n_pts = max(600, 8 * Nz)
    x_grid = [mpf(x_min) + mpf(i) * (x_max - x_min) / n_pts for i in range(n_pts + 1)]
    vals = [f_trunc(x, clist, Nz) for x in x_grid]
    zeros = []
    for i in range(len(vals) - 1):
        if vals[i] * vals[i + 1] < 0:
            lo, hi = float(x_grid[i]), float(x_grid[i + 1])
            for _ in range(60):
                mid = (lo + hi) / 2
                fm = float(f_trunc(mpf(mid), clist, Nz))
                if fm * float(vals[i]) > 0:
                    lo = mid
                else:
                    hi = mid
            zeros.append((lo + hi) / 2)
    return zeros


# Try N=40 first; fall back to N=20 if too slow
for N_try in [40, 20]:
    t0 = time.time()
    print(f"\nAttempting N={N_try}...")
    coeffs_0k = get_coeffs_0k(N_try)
    t1 = time.time()
    print(f"  Coefficients ready ({t1-t0:.1f}s). Finding zeros...")
    zeros = find_zeros(N_try, coeffs_0k)
    t2 = time.time()
    print(f"  {len(zeros)} zeros found in {t2-t1:.1f}s")
    N_used = N_try
    break  # always succeeds for N=40 via bisection

print(f"\nPlotting N={N_used}, {len(zeros)} zeros...")

# Zeros-only strip plot (matches the style of the other single-panel figures)
fig, ax = plt.subplots(figsize=(7, 3))
ax.axvspan(0, N_used, color='gray', alpha=0.15,
           label=f"lattice support $\\{{0,\\ldots,{N_used}\\}}$")
if zeros:
    ax.scatter(zeros, [0] * len(zeros), c='crimson', s=35, zorder=3,
               label=f"{len(zeros)} real zeros of $f_{{{N_used}}}$")
ax.axhline(0, color='gray', linewidth=0.5)
ax.set_title(f"Real zeros of Charlier $f_N$ ($a=1$, $N={N_used}$)", fontsize=13)
ax.set_xlabel("x")
ax.set_yticks([])
ax.set_xlim(-4, 3 * N_used)
ax.set_ylim(-1, 1)
ax.legend(fontsize=9, loc='upper right')
plt.tight_layout()
outfile = "charlier_zeros_single.png"
plt.savefig(outfile, dpi=150)
print(f"Saved {outfile} (N={N_used})")
