"""Row-N zeros for MP lambda=3/4 family.

p_N(x) = sum_{k=0}^{N} c_{N,k} * Qhat_k(x)

where the MP polynomials Q_n satisfy:
  (n+1) Q_{n+1} = 2x Q_n - (n + 1/2) Q_{n-1}
  Q_0 = 1, Q_1 = 2x

(This is the lambda=3/4 version: the recurrence coefficient changes to n+2lambda-1 = n+1/2.)

Weight W(x) = |Gamma(3/4 + ix)|^2.
Matching: c_{N,k} = int Q_N Q_k W^{3/2} dx / (nrm34_N * nrm34_k).

Since zeros are invariant under overall (k-independent) factor 1/nrm34_N:
  eff. coeff: int Q_N Q_k W^{3/2} dx / nrm34_k.

Target: N=20.
"""
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fractions import Fraction
import pickle
import os

mp.mp.dps = 80

NODES = list(range(0, 21, 2)) + [25, 30, 40, 50]
N_TARGET = 20
CACHE_FILE = f"mp34_rowN_{N_TARGET}_cache.pkl"
CACHE_NRMS = f"mp34_nrms_{N_TARGET}.pkl"


def W34(x):
    """Weight for lambda=3/4."""
    g = mp.gamma(mp.mpf(3) / 4 + 1j * x)
    return (g * mp.conj(g)).real


def quadW_half(f):
    """2 * int_0^inf f(x) dx for even integrand."""
    return 2 * mp.quad(f, NODES)


def mp34_poly_coeffs(nmax):
    """Monomial Fraction coefficients of Q_0..Q_nmax (lambda=3/4).
    Recurrence: (n+1) Q_{n+1} = 2x Q_n - (n + 1/2) Q_{n-1}.
    """
    polys = [[Fraction(1)], [Fraction(0), Fraction(2)]]
    for n in range(1, nmax):
        pn, pm = polys[n], polys[n - 1]
        nxt = [Fraction(0)] * (n + 2)
        for j, c in enumerate(pn):
            nxt[j + 1] += 2 * c
        for j, c in enumerate(pm):
            # coefficient is -(n + 2*lambda - 1) = -(n + 1/2) for lambda=3/4
            nxt[j] -= (Fraction(n) + Fraction(1, 2)) * c
        polys.append([x / (n + 1) for x in nxt])
    return polys


def evalQ(n, x):
    """Evaluate Q_n(x) via recurrence (lambda=3/4)."""
    p0, p1 = mp.mpf(1), 2 * x
    if n == 0:
        return p0
    if n == 1:
        return p1
    for j in range(1, n):
        p0, p1 = p1, (2 * x * p1 - (j + mp.mpf(1) / 2) * p0) / (j + 1)
    return p1


print(f"Computing norms nrm34_k for k=0..{N_TARGET}...")
if os.path.exists(CACHE_NRMS):
    with open(CACHE_NRMS, 'rb') as f:
        nrm34 = pickle.load(f)
    print(f"  Loaded {len(nrm34)} norms from cache")
else:
    nrm34 = []
    for k in range(N_TARGET + 1):
        v = quadW_half(lambda x, k=k: evalQ(k, x) ** 2 * W34(x))
        nrm34.append(v)
        if k % 5 == 0:
            print(f"  nrm34[{k}] = {mp.nstr(v, 8)}")
    with open(CACHE_NRMS, 'wb') as f:
        pickle.dump(nrm34, f)
    print(f"  Saved {CACHE_NRMS}")

print(f"\nComputing row-{N_TARGET} coefficients for lambda=3/4...")
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, 'rb') as f:
        b_row = pickle.load(f)
    print(f"  Loaded {len(b_row)} coefficients from cache")
else:
    w32 = mp.mpf(3) / 2
    b_row = []
    for k in range(N_TARGET + 1):
        # Parity: Q_n has parity (-1)^n, W34 is even, so b_{N,k}=0 if N+k odd
        if (N_TARGET + k) % 2 == 1:
            b_row.append(mp.mpf(0))
            continue
        v = quadW_half(lambda x, k=k: evalQ(N_TARGET, x) * evalQ(k, x) * mp.power(W34(x), w32))
        b_row.append(v / nrm34[k])
        if k % 5 == 0:
            print(f"  b_{{{N_TARGET},{k}}} = {mp.nstr(b_row[-1], 8)}")
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(b_row, f)
    print(f"  Saved {CACHE_FILE}")

# Coefficient polynomial: g_N(x) = sum_k bhat_{N,k} x^k with the TRUE
# normalized entries bhat_{N,k} = int Qhat_N Qhat_k W^{3/2}.
# Cached b_row[k] = int Q_N Q_k W^{3/2} / nrm34_k, so up to the
# k-independent factor 1/sqrt(nrm34_N):
#   bhat_{N,k} propto b_row[k] * sqrt(nrm34_k)
print(f"\nComputing zeros of g_{N_TARGET}(x) = sum_k b_{{N,k}} x^k ...")
coeffs = [b_row[k] * mp.sqrt(nrm34[k]) for k in range(N_TARGET + 1)]
while len(coeffs) > 1 and abs(coeffs[-1]) < 1e-100:
    coeffs.pop()
mp.mp.dps = 120
zs = mp.polyroots(coeffs[::-1], maxsteps=300, extraprec=400)
print(f"  {len(zs)} zeros found")

xs = [float(z.real) for z in zs]
ys = [float(z.imag) for z in zs]

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(xs, ys, s=18, color="royalblue", zorder=3,
           label=fr"zeros of $g_{{{N_TARGET}}}(x)=\sum_k b_{{{N_TARGET},k}}x^k$")
ax.axhline(0, color="gray", lw=0.5)
ax.axvline(0, color="gray", lw=0.5)
ax.set_xlabel("Re(x)")
ax.set_ylabel("Im(x)")
# legend removed
plt.tight_layout()
plt.savefig("mp34_rowN_zeros.png", dpi=150)
print("Saved mp34_rowN_zeros.png")
