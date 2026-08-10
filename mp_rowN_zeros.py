"""Row-N zeros for MP lambda=1/4 family.

p_N(x) = sum_{k=0}^{N} c_{N,k} * Phat_k(x)

where c_{N,k} = b_{N,k} / sqrt(nrm_N * nrm_k)
            = int P_N(x) P_k(x) W(x)^{3/2} dx / (nrm_N * nrm_k)

and W(x) = |Gamma(1/4 + ix)|^2, nrm_n = int P_n^2 W dx.

Since we want p_N as a polynomial, and only the zeros matter,
we can absorb 1/nrm_N into overall factor (k-independent), giving:
  p_N(x) proportional to sum_{k=0}^N [int P_N P_k W^{3/2} dx / nrm_k] * P_k(x)

We compute the row coefficients b_{N,k}^{eff} = int P_N(x) P_k(x) W^{3/2} dx / nrm_k
for k = 0..N via mpmath quadrature.

MP recurrence (phi=pi/2): (n+1) P_{n+1} = 2x P_n - (n - 1/2) P_{n-1}
P_0 = 1, P_1 = 2x.

Note: Only even/odd parity terms survive by symmetry (W^{3/2} is even, P_n has parity (-1)^n).
So b_{N,k} = 0 if N+k is odd. We only include terms with N+k even.

Target: N=20. Falls back to N=12 if timing issues.
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
CACHE_FILE = f"mp_rowN_{N_TARGET}_cache.pkl"


def W(x):
    g = mp.gamma(mp.mpf(1) / 4 + 1j * x)
    return (g * mp.conj(g)).real


def quadW_half(f):
    """2 * int_0^inf f(x) dx for even integrand."""
    return 2 * mp.quad(f, NODES)


def mp_poly_coeffs(nmax):
    """Monomial Fraction coefficients of P_0..P_nmax.
    Recurrence: (n+1) P_{n+1} = 2x P_n - (n - 1/2) P_{n-1}.
    """
    polys = [[Fraction(1)], [Fraction(0), Fraction(2)]]
    for n in range(1, nmax):
        pn, pm = polys[n], polys[n - 1]
        nxt = [Fraction(0)] * (n + 2)
        for j, c in enumerate(pn):
            nxt[j + 1] += 2 * c
        for j, c in enumerate(pm):
            nxt[j] -= (Fraction(n) - Fraction(1, 2)) * c
        polys.append([x / (n + 1) for x in nxt])
    return polys


def evalP(n, x, polys=None):
    """Evaluate P_n(x) via recurrence."""
    p0, p1 = mp.mpf(1), 2 * x
    if n == 0:
        return p0
    if n == 1:
        return p1
    for j in range(1, n):
        p0, p1 = p1, (2 * x * p1 - (j - mp.mpf(1) / 2) * p0) / (j + 1)
    return p1


print(f"Computing norms nrm_k for k=0..{N_TARGET}...")
CACHE_NRMS = f"mp_nrms_{N_TARGET}.pkl"
if os.path.exists(CACHE_NRMS):
    with open(CACHE_NRMS, 'rb') as f:
        nrm = pickle.load(f)
    print(f"  Loaded {len(nrm)} norms from cache")
else:
    nrm = []
    for k in range(N_TARGET + 1):
        v = quadW_half(lambda x, k=k: evalP(k, x) ** 2 * W(x))
        nrm.append(v)
        if k % 5 == 0:
            print(f"  nrm[{k}] = {mp.nstr(v, 8)}")
    with open(CACHE_NRMS, 'wb') as f:
        pickle.dump(nrm, f)
    print(f"  Saved {CACHE_NRMS}")

print(f"\nComputing row-{N_TARGET} coefficients b_{{N,k}} for k=0..{N_TARGET}...")
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, 'rb') as f:
        b_row = pickle.load(f)
    print(f"  Loaded {len(b_row)} coefficients from cache")
else:
    b_row = []
    for k in range(N_TARGET + 1):
        if (N_TARGET + k) % 2 == 1:
            b_row.append(mp.mpf(0))
            continue
        w32 = mp.mpf(3) / 2
        v = quadW_half(lambda x, k=k: evalP(N_TARGET, x) * evalP(k, x) * mp.power(W(x), w32))
        b_row.append(v / nrm[k])  # effective: int P_N P_k W^{3/2} / nrm_k
        if k % 5 == 0:
            print(f"  b_{{N,{k}}} = {mp.nstr(b_row[-1], 8)}")
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(b_row, f)
    print(f"  Saved {CACHE_FILE}")

# Coefficient polynomial: g_N(x) = sum_k bhat_{N,k} x^k with the TRUE
# normalized entries bhat_{N,k} = int Phat_N Phat_k W^{3/2}.
# Cached b_row[k] = int P_N P_k W^{3/2} / nrm_k, so up to the
# k-independent factor 1/sqrt(nrm_N):
#   bhat_{N,k} propto b_row[k] * sqrt(nrm_k)
print(f"\nComputing zeros of g_{N_TARGET}(x) = sum_k b_{{N,k}} x^k ...")
coeffs = [b_row[k] * mp.sqrt(nrm[k]) for k in range(N_TARGET + 1)]
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
plt.savefig("mp_rowN_zeros.png", dpi=150)
print("Saved mp_rowN_zeros.png")
