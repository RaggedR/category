"""Row-N zeros for Charlier (a=1) family.

p_N(x) = sum_{k=0}^{N} b_{N,k} / sqrt(h_k) * C_k(x)

where b_{N,k} = e^{-2} * sum_{x=0}^{X} C_N(x; 1) * C_k(x; 1) / (x!)^2

and h_k = a^{-k} k! = k! for a=1.

The Charlier polynomials satisfy the orthogonality
  sum_{x=0}^{inf} C_m(x;a) C_k(x;a) * a^x / x! = e^a * a^{-n} * n! * delta_{mk}
  => h_k = e^1 * 1^{-k} * k! = e * k! (with the e^{-a} weight)
  BUT the matching uses weight w^2 where w(x) = e^{-a/2} a^{x/2} / sqrt(x!),
  i.e. b_{N,k} = e^{-2a} sum_x C_N(x) C_k(x) a^{2x} / (x!)^2 / sqrt(h_N h_k)

From charlier_zeros_single.py:
  normalized_w2(m, k, a) = e^{-2a} sum_{x=0}^X C_m(x) C_k(x) a^{2x} / (x!)^2 / sqrt(h_m h_k)

So b_{N,k} = normalized_w2(N, k, a=1).

The polynomial p_N(x) = sum_{k=0}^N b_{N,k} * C_k(x) / sqrt(h_k)
  = sum_{k=0}^N normalized_w2(N, k, 1) * C_k(x) / sqrt(h_k)

Wait: normalized_w2 already divides by sqrt(h_m h_k), so:
  b_{N,k} (unnormalized) = e^{-2} sum_x C_N(x) C_k(x) / (x!)^2
  b_{N,k} (normalized) = above / sqrt(h_N h_k)

The approximant is sum_k b_{N,k}^{norm} * Chat_k where Chat_k = C_k / sqrt(h_k).
This gives sum_k [b_{N,k}^{norm} / sqrt(h_k)] * C_k
         = sum_k [e^{-2} sum_x C_N(x) C_k(x) / (x!)^2 / (h_N h_k)] * C_k
which is what we want (up to the 1/h_N factor which is k-independent).

For zeros we can drop 1/h_N:
  p_N(x) proportional to sum_{k=0}^N [e^{-2} sum_x C_N(x') C_k(x') / (x'!)^2 / h_k] * C_k(x)

With a=1, h_k = k!.

Zeros plotted vs zeros of C_N(x;1) on real line (positive real axis, discrete).
Since Charlier zeros are on the real positive axis, we show them as a rug plot.
"""
import mpmath as mp
from mpmath import mpf, power, factorial as fac, exp, sqrt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pickle
import os

mp.dps = 60
a = 1
X_sum = 200  # truncation for the infinite sum over x


def charlier_poly(n, x, a_val=1):
    """Charlier C_n(x; a) via 3-term recurrence."""
    if n == 0:
        return mpf(1)
    p, c = mpf(1), 1 - x / mpf(a_val)
    for k in range(1, n):
        c, p = ((k + a_val - x) * c - k * p) / a_val, c
    return c


def h_n(n):
    """h_n = a^{-n} n! = n! for a=1."""
    return fac(n)


def compute_b_Nk_normalized(N, k, X=X_sum):
    """b_{N,k}^{norm} = e^{-2} sum_{x=0}^X C_N(x) C_k(x) / (x!)^2 / sqrt(h_N h_k)."""
    s = mpf(0)
    for x in range(X):
        s += (charlier_poly(N, mpf(x)) *
              charlier_poly(k, mpf(x)) /
              fac(x) ** 2)
    return exp(-2) * s / sqrt(h_n(N) * h_n(k))


# Cache for row N coefficients
CACHE_FILE = f"charlier_rowN_N{40}.pkl"

N = 40

if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, 'rb') as f:
        b_norm = pickle.load(f)
    print(f"Loaded {len(b_norm)} cached b_{{N,k}} values")
else:
    print(f"Computing b_{{{N},k}} for k=0..{N}...")
    b_norm = []
    for k in range(N + 1):
        v = compute_b_Nk_normalized(N, k, X=X_sum)
        b_norm.append(v)
        if k % 5 == 0:
            print(f"  k={k}: b = {mp.nstr(v, 8)}")
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(b_norm, f)
    print("Saved cache")


def charlier_monomial_coeffs(k, x_max=None):
    """Monomial coefficients of C_k(x) = sum_j c_j x^j.
    Uses the fact that C_k(x) satisfies 3-term recurrence, and accumulate
    polynomial coefficients symbolically.

    Actually, use mpmath polynomial multiplication approach:
    C_0 = 1
    C_1 = 1 - x  (a=1: C_1 = 1 - x/1 = 1 - x)
    C_{n+1} = ((n + a - x) C_n - n C_{n-1}) / a
            = ((n+1-x) C_n - n C_{n-1})   for a=1

    Represent as list of mpf: [c_0, c_1, ..., c_k] for coefficient of x^j.
    """
    if k == 0:
        return [mpf(1)]
    if k == 1:
        return [mpf(1), mpf(-1)]  # 1 - x
    p = [mpf(1)]         # C_0
    c = [mpf(1), mpf(-1)]  # C_1
    for n in range(1, k):
        # C_{n+1} = ((n+1-x)*C_n - n*C_{n-1}) / 1
        # = (n+1)*C_n - x*C_n - n*C_{n-1}
        new_len = len(c) + 1
        nxt = [mpf(0)] * new_len
        # (n+1)*C_n:
        for j, cj in enumerate(c):
            nxt[j] += (n + 1) * cj
        # -x*C_n: shift by 1
        for j, cj in enumerate(c):
            nxt[j + 1] -= cj
        # -n*C_{n-1}:
        for j, cj in enumerate(p):
            nxt[j] -= n * cj
        p, c = c, nxt
    return c


# Coefficient polynomial: g_N(x) = sum_k b_{N,k} x^k with the TRUE
# normalized entries (b_norm is already e^{-2} sum C_N C_k /(x!)^2 / sqrt(h_N h_k)).
print(f"Computing zeros of g_{N}(x) = sum_k b_{{N,k}} x^k ...")
coeffs = list(b_norm)
while len(coeffs) > 1 and coeffs[-1] == 0:
    coeffs.pop()
mp.dps = 80
zs = mp.polyroots(coeffs[::-1], maxsteps=300, extraprec=400)
print(f"  {len(zs)} zeros found")

xs = [float(z.real) for z in zs]
ys = [float(z.imag) for z in zs]

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(xs, ys, s=18, color="royalblue", zorder=3,
           label=fr"zeros of $g_{{{N}}}(x)=\sum_k b_{{{N},k}}x^k$")
ax.axhline(0, color="gray", lw=0.5)
ax.axvline(0, color="gray", lw=0.5)
ax.set_xlabel("Re(z)")
ax.set_ylabel("Im(z)")
# legend removed
plt.tight_layout()
plt.savefig("charlier_rowN_zeros.png", dpi=150)
print("Saved charlier_rowN_zeros.png")
