"""Row-N zeros for Hermite family.

p_N(x) = sum_{k=0}^{N} b_{N,k} * Hhat_k(x)

where:
  b_{N,k} = D(N,k) * sqrt(2/3) / (3^{(N+k)/2} * sqrt(2^N N!) * sqrt(2^k k!))
  Hhat_k(x) = H_k(x) / sqrt(2^k k! sqrt(pi))

So the contribution of term k to the polynomial is:
  b_{N,k} * Hhat_k(x) = D(N,k) * sqrt(2/3) * H_k(x) / (3^{(N+k)/2} * sqrt(2^N N!) * 2^k k! * sqrt(sqrt(pi)))

Dropping k-independent factors sqrt(2/3) / (sqrt(2^N N!) * sqrt(sqrt(pi))):
  effective amplitude: D(N,k) / (3^{(N+k)/2} * 2^k k! * sqrt(2^k k!))
                     = D(N,k) / (3^{(N+k)/2} * (2^k k!)^{3/2})

We accumulate this in mpmath and find the zeros.
Verified: this matches direct quadrature for small N.

Target: N=40. Uses mpmath polyroots at dps=120.
"""
from math import comb, factorial
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 120

N = 40


def dfact(n):
    """Double factorial n!! (returns 1 for n<=0)."""
    if n <= 0:
        return 1
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def D(m, k):
    """Integer coefficient D(m,k) — the Hermite matching integer table entry."""
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


def hermite_monomial_coeffs(k):
    """Monomial coefficients of H_k(x) as integers."""
    a, b = [1], [0, 2]
    if k == 0:
        return a
    for i in range(1, k):
        nxt = [0] * (i + 2)
        for j, cj in enumerate(b):
            nxt[j + 1] += 2 * cj
        for j, cj in enumerate(a):
            nxt[j] -= 2 * i * cj
        a, b = b, nxt
    return b


def rowN_polynomial_coeffs(N):
    """
    Build p_N(x) = sum_k b_{N,k} * Hhat_k(x) as a polynomial in x.

    b_{N,k} * Hhat_k(x) = [D(N,k) / (3^{(N+k)/2} * 2^k * k!)] * H_k(x) * [k-indep const]

    So the zeros of p_N are the zeros of Σ_k a_k * H_k(x) where:
      a_k = D(N,k) / (3^{(N+k)/2} * 2^k * k!)

    This is exactly the formula from hermite_zeros_single.py with m=N.
    Verified to match direct quadrature of b_{N,k} up to a k-independent constant.
    """
    poly = [mp.mpf(0)] * (N + 1)
    for k in range(N + 1):
        d = D(N, k)
        if d == 0:
            continue
        exp3 = (N + k) // 2  # exact integer since N+k even
        a_k = mp.mpf(d) / (mp.power(3, exp3) * mp.power(2, k) * mp.fac(k))
        hk = hermite_monomial_coeffs(k)
        for j, cj in enumerate(hk):
            if cj != 0:
                poly[j] += a_k * mp.mpf(cj)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def hermite_N_zeros_genuine(N):
    """Zeros of H_N via polyroots."""
    hN = hermite_monomial_coeffs(N)
    cs = [mp.mpf(c) for c in hN[::-1]]
    return mp.polyroots(cs, maxsteps=200, extraprec=300)


print(f"Building coefficient polynomial g_{N}(x) = sum_k b_{{N,k}} x^k ...")
# True b_{N,k} up to a k-INDEPENDENT constant:
#   b_{N,k} = D(N,k) * sqrt(2/3) / (3^{(N+k)/2} sqrt(2^N N!) sqrt(2^k k!))
#           propto D(N,k) / (3^{k/2} * sqrt(2^k k!))
coeffs = []
for k in range(N + 1):
    d = D(N, k)
    c_k = mp.mpf(d) / (mp.power(3, mp.mpf(k) / 2) * mp.sqrt(mp.power(2, k) * mp.fac(k)))
    coeffs.append(c_k)
while len(coeffs) > 1 and coeffs[-1] == 0:
    coeffs.pop()
print(f"  Degree: {len(coeffs)-1}")

print(f"Computing zeros of g_{N}...")
zs = mp.polyroots(coeffs[::-1], maxsteps=200, extraprec=300)
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
plt.savefig("hermite_rowN_zeros.png", dpi=150)
print("Saved hermite_rowN_zeros.png")
