"""Row-N zeros for Laguerre (alpha=0) family.

The b_{N,k} coefficients come from the GF 2/(3-w-u-wu):
  b_{N,k} = (2/3) * sum_{c=0}^{min(N,k)} C(N+k-c, N-c, k-c) * 3^c / 3^{N+k}
           = (2/3) * sum_c C(N,c) * C(N+k-c, N) * 3^c / 3^{N+k}

The orthonormal basis is L̂_k = L_k / sqrt(h_k) where h_k = 1 for alpha=0
(Laguerre is already orthonormal on [0,inf) w.r.t. e^{-x} dx when h_n = 1).

Wait: for alpha=0, the Laguerre poly L_k has norm h_k = integral_0^inf L_k(x)^2 e^{-x} dx = 1
(they are orthonormal already since L_k = L_k^(0) and <L_m,L_n>_{e^-x} = delta_{mn}).

So p_N(x) = sum_{k=0}^{N} b_{N,k} * L_k(x).

The weight factor W^{3/2} = e^{-3x/2} on [0,inf), and the matching matrix is
  b_{N,k} = int_0^inf L_N(x) L_k(x) e^{-3x/2} dx  (unnormalized)

From the GF formula (laguerre_zeros_single.py):
  E_{N,k} = sum_{c=0}^{min(N,k)} C(N,c)*C(N+k-c,N) * 3^c
  b_{N,k} = (2/3) * E_{N,k} / 3^{N+k}

where L_k monomial: L_k(x) = sum_j C(k,j)*(-1)^j/j! * x^j.

Zeros of p_N (crimson) vs zeros of L_N (open circles).
Uses exact Fraction arithmetic + mpmath polyroots at dps=120.
"""
from fractions import Fraction
from math import comb, factorial
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 120

N = 40


def E_Nk(N, k):
    """Exact integer E_{N,k} = sum_c C(N,c)*C(N+k-c,N)*3^c."""
    tot = 0
    for c in range(min(N, k) + 1):
        tot += comb(N, c) * comb(N + k - c, N) * (3 ** c)
    return tot


def b_Nk(N, k):
    """Exact Fraction b_{N,k} = 2/3 * E_{N,k} / 3^{N+k}."""
    e = E_Nk(N, k)
    return Fraction(2 * e, 3 ** (N + k + 1))


def laguerre_monomial_coeffs(k):
    """Monomial coefficients of L_k(x) = sum_j C(k,j)*(-1)^j/j! * x^j.
    Returns list of Fraction coeffs[j] for j=0..k."""
    return [Fraction(comb(k, j) * ((-1) ** j), factorial(j)) for j in range(k + 1)]


def rowN_polynomial_coeffs(N):
    """Build p_N(x) = sum_{k=0}^N b_{N,k} * L_k(x) in mpmath."""
    poly = [mp.mpf(0)] * (N + 1)
    for k in range(N + 1):
        bNk = b_Nk(N, k)
        if bNk == 0:
            continue
        bNk_mp = mp.mpf(bNk.numerator) / mp.mpf(bNk.denominator)
        lk = laguerre_monomial_coeffs(k)
        for j, cj in enumerate(lk):
            if cj != 0:
                poly[j] += bNk_mp * mp.mpf(cj.numerator) / mp.mpf(cj.denominator)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def laguerre_N_zeros(N):
    """Zeros of L_N via polyroots."""
    lN = laguerre_monomial_coeffs(N)
    cs = [mp.mpf(c.numerator) / mp.mpf(c.denominator) for c in lN[::-1]]
    return mp.polyroots(cs, maxsteps=200, extraprec=300)


print(f"Building coefficient polynomial g_{N}(x) = sum_k b_{{N,k}} x^k ...")
# b_{N,k} = (2/3) E_{N,k} / 3^{N+k}; up to a k-independent constant:
#   c_k = E_{N,k} / 3^k
coeffs = [mp.mpf(E_Nk(N, k)) / mp.power(3, k) for k in range(N + 1)]
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
plt.savefig("laguerre_rowN_zeros.png", dpi=150)
print("Saved laguerre_rowN_zeros.png")
