"""Row-N zeros for Chebyshev-T family.

p_N(x) = sum_{k=0}^{N} b_{N,k} * That_k(x)

where b_{N,k} = [I(N-k) + I(N+k)] / (2 * sqrt(h_N * h_k))

and I(q) = pi * Gamma(1/2) * 2^{1/2} / (Gamma(q/2 + 3/4) * Gamma(3/4 - q/2))
         = integral_{-1}^{1} (1-x^2)^{-1/2} T_m(x) T_k(x) (1-x^2)^{1/2} dx
           (matching integral with weight W = (1-x^2)^{1/2} for T)

Wait: the actual matching integral is
  b_{N,k} = int_{-1}^{1} That_N(x) That_k(x) W^{3/2}(x) dx

where W(x) = (1-x^2)^{-1/2} is the Chebyshev-T weight, so W^{3/2} = (1-x^2)^{-3/4}.

From cheb_zeros_single.py, the formula used is:
  I_sin_cos(p, q) = pi * Gamma(p+1) / (2^p * Gamma((p+q+2)/2) * Gamma((p-q+2)/2))
  b_T_0k = (1/sqrt(pi)) * (1/sqrt(pi) if k==0 else sqrt(2/pi)) * I_sin_cos(-1/2, k)

The general b_{N,k}: by parity, b_{N,k} = 0 if N+k odd.
For general m,k: b_{N,k} = (1/sqrt(h_N * h_k)) * [I(N-k) + I(N+k)] / 2
where I(q) = I_sin_cos(-1/2, q) = pi * Gamma(1/2) / (2^{-1/2} * Gamma((q+1)/2) * Gamma((1-q+1)/2))

Actually re-deriving from first principles:
  b_{N,k} = int_{-1}^1 That_N That_k (1-x^2)^{-3/4} dx
           = (h_N h_k)^{-1/2} int T_N T_k (1-x^2)^{-3/4} dx

Using T_N(x) = cos(N arccos x), substitute x = cos(theta):
  int_{-1}^1 T_N T_k (1-x^2)^{-3/4} dx = int_0^pi cos(N theta) cos(k theta) sin^{-1/2}(theta) d(theta)
  = (1/2)[int_0^pi cos((N-k)theta) sin^{-1/2}(theta) dtheta + int_0^pi cos((N+k)theta) sin^{-1/2}(theta) dtheta]

and int_0^pi cos(q theta) sin^{-1/2}(theta) dtheta = I_sin_cos(-1/2, q) / ... let's use mpmath directly.

h_0 = pi, h_k = pi/2 for k >= 1 (unnormalized).
That_0 = 1/sqrt(pi), That_k = sqrt(2/pi) T_k for k>=1.

So b_{N,k} = (1/sqrt(h_N) * 1/sqrt(h_k)) * (1/2) * [J(N-k) + J(N+k)]
where J(q) = int_0^pi cos(q theta) (sin theta)^{-1/2} dtheta.

J(q) = pi * Gamma(1/2) / (2^{-1/2} * Gamma((|q|+1)/2) * ... )
Let's compute J(q) via mpmath:
  J(q) = Re[ int_0^pi e^{i q theta} sin^{-1/2}(theta) dtheta ]
Note: J(q) = J(-q) since cos is even.
The formula: J(q) = 2 * int_0^{pi/2} cos(2 q phi) ... no, let me just use mpmath.

Actually from the existing script's formula:
  I_sin_cos(p, q) = pi * Gamma(p+1) / (2^p * Gamma((p+q+2)/2) * Gamma((p-q+2)/2))
with p = -1/2:
  J(q) = I_sin_cos(-1/2, |q|) = pi * Gamma(1/2) / (2^{-1/2} * Gamma((|q|+1)/2) * Gamma((1-|q|)/2))

But this has singularities when q is a positive even integer... let's just use mpmath.quad.

For speed at N=40 with ~40 evaluations of J: direct mpmath.quad integration is fine.
"""
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 80

N = 40


def J_integral(q):
    """J(q) = int_0^pi cos(q*theta) * sin(theta)^{-1/2} dtheta.
    This is even in q, so J(q) = J(-q). Uses mpmath.quad.
    """
    q = abs(q)
    if q == 0:
        # int_0^pi sin^{-1/2}(theta) dtheta = sqrt(pi)*Gamma(1/4)/Gamma(3/4) -- but let's compute
        return mp.quad(lambda t: mp.power(mp.sin(t), mp.mpf('-1/2')),
                       [0, mp.pi/4, mp.pi/2, 3*mp.pi/4, mp.pi])
    return mp.quad(lambda t: mp.cos(q * t) * mp.power(mp.sin(t), mp.mpf('-1/2')),
                   [0, mp.pi/4, mp.pi/2, 3*mp.pi/4, mp.pi])


# Precompute h values
h_N = mp.pi if N == 0 else mp.pi / 2

# Build b_{N,k} for k=0..N
# b_{N,k} = (1/sqrt(h_N * h_k)) * (1/2) * [J(N-k) + J(N+k)]
print(f"Precomputing J integrals for row N={N}...")
# We need J(N-k) for k=0..N: i.e. J(N), J(N-1), ..., J(0), J(1), ..., J(N)
# So J(0), J(1), ..., J(2N)
J_cache = {}
for q in range(0, 2 * N + 1):
    J_cache[q] = J_integral(q)
    if q % 10 == 0:
        print(f"  J({q}) = {float(J_cache[q]):.6f}")

print("Computing b_{N,k} values...")
b_values = []
for k in range(N + 1):
    h_k = mp.pi if k == 0 else mp.pi / 2
    q1 = abs(N - k)
    q2 = N + k
    b_Nk = (J_cache[q1] + J_cache[q2]) / (2 * mp.sqrt(h_N * h_k))
    b_values.append(b_Nk)

# Verify b_{N,k}=0 for N+k odd (parity check)
# T_N has parity (-1)^N, so b_{N,k}=0 when N+k odd (odd integrand on symmetric domain)
n_nonzero = sum(1 for k, b in enumerate(b_values) if abs(b) > 1e-50)
print(f"  Nonzero b_{{N,k}}: {n_nonzero}")


def chebT_monomial_coeffs(k, dps=80):
    """Monomial coefficients of T_k(x) as mpmath numbers."""
    # T_k(x) = 2^{k-1} prod_{j: 0<j<k, j step 2} (x - cos(...)): harder to expand
    # Better: use the recurrence T_0=1, T_1=x, T_{n+1}=2x T_n - T_{n-1}
    # and accumulate monomial coefficients.
    if k == 0:
        return [mp.mpf(1)]
    if k == 1:
        return [mp.mpf(0), mp.mpf(1)]
    p, c = [mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]
    for _ in range(1, k):
        nxt = [mp.mpf(0)] * (len(c) + 1)
        for j, cj in enumerate(c):
            nxt[j + 1] += 2 * cj
        for j, cj in enumerate(p):
            nxt[j] -= cj
        p, c = c, nxt
    return c


print(f"Computing zeros of g_{N}(x) = sum_k b_{{N,k}} x^k ...")
# The b_values ARE the true normalized b_{N,k}; use them directly as
# monomial coefficients.
coeffs = list(b_values)
while len(coeffs) > 1 and coeffs[-1] == 0:
    coeffs.pop()
mp.mp.dps = 120
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
plt.savefig("cheb_rowN_zeros.png", dpi=150)
print("Saved cheb_rowN_zeros.png")
