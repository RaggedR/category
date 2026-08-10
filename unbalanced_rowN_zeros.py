"""Row-N zeros for unbalanced family.

The unbalanced weight is W(x) = Gamma(1/4 + ix)^2 * Gamma(3/4 - ix).
The unbalanced polynomials are Q̃_n(x) = P^{(1/2)}_n(x + i/4),
i.e. the standard MP lambda=1/2 polynomials evaluated at x + i/4.

MP lambda=1/2 polynomials R_n satisfy:
  (n+1) R_{n+1} = 2x R_n - n R_{n-1}
  R_0 = 1, R_1 = 2x
(lambda=1/2: coefficient is n + 2*1/2 - 1 = n, so: -(n+2lambda-1) = -n.)

Then Q̃_n(x) = R_n(x + i/4) evaluated symbolically.

Wait: the existing unbalanced_zeros_single.py uses evalP_std and poly_coeffs_shifted.
Let me use their convention:
  evalP_std(n, y): standard MP lambda=1/2 evaluated at y
    (n+1) P_{n+1} = 2y P_n - (n+1/2) P_{n-1}  ???

Actually checking unbalanced_zeros_single.py:
  evalP_std: p0=1, p1=2y, p0,p1 = p1, (2y*p1 - j*p0)/(j+1)
  so: (j+1)*P_{j+1} = 2y*P_j - j*P_{j-1}
  => (n+1) P_{n+1} = 2x P_n - n P_{n-1}

This is lambda=1/2: n+2lambda-1 = n+2*1/2-1 = n. Correct.

And poly_coeffs_shifted builds these as polynomials in x with the shift x -> x + i/4:
  result[0] = [1]
  result[1] = [i/2, 2]    (i.e., 1*(i/4*2) + 0 for const, +2x... wait)

Looking at poly_coeffs_shifted: result[1] = [1j/2, mpc(2)].
So Q̃_1(x) = 2x + i/2? Let's verify: R_1(y) = 2y, Q̃_1(x) = R_1(x+i/4) = 2(x+i/4) = 2x + i/2. ✓

For the unbalanced matching:
  J_{N,k} = int Q̃_N(x) Q̃_k(x) W(x) dx
where W(x) = Gamma(1/4+ix)^2 * Gamma(3/4-ix).

The approximant is p_N(x) = sum_{k=0}^N J_{N,k}/nrm_k * Q̃_k(x)
where nrm_k = int |Q̃_k(x)|^2 * W_std dx for some standard weight...

Actually looking more carefully at unbalanced_dual_matching.pkl:
The pkl stores J_{m,k} directly as the raw integrals (not normalized).
The approximant in unbalanced_zeros_single.py uses J_ext[(m_row, k)] directly
without any normalization, and builds:
  f_coeffs[j] += J_{m,k} * poly_k[j]

So the polynomial is sum_k J_{m,k} Q̃_k(x) — the normalization is absorbed into J.
This means p_N(x) = sum_k J_{N,k} Q̃_k(x), using J_{N,k} as the raw matching integrals.

For large N, we need to compute J_{N,k} for k=0..N via quadrature.
The integral: J_{N,k} = int_{-inf}^{inf} Q̃_N(x) Q̃_k(x) * Gamma(1/4+ix)^2 * Gamma(3/4-ix) dx

This is a complex integrand.

Target: N=15 (fall back to N=10 if too slow).
"""
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pickle
import os

mp.mp.dps = 60

NODES = [-60, -40, -25, -15, -8, -4, -2, 0, 2, 4, 8, 15, 25, 40, 60]
N_TARGET = 15
CACHE_FILE = f"unbalanced_rowN_{N_TARGET}_cache.pkl"


def K_unbalanced(x):
    """Kernel: Gamma(1/4+ix)^2 * Gamma(3/4-ix)."""
    return mp.gamma(mp.mpf(1) / 4 + 1j * x) ** 2 * mp.gamma(mp.mpf(3) / 4 - 1j * x)


def poly_coeffs_shifted(N):
    """Polynomial coefficients of Q̃_0..Q̃_N in x (shifted by i/4)."""
    result = [[mp.mpc(1)]]
    if N == 0:
        return result
    result.append([1j / 2, mp.mpc(2)])
    p_prev, p_curr = result[0], result[1]
    for j in range(1, N):
        new_len = j + 2
        p_new = [mp.mpc(0)] * new_len
        # Q̃_{j+1} = ((j+1) terms from (j+1)*Q = 2y*Q_j - j*Q_{j-1})
        # but y = x + i/4, so 2y*Q̃_j = 2(x+i/4)*Q̃_j
        for idx, c in enumerate(p_curr):
            p_new[idx + 1] += 2 * c       # 2x term
            p_new[idx] += (1j / 2) * c    # 2*(i/4)*Q_j term
        for idx, c in enumerate(p_prev):
            p_new[idx] -= j * c
        p_new = [c / (j + 1) for c in p_new]
        result.append(p_new)
        p_prev = p_curr
        p_curr = p_new
    return result


def evalQ_shifted(n, x, pcs=None):
    """Evaluate Q̃_n(x) via recurrence."""
    # (n+1) Q̃_{n+1}(x) = 2(x+i/4) Q̃_n(x) - n Q̃_{n-1}(x)
    p0 = mp.mpc(1)
    if n == 0:
        return p0
    p1 = 2 * (x + 1j / 4)
    if n == 1:
        return p1
    for j in range(1, n):
        p0, p1 = p1, (2 * (x + 1j / 4) * p1 - j * p0) / (j + 1)
    return p1


print(f"Computing row-{N_TARGET} unbalanced coefficients J_{{{N_TARGET},k}} for k=0..{N_TARGET}...")

if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, 'rb') as f:
        J_row = pickle.load(f)
    print(f"  Loaded {len(J_row)} entries from cache")
else:
    # Check if existing pkl has some entries
    with open("unbalanced_dual_matching.pkl", "rb") as f:
        J_cached = pickle.load(f)

    J_row = {}
    for k in range(N_TARGET + 1):
        key = (N_TARGET, k)
        alt = (k, N_TARGET)
        if key in J_cached:
            J_row[k] = J_cached[key]
            print(f"  k={k}: from cache = {mp.nstr(J_row[k], 8)}")
        elif alt in J_cached:
            J_row[k] = J_cached[alt]
            print(f"  k={k}: from cache (symm) = {mp.nstr(J_row[k], 8)}")
        else:
            val = mp.quad(
                lambda x, k=k: evalQ_shifted(N_TARGET, x) * evalQ_shifted(k, x) * K_unbalanced(x),
                NODES)
            J_row[k] = val
            print(f"  k={k}: computed = {mp.nstr(val, 8)}")
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(J_row, f)
    print(f"  Saved {CACHE_FILE}")

# Coefficient polynomial: g_N(x) = sum_k J_{N,k} x^k using the raw
# matching integrals as monomial coefficients.
print(f"\nComputing zeros of g_{N_TARGET}(x) = sum_k J_{{N,k}} x^k ...")
f_coeffs = [J_row.get(k, mp.mpc(0)) for k in range(N_TARGET + 1)]
while len(f_coeffs) > 1 and abs(f_coeffs[-1]) < 1e-40:
    f_coeffs.pop()
print(f"  Degree: {len(f_coeffs) - 1}")

mp.mp.dps = 80
zs = mp.polyroots(f_coeffs[::-1], maxsteps=300, extraprec=300)
print(f"  {len(zs)} zeros found")

xs = [float(z.real) for z in zs]
ys = [float(z.imag) for z in zs]

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(xs, ys, s=18, color="royalblue", zorder=3,
           label=fr"zeros of $g_{{{N_TARGET}}}(x)=\sum_k J_{{{N_TARGET},k}}x^k$")
ax.axhline(0, color="gray", lw=0.5)
ax.axvline(0, color="gray", lw=0.5)
ax.set_xlabel("Re(x)")
ax.set_ylabel("Im(x)")
# legend removed
plt.tight_layout()
plt.savefig("unbalanced_rowN_zeros.png", dpi=150)
print("Saved unbalanced_rowN_zeros.png")
