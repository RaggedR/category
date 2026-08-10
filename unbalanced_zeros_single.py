"""Single-panel: zeros of unbalanced dual F_4^(N) for largest feasible N.
Uses GF to extend the J_mk table beyond the cached MMAX=11.
Keeps the pole-line Im x = -3/4 reference."""
import mpmath as mp
import pickle
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 60

NODES = [-60, -40, -25, -15, -8, -4, -2, 0, 2, 4, 8, 15, 25, 40, 60]

# Load existing J table
with open("unbalanced_dual_matching.pkl", "rb") as f:
    J_cached = pickle.load(f)
print(f"Loaded cache: max index = {max(max(k) for k in J_cached.keys())}")


def evalP_std(n, y):
    p0, p1 = mp.mpf(1), 2*y
    if n == 0: return p0
    if n == 1: return p1
    for j in range(1, n):
        p0, p1 = p1, (2*y*p1 - j*p0) / (j+1)
    return p1


def GF_master(w, u):
    eta = (1 + 1j*w) * (1 + 1j*u) / ((1 - 1j*w) * (1 - 1j*u))
    return 2 * mp.pi * mp.exp(eta) * mp.e1(eta) / ((1 - 1j*w) * (1 - 1j*u))


def K3(x):
    return mp.gamma(mp.mpf(1)/4 + 1j*x)**2 * mp.gamma(mp.mpf(3)/4 - 1j*x)


def compute_J_quad(m, k):
    """Compute J_mk^{std} via quadrature."""
    val = mp.quad(
        lambda x: evalP_std(m, x + 1j/4) * evalP_std(k, x + 1j/4) * K3(x),
        NODES)
    return val


def get_J(m, k):
    key = (m, k)
    alt = (k, m)
    if key in J_cached:
        return J_cached[key]
    if alt in J_cached:
        return J_cached[alt]
    return None


# Target: m=4, N=30. We need J_{4,k} for k=0..30.
# The cache has k up to 11. We need to extend to k=30.
# Use quadrature for missing entries.
m_row = 4
N_target = 30

CACHE_EXT = "unbalanced_J4_ext.pkl"
if os.path.exists(CACHE_EXT):
    with open(CACHE_EXT, "rb") as f:
        J_ext = pickle.load(f)
    print(f"Loaded extended cache: {len(J_ext)} entries")
else:
    print(f"Computing J_{{4,k}} for k=0..{N_target} via quadrature...")
    J_ext = {}
    for k in range(N_target + 1):
        existing = get_J(m_row, k)
        if existing is not None:
            J_ext[(m_row, k)] = existing
            print(f"  k={k}: from cache = {mp.nstr(existing, 8)}")
        else:
            val = compute_J_quad(m_row, k)
            J_ext[(m_row, k)] = val
            print(f"  k={k}: computed = {mp.nstr(val, 8)}")
    with open(CACHE_EXT, "wb") as f:
        pickle.dump(J_ext, f)
    print(f"Saved {CACHE_EXT}")


def poly_coeffs_shifted(N):
    result = [[mp.mpc(1)]]
    if N == 0: return result
    result.append([1j/2, mp.mpc(2)])
    p_prev, p_curr = result[0], result[1]
    for j in range(1, N):
        new_len = j + 2
        p_new = [mp.mpc(0)] * new_len
        for idx, c in enumerate(p_curr):
            p_new[idx + 1] += 2 * c
            p_new[idx] += 1j/2 * c
        for idx, c in enumerate(p_prev):
            p_new[idx] -= j * c
        p_new = [c / (j + 1) for c in p_new]
        result.append(p_new)
        p_prev = p_curr
        p_curr = p_new
    return result


print(f"\nBuilding polynomial coefficients for N={N_target}...")
pcs = poly_coeffs_shifted(N_target)

# Build truncated function F_4^(N)
f_coeffs = [mp.mpc(0)] * (N_target + 1)
for k in range(N_target + 1):
    Jmk = J_ext.get((m_row, k), mp.mpc(0))
    pk = pcs[k]
    for idx, c in enumerate(pk):
        if idx < len(f_coeffs):
            f_coeffs[idx] += Jmk * c

while len(f_coeffs) > 1 and abs(f_coeffs[-1]) < 1e-30:
    f_coeffs.pop()

print(f"Finding roots of degree-{len(f_coeffs)-1} polynomial...")
mp.mp.dps = 80
roots = mp.polyroots(f_coeffs[::-1], maxsteps=300, extraprec=300)
mp.mp.dps = 60
good = [(float(z.real), float(z.imag)) for z in roots if abs(z.imag) < 5]
print(f"  {len(good)} zeros with |Im| < 5")

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter([r[0] for r in good], [r[1] for r in good],
           s=14, c='steelblue', alpha=0.8, zorder=5,
           label=r"zeros of $F_4^{(30)}$")
ax.axhline(y=0, color='gray', linewidth=0.5)
ax.axhline(y=-0.75, color='darkred', linewidth=1.2, linestyle='--',
           label=r'Im$\,x = -\frac{3}{4}$', alpha=0.85)
ax.axhline(y=-0.25, color='crimson', linewidth=0.8, linestyle=':',
           label=r'Im$\,x = -\frac{1}{4}$', alpha=0.7)
ax.set_title(
    r"Zeros of $F_4^{(30)}$: $W = \Gamma(\frac{1}{4}+ix)\Gamma(\frac{3}{4}-ix)$",
    fontsize=12)
ax.set_xlabel("Re(x)")
ax.set_ylabel("Im(x)")
ax.set_ylim(-2.0, 1.5)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("unbalanced_zeros_single.png", dpi=150)
print("Saved unbalanced_zeros_single.png")
