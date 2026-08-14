"""
Artin's conjecture: S(r) = { q prime : r generates (ℤ/qℤ)* }.

Bin consecutive primes into blocks of 100 and count how many
in each block belong to S(r).  Plot the histogram of counts
and the Shannon entropy of fluctuations around the Gaussian fit.
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt


def sieve(n):
    """Sieve of Eratosthenes up to n."""
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = False
    return np.nonzero(is_prime)[0]


def prime_factors(n):
    """Distinct prime factors of n by trial division."""
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def is_primitive_root(r, q):
    """Check if r is a primitive root mod q."""
    if q <= 2 or r % q == 0:
        return False
    order = q - 1
    for p in prime_factors(order):
        if pow(r, order // p, q) == 1:
            return False
    return True


r = 3
num_primes = 50000
block_size = 100
cache_file = f'S{r}_{num_primes}.npz'

import os
if os.path.exists(cache_file):
    data = np.load(cache_file)
    primes = data['primes']
    membership = data['membership']
    print(f"Loaded {len(primes)} primes from {cache_file}")
else:
    # sieve enough primes (p_50000 ≈ 611953)
    primes = sieve(700000)[:num_primes]
    membership = np.array([is_primitive_root(r, int(q)) for q in primes])
    np.savez(cache_file, primes=primes, membership=membership)
    print(f"Computed and saved to {cache_file}")

print(f"Got {len(primes)} primes up to {primes[-1]}")
print(f"|S({r})| = {membership.sum()} out of {len(primes)}  ({membership.mean():.4f})")

# count S(r) members per block
counts = []
for i in range(0, len(primes), block_size):
    block_mask = membership[i:i + block_size]
    if len(block_mask) < block_size:
        break
    counts.append(block_mask.sum())

counts = np.array(counts)
artin_const = 0.3739558136
mean_c = counts.mean()
std_c = counts.std()

print(f"Blocks: {len(counts)}")
print(f"Mean count per block: {mean_c:.2f}  (Artin predicts ~{100 * artin_const:.1f})")
print(f"Std: {std_c:.2f}")
print(f"Range: [{counts.min()}, {counts.max()}]")

kurt = ((counts - mean_c) ** 4).mean() / std_c ** 4 - 3
print(f"Excess kurtosis: {kurt:.3f}  (Gaussian = 0, positive = tall & skinny)")

# --- Histogram ---
fig, ax = plt.subplots(figsize=(8, 5))
bins = np.arange(counts.min() - 0.5, counts.max() + 1.5, 1)
ax.hist(counts, bins=bins, density=True, alpha=0.7, color='steelblue',
        edgecolor='white', label=f'S({r}) block counts')

x = np.linspace(counts.min() - 2, counts.max() + 2, 200)
gauss = np.exp(-0.5 * ((x - mean_c) / std_c) ** 2) / (std_c * np.sqrt(2 * np.pi))
ax.plot(x, gauss, 'r-', lw=2, label=f'Gaussian (μ={mean_c:.1f}, σ={std_c:.1f})')

ax.set_xlabel('Primes in S(3) per block of 100 consecutive primes')
ax.set_ylabel('Density')
ax.set_title(f"Artin's conjecture: distribution of S({r}) density across blocks")
ax.legend()
ax.text(0.98, 0.95, f'κ = {kurt:.3f}', transform=ax.transAxes,
        ha='right', va='top', fontsize=11)
plt.tight_layout()
plt.savefig('artin_histogram.png', dpi=150)
plt.show()
print("Saved artin_histogram.png")

# --- Shannon entropy of fluctuations around Gaussian ---
print("\n=== Fluctuation entropy ===")

vals, emp_counts = np.unique(counts, return_counts=True)
p_emp = emp_counts / emp_counts.sum()

p_gauss = np.exp(-0.5 * ((vals - mean_c) / std_c) ** 2) / (std_c * np.sqrt(2 * np.pi))
p_gauss = p_gauss / p_gauss.sum()

H_emp = -np.sum(p_emp * np.log2(p_emp))
H_gauss = -np.sum(p_gauss * np.log2(p_gauss))

mask = p_emp > 0
D_KL = np.sum(p_emp[mask] * np.log2(p_emp[mask] / p_gauss[mask]))

print(f"H(empirical)  = {H_emp:.4f} bits")
print(f"H(Gaussian)   = {H_gauss:.4f} bits")
print(f"D_KL(emp||G)  = {D_KL:.6f} bits")
print(f"  → empirical is {'lower' if H_emp < H_gauss else 'higher'} entropy than Gaussian")
print(f"  → {'more concentrated (tall & skinny)' if H_emp < H_gauss else 'more spread out'}")

# fluctuation plot
fig2, ax2 = plt.subplots(figsize=(8, 4))
delta = p_emp - p_gauss
ax2.bar(vals, delta, color=['firebrick' if d > 0 else 'steelblue' for d in delta],
        alpha=0.8, edgecolor='white')
ax2.axhline(0, color='black', lw=0.5)
ax2.set_xlabel('Count k')
ax2.set_ylabel('P_emp(k) − P_Gauss(k)')
ax2.set_title(f'Fluctuations around Gaussian  (D_KL = {D_KL:.5f} bits)')
plt.tight_layout()
plt.savefig('artin_fluctuations.png', dpi=150)
plt.show()
print("Saved artin_fluctuations.png")
