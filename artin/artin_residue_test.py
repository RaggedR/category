"""
A statistical test that S(3) passes and a PRNG fails.

The test: compute P(in S(r) | p ≡ a mod m) for small moduli m.
Number theory predicts these split by residue class (Chebotarev).
A PRNG has no access to the underlying prime, so the conditional
probabilities are flat.
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Load S(3) data
data = np.load('S3_50000.npz')
primes = data['primes']
bits = data['membership'].astype(int)
n = len(bits)
p_overall = bits.mean()

print(f"S(3): {n} primes, density {p_overall:.4f}\n")

# ============================================================
# Part 1: Residue-class splitting in S(3)
# ============================================================
moduli = [3, 4, 8, 12, 24]

print("=" * 70)
print("RESIDUE-CLASS SPLITTING IN S(3)")
print("=" * 70)

splits = {}
for m in moduli:
    residues = primes % m
    classes = sorted(set(residues))
    print(f"\nmod {m}:")
    class_data = []
    for a in classes:
        mask = residues == a
        count = mask.sum()
        if count < 10:
            continue
        hit_rate = bits[mask].mean()
        # z-test against overall density
        se = np.sqrt(p_overall * (1 - p_overall) / count)
        z = (hit_rate - p_overall) / se
        p_val = 2 * (1 - stats.norm.cdf(abs(z)))
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        print(f"  p ≡ {a:2d} mod {m}: n={count:6d}  P(S(3))={hit_rate:.4f}  z={z:+6.2f}  {sig}")
        class_data.append((a, count, hit_rate))
    splits[m] = class_data

    # chi-squared test for uniformity across classes
    observed_hits = [bits[residues == a].sum() for a in classes if (residues == a).sum() >= 10]
    expected_hits = [p_overall * (residues == a).sum() for a in classes if (residues == a).sum() >= 10]
    chi2 = sum((o - e)**2 / e for o, e in zip(observed_hits, expected_hits))
    df = len(observed_hits) - 1
    p_chi = 1 - stats.chi2.cdf(chi2, df)
    print(f"  χ² = {chi2:.1f}, df = {df}, p = {p_chi:.2e}")

# ============================================================
# Part 2: Same test on a bias-matched PRNG
# ============================================================
print("\n" + "=" * 70)
print("SAME TEST ON BIAS-MATCHED PRNG")
print("=" * 70)

rng = np.random.default_rng(42)
# Generate a PRNG bitstream with same density as S(3)
# Also add lag-1 autocorrelation to match S(3)'s r1 ≈ -0.06
fake_bits = np.zeros(n, dtype=int)
fake_bits[0] = int(rng.random() < p_overall)
r1_target = -0.06
p_after_1 = p_overall + r1_target * p_overall * (1 - p_overall) / p_overall  # approximate
p_after_0 = p_overall - r1_target * p_overall * (1 - p_overall) / (1 - p_overall)
for i in range(1, n):
    if fake_bits[i-1] == 1:
        fake_bits[i] = int(rng.random() < p_after_1)
    else:
        fake_bits[i] = int(rng.random() < p_after_0)

print(f"PRNG: density = {fake_bits.mean():.4f}, lag-1 r = {np.corrcoef(fake_bits[:-1], fake_bits[1:])[0,1]:.4f}")
print(f"(matched to S(3): density {p_overall:.4f}, lag-1 r = -0.062)\n")

for m in moduli:
    residues = primes % m
    classes = sorted(set(residues))
    print(f"mod {m}:")
    for a in classes:
        mask = residues == a
        count = mask.sum()
        if count < 10:
            continue
        hit_rate = fake_bits[mask].mean()
        se = np.sqrt(p_overall * (1 - p_overall) / count)
        z = (hit_rate - fake_bits.mean()) / se
        p_val = 2 * (1 - stats.norm.cdf(abs(z)))
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        print(f"  p ≡ {a:2d} mod {m}: n={count:6d}  P(hit)={hit_rate:.4f}  z={z:+6.2f}  {sig}")

    observed_hits = [fake_bits[residues == a].sum() for a in classes if (residues == a).sum() >= 10]
    expected_hits = [fake_bits.mean() * (residues == a).sum() for a in classes if (residues == a).sum() >= 10]
    chi2 = sum((o - e)**2 / e for o, e in zip(observed_hits, expected_hits))
    df = len(observed_hits) - 1
    p_chi = 1 - stats.chi2.cdf(chi2, df)
    print(f"  χ² = {chi2:.1f}, df = {df}, p = {p_chi:.2e}")

# ============================================================
# Part 3: Visualise the split
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# mod 12 is the most informative for r=3
m = 12
residues = primes % m
classes = [a for a in sorted(set(residues)) if (residues == a).sum() >= 10]

# S(3)
rates_real = [bits[residues == a].mean() for a in classes]
rates_fake = [fake_bits[residues == a].mean() for a in classes]

x = np.arange(len(classes))
w = 0.35
axes[0].bar(x - w/2, rates_real, w, color='steelblue', alpha=0.8, label='S(3)')
axes[0].bar(x + w/2, rates_fake, w, color='firebrick', alpha=0.8, label='PRNG')
axes[0].axhline(p_overall, color='black', ls='--', lw=0.8, label=f'overall = {p_overall:.3f}')
axes[0].set_xticks(x)
axes[0].set_xticklabels([f'{a}' for a in classes])
axes[0].set_xlabel(f'p mod {m}')
axes[0].set_ylabel('P(primitive root)')
axes[0].set_title(f'Residue-class splitting mod {m}')
axes[0].legend(fontsize=9)

# mod 24
m = 24
residues = primes % m
classes = [a for a in sorted(set(residues)) if (residues == a).sum() >= 10]
rates_real = [bits[residues == a].mean() for a in classes]
rates_fake = [fake_bits[residues == a].mean() for a in classes]

x = np.arange(len(classes))
axes[1].bar(x - w/2, rates_real, w, color='steelblue', alpha=0.8, label='S(3)')
axes[1].bar(x + w/2, rates_fake, w, color='firebrick', alpha=0.8, label='PRNG')
axes[1].axhline(p_overall, color='black', ls='--', lw=0.8, label=f'overall = {p_overall:.3f}')
axes[1].set_xticks(x)
axes[1].set_xticklabels([f'{a}' for a in classes], fontsize=8)
axes[1].set_xlabel(f'p mod {m}')
axes[1].set_ylabel('P(primitive root)')
axes[1].set_title(f'Residue-class splitting mod {m}')
axes[1].legend(fontsize=9)

plt.suptitle('Test: S(3) "knows" the prime\'s residue class. A PRNG does not.', fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig('artin_residue_test.png', dpi=150, bbox_inches='tight')
print("\nSaved artin_residue_test.png")
