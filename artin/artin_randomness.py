"""
Randomness tests on the S(3) membership sequence.

Treat the boolean sequence (is p_n in S(3)?) as a bitstream and run
standard PRNG tests. We expect failures — the interesting question
is *which* tests fail and what that says about the number theory.
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Load cached data
data = np.load('S3_50000.npz')
bits = data['membership'].astype(int)
n = len(bits)
p_hat = bits.mean()
print(f"Bitstream: {n} bits, density {p_hat:.4f}")
print(f"Expected for fair coin: 0.5000")
print(f"Expected for Artin: 0.3740\n")

results = []

def report(name, p_value, detail=""):
    verdict = "PASS" if p_value > 0.01 else "FAIL"
    results.append((name, p_value, verdict))
    print(f"{'PASS' if p_value > 0.01 else '** FAIL **':>10}  {name:40s}  p = {p_value:.6f}")
    if detail:
        print(f"{'':>10}  {detail}")


# ============================================================
# 1. Frequency (monobit) test
#    Is the proportion of 1s close to 0.5?
# ============================================================
S = 2 * bits.sum() - n  # sum of (+1/-1)
stat = abs(S) / np.sqrt(n)
p_val = 2 * (1 - stats.norm.cdf(stat))
report("1. Frequency (monobit) vs fair coin", p_val,
       f"proportion = {p_hat:.4f}, deviation from 0.5")

# Same test but against Artin's constant
artin = 0.3739558136
S_artin = bits.sum() - n * artin
stat_artin = abs(S_artin) / np.sqrt(n * artin * (1 - artin))
p_artin = 2 * (1 - stats.norm.cdf(stat_artin))
report("1b. Frequency vs Artin's constant", p_artin,
       f"proportion = {p_hat:.4f}, expected {artin:.4f}")

# ============================================================
# 2. Runs test
#    Are there too many or too few runs (consecutive 0s or 1s)?
# ============================================================
runs = 1 + np.sum(bits[1:] != bits[:-1])
expected_runs = 1 + 2 * n * p_hat * (1 - p_hat)  # under independence with observed p
var_runs = 2 * n * p_hat * (1 - p_hat) * (2 * n * p_hat * (1 - p_hat) - 1) / (n - 1)
z_runs = (runs - expected_runs) / np.sqrt(var_runs)
p_runs = 2 * (1 - stats.norm.cdf(abs(z_runs)))
report("2. Runs test (independence)", p_runs,
       f"observed {runs}, expected {expected_runs:.0f}")

# ============================================================
# 3. Serial test (pairs)
#    Count 00, 01, 10, 11 bigrams. Under independence these
#    should follow p(ij) = p(i)*p(j).
# ============================================================
pairs = bits[:-1] * 2 + bits[1:]  # encode as 0,1,2,3
observed = np.bincount(pairs, minlength=4)
p0, p1 = 1 - p_hat, p_hat
expected = np.array([p0*p0, p0*p1, p1*p0, p1*p1]) * (n - 1)
chi2 = np.sum((observed - expected)**2 / expected)
p_serial = 1 - stats.chi2.cdf(chi2, df=1)  # df = 4 - 2 - 1 = 1
report("3. Serial test (bigrams)", p_serial,
       f"χ² = {chi2:.2f}, observed {observed}, expected {expected.astype(int)}")

# ============================================================
# 4. Autocorrelation test
#    Test correlation at lags 1..50
# ============================================================
centered = bits - p_hat
max_lag = 50
autocorr = np.correlate(centered, centered, mode='full')
autocorr = autocorr[n-1:] / autocorr[n-1]  # normalise by lag-0
lag_pvals = []
for lag in range(1, max_lag + 1):
    # under independence, r(lag) ~ N(0, 1/n)
    z = autocorr[lag] * np.sqrt(n)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    lag_pvals.append((lag, autocorr[lag], p))

sig_lags = [(lag, r, p) for lag, r, p in lag_pvals if p < 0.01]
if sig_lags:
    report("4. Autocorrelation (lags 1-50)", min(p for _, _, p in sig_lags),
           f"{len(sig_lags)} significant lags: {[(l,f'{r:.4f}') for l,r,_ in sig_lags[:10]]}")
else:
    report("4. Autocorrelation (lags 1-50)", max(p for _, _, p in lag_pvals),
           "no significant lags")

# ============================================================
# 5. Block frequency test
#    Divide into blocks of M bits. Is the proportion in each
#    block consistent with a common p?
# ============================================================
M = 100
num_blocks = n // M
block_props = np.array([bits[i*M:(i+1)*M].mean() for i in range(num_blocks)])
chi2_block = (M / (p_hat * (1 - p_hat))) * np.sum((block_props - p_hat)**2)
p_block = 1 - stats.chi2.cdf(chi2_block, df=num_blocks - 1)
report("5. Block frequency (M=100)", p_block,
       f"χ² = {chi2_block:.1f}, df = {num_blocks-1}")

# ============================================================
# 6. Longest run of ones test
# ============================================================
M_lr = 1000
num_blocks_lr = n // M_lr
longest_runs = []
for i in range(num_blocks_lr):
    block = bits[i*M_lr:(i+1)*M_lr]
    max_run = 0
    current = 0
    for b in block:
        if b == 1:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 0
    longest_runs.append(max_run)
longest_runs = np.array(longest_runs)
mean_lr = np.mean(longest_runs)
# expected longest run of 1s in M Bernoulli(p) trials ≈ log(M*p) / log(1/p)
expected_lr = np.log(M_lr * p_hat) / np.log(1 / p_hat)
report("6. Longest run of 1s (blocks of 1000)",
       stats.ttest_1samp(longest_runs, expected_lr).pvalue,
       f"mean longest = {mean_lr:.1f}, expected ≈ {expected_lr:.1f}")

# ============================================================
# 7. Cumulative sums (CUSUM) test
# ============================================================
walk = np.cumsum(2 * bits - 1)  # ±1 random walk (biased)
# adjust for bias
walk_adj = np.cumsum(bits - p_hat) / np.sqrt(n * p_hat * (1 - p_hat))
ks_stat, ks_p = stats.kstest(walk_adj / np.sqrt(np.arange(1, n+1) / n), 'norm')
report("7. CUSUM (Kolmogorov-Smirnov)", ks_p,
       f"KS stat = {ks_stat:.4f}")

# ============================================================
# 8. Spectral (DFT) test
#    Look for periodic structure in the bitstream
# ============================================================
centered_bits = bits - p_hat
fft = np.fft.rfft(centered_bits)
magnitudes = np.abs(fft[1:])  # drop DC component
threshold = np.sqrt(n * p_hat * (1 - p_hat) * np.log(1/0.05))  # 95% threshold
num_peaks = np.sum(magnitudes > threshold)
expected_peaks = 0.05 * len(magnitudes)
# test if number of peaks exceeding threshold is consistent with random
z_spec = (num_peaks - expected_peaks) / np.sqrt(expected_peaks * 0.95)
p_spec = 2 * (1 - stats.norm.cdf(abs(z_spec)))
report("8. Spectral (DFT) peaks", p_spec,
       f"{num_peaks} peaks above 95% threshold (expected {expected_peaks:.0f})")

# find the dominant frequencies
top_freqs = np.argsort(magnitudes)[-5:][::-1]
print(f"{'':>10}  Top 5 frequencies: {top_freqs + 1}")
print(f"{'':>10}  Corresponding periods: {[f'{n/(f+1):.1f}' for f in top_freqs]} primes")

# ============================================================
# Summary
# ============================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
failures = [r for r in results if r[2] == "FAIL"]
passes = [r for r in results if r[2] == "PASS"]
print(f"Passed: {len(passes)}/{len(results)}")
print(f"Failed: {len(failures)}/{len(results)}")
for name, pval, _ in failures:
    print(f"  FAIL: {name}  (p = {pval:.2e})")

# ============================================================
# Plot: autocorrelation and spectral
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))

lags = [l for l, _, _ in lag_pvals]
acorrs = [r for _, r, _ in lag_pvals]
acpvals = [p for _, _, p in lag_pvals]
colors = ['firebrick' if p < 0.01 else 'steelblue' for p in acpvals]
ax1.bar(lags, acorrs, color=colors, alpha=0.8)
ax1.axhline(0, color='black', lw=0.5)
ax1.axhline(2/np.sqrt(n), color='grey', ls='--', lw=0.8, label='±2σ')
ax1.axhline(-2/np.sqrt(n), color='grey', ls='--', lw=0.8)
ax1.set_xlabel('Lag')
ax1.set_ylabel('Autocorrelation')
ax1.set_title('S(3) membership autocorrelation (red = significant at 1%)')
ax1.legend()

freqs = np.arange(1, len(magnitudes) + 1)
ax2.semilogy(freqs[:500], magnitudes[:500], color='steelblue', alpha=0.5, lw=0.5)
ax2.axhline(threshold, color='firebrick', ls='--', label=f'95% threshold')
ax2.set_xlabel('Frequency index')
ax2.set_ylabel('|FFT| magnitude')
ax2.set_title('Spectral analysis of S(3) bitstream (first 500 frequencies)')
ax2.legend()

plt.tight_layout()
plt.savefig('artin_randomness.png', dpi=150)
print("\nSaved artin_randomness.png")
