# Artin's Conjecture — Block Density Fluctuations

## Setup

S(r) = { q prime : r generates (ℤ/qℤ)* }. Artin's conjecture says |S(r)| is infinite
with natural density equal to Artin's constant C ≈ 0.3739558136 for any non-square r ≠ −1.

We computed S(3) over the first 50,000 primes (up to 611,953), binned them into 500 blocks
of 100 consecutive primes, and counted how many primes in each block belong to S(3).

## Results

| Quantity | Value |
|----------|-------|
| \|S(3)\| in first 50K primes | 18,760 (density 0.3752) |
| Artin's constant | 0.3740 |
| Mean block count | 37.52 |
| Std deviation | 3.95 |
| Range | [27, 48] |
| Excess kurtosis | −0.512 |

## Fluctuation entropy

| Quantity | Value |
|----------|-------|
| H(empirical) | 3.9909 bits |
| H(Gaussian) | 3.9887 bits |
| D_KL(empirical \|\| Gaussian) | 0.0308 bits |

The empirical distribution has *higher* entropy than the best-fit Gaussian — it is more
spread out, not more concentrated. The excess kurtosis is negative (platykurtic): the
distribution is flat-topped relative to a Gaussian, with excess mass at k ≈ 32 and k ≈ 41
and deficits at the peak (k ≈ 37) and shoulders.

## Randomness tests

We treated the S(3) membership sequence as a bitstream and ran a battery of standard
PRNG tests. 5 of 9 tests failed.

### Failures

| Test | p-value | What it means |
|------|---------|---------------|
| Frequency (monobit vs 0.5) | ≈ 0 | Trivial: density is 0.375, not 0.5 |
| Runs | ≈ 0 | Too many runs (24,893 vs 23,443 expected). Sequence switches 0↔1 more than independence predicts. |
| Serial (bigrams) | ≈ 0 | Too many 01/10 pairs, too few 00/11. χ² = 191.2. Anti-clustering. |
| Autocorrelation | ≈ 0 | Lag 1: r = −0.062. Lag 3: r = −0.025. Lag 5: r = −0.015. Negative at odd lags. |
| CUSUM | ≈ 0 | Cumulative sum drifts (consequence of bias + correlation) |

### Passes

| Test | p-value | What it means |
|------|---------|---------------|
| Frequency vs Artin's constant | 0.565 | Density is spot-on Artin's constant |
| Block frequency (M=100) | 1.000 | Density is stable across blocks |
| Longest run of 1s | 0.018 | Marginal but passes (mean 5.7, expected 6.0) |
| Spectral (DFT) | 0.296 | No periodic structure |

### Interpretation

The sequence has correct long-range statistics (right density, stable across blocks, no
periodicity) but wrong short-range correlations. The dominant signal is **negative lag-1
autocorrelation**: consecutive primes repel each other in S(3).

If p is in S(3), then 3 has maximal order mod p, meaning p−1 has no small factor that
kills 3's order. The next prime p' = p + 2k often shares factors of p'−1 with p−1 (both
are even, often both divisible by 3, etc.). If p−1 has the "right" factorisation for 3 to
be a primitive root, p'−1 is slightly more likely to have a "wrong" factorisation — because
the small primes dividing p'−1 are correlated with those dividing p−1. This is Chebotarev
repulsion: the Frobenius elements in Gal(ℚ(ζ_q, q^{1/3})/ℚ) are equidistributed but not
independent.

The negative kurtosis (−0.512) in the block-count distribution is a macroscopic consequence
of this microscopic repulsion: blocks are more uniform than a binomial model predicts.

## Repulsion across r

The repulsion effect varies dramatically with r. Computed over 50K primes for non-square
r from 2 to 20:

| r | P(in S(r)) | P after hit | P after miss | Drop | Lag-1 autocorr |
|---|-----------|-------------|--------------|------|----------------|
| 2 | 0.3740 | 0.3290 | 0.4010 | +12.0% | −0.072 |
| 3 | 0.3752 | 0.3366 | 0.3984 | +10.3% | −0.062 |
| 5 | 0.3940 | 0.3376 | 0.4306 | +14.3% | −0.093 |
| 6 | 0.3736 | 0.3505 | 0.3874 | +6.2% | −0.037 |
| 7 | 0.3737 | 0.3690 | 0.3765 | +1.3% | −0.008 |
| 10 | 0.3756 | 0.3623 | 0.3837 | +3.6% | −0.021 |
| 11 | 0.3754 | 0.3590 | 0.3853 | +4.4% | −0.026 |
| 12 | 0.3755 | 0.3326 | 0.4013 | +11.4% | −0.069 |
| 13 | 0.3772 | 0.3459 | 0.3962 | +8.3% | −0.050 |
| 14 | 0.3735 | 0.3738 | 0.3734 | −0.1% | +0.000 |
| **15** | **0.3765** | **0.4195** | **0.3505** | **−11.4%** | **+0.069** |
| 17 | 0.3759 | 0.3503 | 0.3913 | +6.8% | −0.041 |
| 18 | 0.3730 | 0.3294 | 0.3989 | +11.7% | −0.070 |
| 19 | 0.3759 | 0.3733 | 0.3775 | +0.7% | −0.004 |
| 20 | 0.3937 | 0.3352 | 0.4316 | +14.8% | −0.096 |

Most values of r show repulsion (negative autocorrelation), but **r = 15 attracts** —
consecutive primes cluster in S(15), with an 11.4% boost after a hit. r = 14 is neutral.

The strongest repulsors (r = 2, 5, 12, 18, 20) and the lone attractor (r = 15) suggest the
effect depends on the arithmetic of r, not just its size. The pattern deserves further
investigation: what property of r controls the sign and magnitude of the lag-1 autocorrelation?

## Files

- `artin_histogram.py` — block density histogram + Shannon entropy analysis
- `artin_randomness.py` — PRNG test battery (monobit, runs, serial, autocorrelation, block freq, longest run, CUSUM, spectral)
- `S3_50000.npz` — cached data: `primes` (first 50K primes) and `membership` (boolean mask for S(3))
- `artin_histogram.png` — histogram with Gaussian overlay
- `artin_fluctuations.png` — residuals P_emp(k) − P_Gauss(k)
- `artin_randomness.png` — autocorrelation and spectral plots
