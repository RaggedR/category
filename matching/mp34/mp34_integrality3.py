"""Finding the correct integer normalization for lambda=3/4.

Key observation: The denominators appearing in T_mk = p(m,k)*A' + q(m,k)*B'
are products of (2*1+1)(2*2+1)... type factors (from Pochhammer symbols).

For lambda=3/4: beta(n) = n! Gamma(3/2)/Gamma(n+3/2) = n!/((3/2)_n).
The inverse: (3/2)_n = Gamma(n+3/2)/Gamma(3/2).

For lambda=1/4: the correct normalization was
  S_mk = 2^m m! * 2^k k! / h_0 * J_mk
which in terms of the GF coefficients meant multiplying by [coeff^{-1} of w^n in 1/sqrt(1-iw)^{-1}]
= 1 / (i^n * C(2n,n)/4^n) = (-i)^n * 4^n / C(2n,n).

For lambda=3/4, the analogous inverse of [coeff of u^k in 1/(1+iu)^{3/2}]:
[u^k] (1+iu)^{-3/2} = (-i)^k * (3/2)_k / k!
Inverse normalization: k! * (-i)^{-k} / (3/2)_k = (-1)^k * i^k * k! / (3/2)_k = i^k * k! / (3/2)_k.

So the right definition of the integer table should be:
  S_mk = i^{m+k} * m! * k! / ((3/2)_m * (3/2)_k) * J_mk / pref

where pref = pi^2/2.

Let's test this.
"""
import mpmath as mp
import pickle
from fractions import Fraction

mp.mp.dps = 50

with open("mp34_dual_matching.pkl", "rb") as f:
    nrm, C = pickle.load(f)

def c(m, k):
    return C[(min(m,k), max(m,k))]

def J(m, k):
    return c(m, k) * mp.sqrt(nrm[m] * nrm[k])

h0 = nrm[0]
J00 = J(0, 0)

Ap = mp.hyperu(mp.mpf(3)/2, 1, 1)   # U(3/2,1,1)
Bp = mp.hyperu(mp.mpf(5)/2, 2, 1)   # U(5/2,2,1)
pref = mp.pi**2 / 2

def Poch(n):
    """(3/2)_n = Gamma(n+3/2)/Gamma(3/2)."""
    return mp.gamma(n + mp.mpf(3)/2) / mp.gamma(mp.mpf(3)/2)

print("=== Testing S_mk = i^{m+k} * m! * k! / (Poch_m * Poch_k) * J_mk / pref ===")
for m in range(6):
    for k in range(m, 6):
        Jmk = J(m, k)
        S_raw = 1j**(m+k) * mp.fac(m) * mp.fac(k) / (Poch(m) * Poch(k)) * Jmk / pref
        S = mp.re(S_raw)  # should be real
        err = abs(mp.im(S_raw))
        r = mp.pslq([S, Ap, Bp], tol=1e-35, maxcoeff=10000000)
        if r is not None:
            cA = mp.mpf(-r[1]) / mp.mpf(r[0])
            cB = mp.mpf(-r[2]) / mp.mpf(r[0])
            fA = Fraction(float(cA)).limit_denominator(10000)
            fB = Fraction(float(cB)).limit_denominator(10000)
            print(f"S_{m},{k} = {fA}*A' + {fB}*B'")
        else:
            print(f"S_{m},{k} = {mp.nstr(S, 8)} (PSLQ failed)")

# The correct T was: T_mk = beta_m * beta_k * (-i)^{m+k} * J_mk / pref
# where beta(n) = n!/Poch(n). So:
# S_mk = i^{m+k} * m! * k! / Poch_m Poch_k * J_mk / pref
# = i^{m+k} * (beta_m * Poch_m / Poch_m) ... same as (-i)^{m+k} but with opposite sign for i!
# Actually i^{m+k} vs (-i)^{m+k}: for m+k even, same; for m+k odd, opposite sign.
# But T_mk was taking re(), so the sign matters for odd.
# Let's check: T_01 = A' - 2B' (negative), S_01 should be -(A'-2B') = 2B'-A'.

print("\n--- Checking S_01 ---")
J01 = J(0, 1)
S01_raw = 1j * mp.fac(0) * mp.fac(1) / (Poch(0) * Poch(1)) * J01 / pref
print(f"S_01 = {mp.nstr(S01_raw, 12)}")
print(f"2B' - A' = {mp.nstr(2*Bp - Ap, 12)}")

# Also try with the T normalization (minus sign for odd):
# S_mk = (-1)^{sign} * T_mk, where we want S_mk > 0 and integers.

# Maybe the right normalization should absorb the pi^2/2 differently.
# For lambda=1/4: h_0 = sqrt(2) pi^{3/2}, so 2pi*Gamma(1/2)^2 = 2pi*pi = 2pi^2,
# and J_00 = 2pi^{3/2} sqrt(e) K_0(1/2) = 2pi^2 / pi * sqrt(e) K_0(1/2) / sqrt(pi) = ...
# Actually J_00^{1/4} = 2pi^{3/2} sqrt(e) K_0(1/2) and h_0^{1/4} = sqrt(2) pi^{3/2}.
# c_00 = J_00/h_0 = sqrt(e)*K_0(1/2)*sqrt(2) = A.

# For lambda=3/4: J_00 = pi^2/2 * U(3/2,1,1), h_0 = J_00/c_00,
# c_00 = J_00/h_0 = 0.897... (from table).

# The issue: for lambda=1/4 the integer table used the NORMALIZATION
# T_mk = m! k! sqrt(C(2m,m) C(2k,k)) * c_mk.
# These were integers times A (or i*integers times A, B).
# The denominator structure: for l=1/4, h_n = h_0 C(2n,n)/4^n, and
# T_mk = m! k! * sqrt(h_m/h_0 * h_k/h_0) * 4^{(m+k)/2} * c_mk
# = m! k! * (C(2m,m) C(2k,k))^{1/2} * c_mk.  -- these were in Z*A+Z*B.

# For lambda=3/4, try: T_mk = m! k! * sqrt(Poch_m/1 * Poch_k/1) * c_mk (analogous)
# = m! k! * sqrt(Poch_m * Poch_k) * c_mk.
# Let me test this.
print("\n=== Testing T_mk = m! k! sqrt(Poch_m * Poch_k) * c_mk ===")
for m in range(6):
    for k in range(m, 6):
        cmk = c(m, k)
        Pm = Poch(m); Pk = Poch(k)
        S_raw = mp.fac(m) * mp.fac(k) * mp.sqrt(Pm * Pk) * cmk
        S = mp.re(S_raw) if (m+k)%2==0 else mp.im(S_raw)
        # S is real or purely imaginary based on parity
        # For even: c_mk real, so S real. For odd: c_mk imaginary, S imaginary.
        r = mp.pslq([S, Ap, Bp], tol=1e-35, maxcoeff=10000000)
        if r is not None:
            cA = mp.mpf(-r[1]) / mp.mpf(r[0])
            cB = mp.mpf(-r[2]) / mp.mpf(r[0])
            fA = Fraction(float(cA)).limit_denominator(10000)
            fB = Fraction(float(cB)).limit_denominator(10000)
            print(f"T_{m},{k} = {fA}*A' + {fB}*B'  (S = {mp.nstr(S,8)})")
        else:
            print(f"T_{m},{k} = PSLQ failed  (S = {mp.nstr(S,8)})")

# Also try: maybe A', B' need to be scaled for integrality.
# For lambda=1/4: A = sqrt(2e) K_0(1/2) and B = sqrt(2e) K_1(1/2).
# The "sqrt(2e)" factor arose naturally from the formula:
# U(1/2,1,1) = e^{1/2} K_0(1/2)/sqrt(pi), and h_0 = sqrt(2) pi^{3/2}.
# A = T_00 = c_00 = J_00/h_0 = 2pi^{3/2} sqrt(e) K_0(1/2) / (sqrt(2) pi^{3/2})
# = sqrt(2e) K_0(1/2). So A = sqrt(2e) K_0(1/2) naturally (it's c_00 itself).

# For lambda=3/4: "A'" = c_00 = 0.897...
# This is not a "nice" constant like sqrt(2e) K_0(1/2).
# But U(3/2,1,1) and U(5/2,2,1) are the natural generators.
# The denominators in the rational table are products of Pochhammer symbols:
# For T_02: 13/5 A' - 4 B'. Denominator = 5.
# For T_11: 11/3 A' - 20/3 B'. Denominator = 3.
# For T_22: 2203/75 A' - 728/15 B'. Denominator = 75 = 3*25.
# Pattern: denominators come from the Taylor coefficients of the GF.
# The Taylor coefficient [u^n] (1+iu)^{-3/2} = (-i)^n Gamma(n+3/2)/(n! Gamma(3/2)) = (-i)^n (3/2)_n/n!
# Inverse: n!/(3/2)_n. And (3/2)_n = 3/2 * 5/2 * ... * (2n+1)/2 = (2n+1)!!/(2^n).
# So the "natural" integer normalization absorbs these: multiply T_mk by (3/2)_m * (3/2)_k.

print("\n=== Testing: (3/2)_m * (3/2)_k * T_mk (with T from (-i)^{m+k} ... / pref) ===")
for m in range(6):
    for k in range(m, 6):
        Jmk = J(m, k)
        # T_raw = (-i)^{m+k} * beta_m * beta_k * J_mk / pref (real)
        # multiply by Poch_m * Poch_k / (beta_m * beta_k) = ((3/2)_m/m!)^2 * m!k! * ...
        # Actually: Poch_m * Poch_k * T_mk = Poch_m * Poch_k * beta_m * beta_k * (-i)^{m+k} * J_mk / pref
        # = Poch_m * (m!/Poch_m) * Poch_k * (k!/Poch_k) * (-i)^{m+k} * J_mk / pref
        # = m! * k! * (-i)^{m+k} * J_mk / pref  ... different from S above
        S_raw = mp.fac(m) * mp.fac(k) * (-1j)**(m+k) * Jmk / pref
        S = mp.re(S_raw)
        r = mp.pslq([S, Ap, Bp], tol=1e-35, maxcoeff=10000000)
        if r is not None:
            cA = mp.mpf(-r[1]) / mp.mpf(r[0])
            cB = mp.mpf(-r[2]) / mp.mpf(r[0])
            fA = Fraction(float(cA)).limit_denominator(100000)
            fB = Fraction(float(cB)).limit_denominator(100000)
            print(f"S_{m},{k} = {fA}*A' + {fB}*B'")
        else:
            print(f"S_{m},{k} = {mp.nstr(S, 8)} (PSLQ failed)")

# From my T_{01} computation: T_01 = A' - 2B' = U(3/2,1,1) - 2*U(5/2,2,1) = 0.35793 - 0.43065 = -0.0727.
# And (3/2)_0 * (3/2)_1 * T_01 = 1 * (3/2) * T_01 = 1.5 * (-0.0727) = -0.1091.
# vs m!k! (-i)^{m+k} J_mk / pref = 1*1 * (-i) * J_01 / pref
# J_01 = -0.5383i, so 1*1*(-i)*(-0.5383i) = (-i)(-0.5383i) = -0.5383*i^2 = 0.5383... wait:
# (-i)*(-0.5383 i) = (-i)*(-i)*0.5383 = i^2 * (-0.5383) = -(-0.5383) = 0.5383? No:
# (-i)*(-0.5383 i) = 0.5383 * i * i = 0.5383 * (-1) = -0.5383. Hmm.
# Wait: J_01 = c(0,1) * sqrt(h0*h1) = (numerically imaginary) * (real positive).
# c(0,1) = 2.101e-48 - 0.2232j (from output). So J_01 = -0.2232 * sqrt(h0*h1) * i ≈ -0.5383 i.
# (-i) * (-0.5383 i) = (-i)(-0.5383)(i) = 0.5383 * i^2 * (-1)...
# Actually (-i)*(-0.5383 i) = 0.5383 * ((-i)(i)) = 0.5383 * (-(i^2)) = 0.5383 * 1 = 0.5383.
# Then / pref = 0.5383 / (pi^2/2) = 0.5383 / 4.9348 = 0.1091 = 3/2 * (A' - 2B').
# So m!k! (-i)^{m+k} J_mk / pref = 0.1091 ≈ 3/2 * T_01. Makes sense:
# (3/2)_1 = 3/2, so this equals (3/2)_1 * T_01. ✓

# Let's look at whether m!k! (-i)^{m+k} J_mk / pref has rational A', B' coefficients.
print("\n--- Actually checking more carefully ---")
print("m!k! (-i)^{m+k} J_mk / pref = (3/2)_m (3/2)_k * T_mk")
print("Since T_mk = rational * A' + rational * B', these will also be rational * A' + B'.")
print("But the denominators from T times Pochhammer may cancel.\n")

for m in range(5):
    for k in range(m, 5):
        # Get the rational coeff from T_mk
        Jmk = J(m, k)
        beta_m = mp.fac(m) / Poch(m)
        beta_k = mp.fac(k) / Poch(k)
        T_raw = (-1j)**(m+k) * beta_m * beta_k * Jmk / pref
        T = mp.re(T_raw)
        r = mp.pslq([T, Ap, Bp], tol=1e-35, maxcoeff=100000000)
        if r is None:
            print(f"  T_{m},{k}: PSLQ failed (T = {mp.nstr(T,8)})")
            continue
        cA = Fraction(float(mp.mpf(-r[1]) / mp.mpf(r[0]))).limit_denominator(100000)
        cB = Fraction(float(mp.mpf(-r[2]) / mp.mpf(r[0]))).limit_denominator(100000)
        # T_mk = cA * A' + cB * B'
        # Poch_m * Poch_k * T_mk = Poch_m * Poch_k * (cA * A' + cB * B')
        pm = Fraction(float(Poch(m))).limit_denominator(1000)
        pk = Fraction(float(Poch(k))).limit_denominator(1000)
        print(f"  T_{m},{k}: T = {cA}*A' + {cB}*B'")
        print(f"    (3/2)_m = {pm}, (3/2)_k = {pk}")
        print(f"    Poch*Poch*T = {pm*pk*cA}*A' + {pm*pk*cB}*B'")

# Let me think differently: what are the "correct" A', B' for this problem?
# For lambda=1/4, the T_00 = A = sqrt(2e) K_0(1/2) was c_00 itself!
# And A = c_00, B = ? was identified via the GF.
# For lambda=3/4, c_00 = 0.897... = J_00/h_0 = pref*Ap/h_0.
# The question is: what normalization N makes N*c_mk "nice"?
# For l=1/4: N_m = m! * sqrt(C(2m,m)) = m! * sqrt(h_m^{1/4}/h_0^{1/4}) * 2^m.
# Equivalently: N_m^2 = (m!)^2 * C(2m,m) = (2m)!/m!.
# For l=3/4: N_m^2 = ? To make (N_m N_k c_mk)^2 integer.
# From h_n^{3/4} = h_0 * Poch_n / n!:
# N_m^{3/4} = m! / sqrt(Poch_m): try N_m = m! / (Poch_m)^{1/2} and T_mk = N_m N_k c_mk.
# T_mk = m! k! * J_mk / sqrt(Poch_m * Poch_k) / sqrt(h_m h_k)
# = m! k! / sqrt(Poch_m * Poch_k * h_m h_k) * J_mk
# h_m = h_0 * Poch_m / m!, so sqrt(h_m h_k) = h_0 * sqrt(Poch_m Poch_k / (m! k!))
# T_mk = m! k! / sqrt(Poch_m Poch_k) / (h_0 * sqrt(Poch_m Poch_k/(m!k!))) * J_mk
# = m! k! / (Poch_m Poch_k) * sqrt(m! k!) / h_0 * J_mk ... messy.

# Let me just try ALL normalizations of the form T_mk = m!^a k!^b Poch_m^c Poch_k^d * c_mk
# for a,b,c,d = 0,1,2 and see which gives integer coefficients.
print("\n=== Brute force: testing normalizations N(m)*N(k)*c_mk ===")
# Try N(n) = n! / (3/2)_n^p for p = 0, 1/2, 1, 3/2, 2
import itertools
best = None
for p in [0, 0.5, 1, 1.5, 2]:
    def N(n, p=p):
        return mp.fac(n) / Poch(n)**mp.mpf(p)
    # Check T_00: N(0)^2 * c_00 = 1^2 * c_00 = c_00 (always = 1 generator)
    # Check T_01: N(0)*N(1)*c_01
    T01 = N(0)*N(1) * c(0,1)  # purely imaginary
    T01_val = mp.im(T01)  # extract imaginary part
    # Check T_11: N(1)^2 * c_11 (real)
    T11 = mp.re(N(1)**2 * c(1,1))
    # Check T_02: N(0)*N(2)*c_02 (real)
    T02 = mp.re(N(0)*N(2)*c(0,2))
    r01 = mp.pslq([T01_val, Ap, Bp], tol=1e-30, maxcoeff=10000000)
    r11 = mp.pslq([T11, Ap, Bp], tol=1e-30, maxcoeff=10000000)
    r02 = mp.pslq([T02, Ap, Bp], tol=1e-30, maxcoeff=10000000)
    ok01 = r01 is not None and abs(r01[0]) == 1
    ok11 = r11 is not None and abs(r11[0]) == 1
    ok02 = r02 is not None and abs(r02[0]) == 1
    print(f"p={p}: T01 {r01}, T11 {r11}, T02 {r02}")
