"""Fix phase issue and find the right normalization with rational denominators.

The issue: T_mk definition was using (-i)^{m+k} but the phase law says c_mk = (-i)^{|m-k|} |c_mk|.
For m+k odd, m≠k: |m-k| is odd, so the phase is (-i)^{|m-k|} (purely imaginary).
But for m=k (diagonal), all entries are real regardless.

Let me redo: work directly with J_mk.

From the GF derivation:
  sum J_mk w^m u^k = pi^2/2 * U(3/2, 1, q) / [(1+iw)(1+iu)]^{3/2}
  where q = (1-iw)(1-iu)/((1+iw)(1+iu)) (VERIFIED: GF(0,0) = J_00 ✓, [u^1] U(q)|_0 = 3i*B' ✓)

But: T_{0,1} should come from [w^0 u^1] of GF:
  [u^1] pi^2/2 * U(3/2,1,q) / [(1+iu)^{3/2}] |_{w=0}
  = pi^2/2 * d/du [U(3/2,1,(1-iu)/(1+iu)) / (1+iu)^{3/2}]|_{u=0}
  = pi^2/2 * [U'(1)*(-2i)/(1+iu)^{3/2} + U(1)*(-3/2*i)/(1+iu)^{5/2}]|_{u=0}
  = pi^2/2 * [(-3/2 B')(-2i) - (3/2 i) A']
  = pi^2/2 * [3i B' - (3/2) i A']
  = pi^2/2 * 3i/2 * [2B' - A'/1] ... hmm
  = pi^2/2 * i * (3 B' - 3/2 A')

So J_{0,1} = pi^2/2 * i * (3 B' - 3/2 A').
But J_{0,1} should be purely imaginary (m+k=1 odd, phase = -i^1 = -i times real).
J_{0,1} = (-i) * real, so pi^2/2 * i * (3B' - 3/2 A') = (-i) * (something real).
=> 3B' - 3/2 A' should be real. And it is: 3*0.2153 - 3/2*0.3579 = 0.6459 - 0.5369 = 0.1090. ✓

So J_{0,1} = -pi^2/2 * i * (3/2 A' - 3 B')... let me compute numerically.
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

# Verify J_01 prediction
pref = mp.pi**2 / 2
J01_pred = pref * 1j * (3*Bp - mp.mpf(3)/2*Ap)
J01_num = J(0, 1)
print(f"J_01 prediction: {mp.nstr(J01_pred, 12)}")
print(f"J_01 numerical: {mp.nstr(J01_num, 12)}")
print(f"Diff: {mp.nstr(abs(J01_pred - J01_num), 5)}")

# So J_01 = i * pi^2/2 * (3B' - 3/2 A'). Phase = i = (-i)^{-1}...
# Phase law: J_mk = J_mk_real * phase where phase = (-i)^{|m-k|} (for |m-k| odd: phase is -i or i).
# |m-k| = 1: phase = (-i)^1 = -i. So J_01 = (-i) * R for some real R.
# J_01 = pref * i * (3B' - 3/2 A') = -(-i) * pref * (3B' - 3/2 A')
# So R_{0,1} = -pref * (3B' - 3/2 A') = pref * (3/2 A' - 3 B').
# Normalized: R_{0,1} / pref = 3/2 A' - 3 B'.

# The right way to extract real entries:
# Let S_mk = J_mk / ((-i)^{|m-k|} * pref * beta_m * beta_k)
# Wait, we need the Pochhammer factors too. From the GF:
# J_mk = pi^2/2 * [w^m u^k] U(q) / [(1+iw)(1+iu)]^{3/2}
# [w^m u^k] 1/(1+iw)^{3/2} = (-i)^m (3/2)_m / m!
# So the "stripped" coefficient is:
# J_mk * (m! / ((-i)^m (3/2)_m)) * (k! / ((-i)^k (3/2)_k)) / pref
# = J_mk * m! k! / ((-i)^{m+k} (3/2)_m (3/2)_k) / pref
# But this is [w^m u^k] U(q) = T_mk (what we computed before, which was ~0 for odd entries).

# Ah, I see the bug: for odd m+k, the T_mk is purely imaginary, not real.
# I was extracting im(T_mk) but the issue is that J_mk for m+k odd is imaginary
# and beta factors are real, so T_mk = (-i)^{m+k} * (real) * J_mk / pref
# For m+k=1: T_mk = (-i)^1 * (real) * J_{0,1} / pref = (-i) * (real) * J_{0,1} / pref
# J_{0,1} is purely imaginary, so (-i) * imaginary = (-i)*(i*r) = r (real!)
# So T_mk IS real for all (m,k). Let me recheck.

# T_01 = (-i)^1 * beta(0)*beta(1) * J(0,1) / pref
# = (-i) * 1 * 1 * J(0,1) / pref
beta0 = mp.fac(0) * mp.gamma(mp.mpf(3)/2) / mp.gamma(mp.mpf(3)/2)  # = 1
beta1 = mp.fac(1) * mp.gamma(mp.mpf(3)/2) / mp.gamma(mp.mpf(3)/2 + 1)
# gamma(3/2) = sqrt(pi)/2, gamma(5/2) = 3sqrt(pi)/4, so beta(1) = 1 * (sqrt(pi)/2) / (3sqrt(pi)/4) = 2/3
print(f"\nbeta(0) = {mp.nstr(beta0, 6)}")
print(f"beta(1) = {mp.nstr(beta1, 6)} (should be 2/3 = 0.6667)")

T01_raw = (-1j)**1 * beta0 * beta1 * J(0,1) / pref
print(f"T_01_raw = {mp.nstr(T01_raw, 12)}")
# J(0,1) is purely imaginary
print(f"J(0,1) = {mp.nstr(J(0,1), 12)}")
# (-i) * (i * r) = r (real). So T_01 should be real!
# From J_01 ~ i * pref * (3B' - 3/2 A') * (1/beta0/beta1):
# Wait, J_01 = pref * i * (3B' - 3/2 A') (from our derivation).
# But beta factors in T are outside:
# T_01 = (-i) * beta(0)*beta(1) * J(0,1) / pref = beta(1) * (-i*i) * (3B' - 3/2 A') * beta(0)
# = beta(0)*beta(1) * (3B' - 3/2 A') = 1 * 2/3 * (3B' - 3/2 A') = 2B' - A'.

# Let me check: J(0,1) was giving ~0 in the T computation above, but it shouldn't.
# The issue was in the code: T_raw = (-1j)^{m+k} not (-1j)^{|m-k|}.
# For m=0,k=1: m+k=1, (-i)^1 = -i, J(0,1) = imaginary. OK, that should give real.
# But the output said T_0,1 = -6.8e-49 (zero). That means J(0,1) = 0 numerically?!
# Look at the raw data: c_0,1 = 2.101e-48 - 0.2232j
# So c(0,1) is IMAGINARY (m+k=1, imaginary entry). J(0,1) = c(0,1)*sqrt(h0*h1).
# c(0,1) = (2.101e-48 - 0.2232j): the imaginary part is -0.2232.
# T_01 = (-i) * 1 * 2/3 * J(0,1) / pref

# Let's check carefully
J01 = J(0, 1)
print(f"\nJ(0,1) = {mp.nstr(J01, 12)}")
T01 = (-1j) * beta0 * beta1 * J01 / pref
print(f"T_01 = {mp.nstr(T01, 12)}")
print(f"Expected: 2B' - A' = {mp.nstr(2*Bp - Ap, 12)}")

# Hmm, but earlier the code gave T_0,1 = -6.8e-49 ≈ 0. That can't be right.
# Let me look at the raw c_01 value from the CACHE.
print(f"\nc(0,1) from cache = {mp.nstr(C[(0,1)], 12)}")
print(f"nrm[0] = {mp.nstr(nrm[0], 12)}")
print(f"nrm[1] = {mp.nstr(nrm[1], 12)}")
print(f"sqrt(h0*h1) = {mp.nstr(mp.sqrt(nrm[0]*nrm[1]), 12)}")

# beta(1): 1! * Gamma(3/2) / Gamma(1 + 3/2) = Gamma(3/2)/Gamma(5/2) = (sqrt(pi)/2)/(3sqrt(pi)/4) = 2/3
beta1_check = mp.fac(1) * mp.gamma(mp.mpf(3)/2) / mp.gamma(1 + mp.mpf(3)/2)
print(f"beta(1) = {mp.nstr(beta1_check, 6)}")

# Now: the issue in mp34_integrality.py was in the loop:
# for m in range(7): for k in range(m, 7): T_raw = (-1j)^{m+k} * beta(m)*beta(k)*J_mk/pref
# This DOES match T_01 = (-i)^1 * ... So why was T_01 printing ~0?
# ANSWER: because in the loop, I was taking re(T_raw) for m+k even and im(T_raw) for m+k odd!
# T_01 = (-i)^{1} * (real)*J_{0,1}/pref = (-i)*imaginary = real!
# So T_01 IS real. But the loop code says:
# "if (m+k) % 2 == 0: val = re(T_raw) else: val = im(T_raw)"
# For m+k=1: val = im(T_raw). But T_01 is REAL, so im(T_01) ≈ 0. Bug!

# The correct thing: T_mk = re(T_raw) for ALL (m,k) since T_raw should always be real.
# Let me verify.
print("\n=== Verifying T_mk = re((-i)^{m+k} beta_m beta_k J_mk / pref) for all m,k ===")
def beta(n):
    return mp.fac(n) * mp.gamma(mp.mpf(3)/2) / mp.gamma(n + mp.mpf(3)/2)

all_real = True
for m in range(6):
    for k in range(m, 6):
        Jmk = J(m, k)
        T_raw = (-1j)**(m+k) * beta(m) * beta(k) * Jmk / pref
        re_part = mp.re(T_raw)
        im_part = mp.im(T_raw)
        if abs(im_part) > 1e-30 * abs(re_part):
            print(f"  T_{m},{k}: NOT REAL! re={mp.nstr(re_part,5)}, im={mp.nstr(im_part,5)}")
            all_real = False
        else:
            pass  # good
if all_real:
    print("All T_mk are real (as expected)")

# Now compute the correct T table and run PSLQ
print("\n=== CORRECT T_mk values ===")
T_table = {}
for m in range(7):
    for k in range(m, 7):
        Jmk = J(m, k)
        T_raw = (-1j)**(m+k) * beta(m) * beta(k) * Jmk / pref
        T = mp.re(T_raw)
        T_table[(m, k)] = T

# Find rational coefficients A', B' for each entry
print("Searching for rational A', B' coefficients (T_mk = p*A' + q*B'):")
rational_table = {}
for (m, k), T in T_table.items():
    if m > 5 or k > 5:
        continue
    r = mp.pslq([T, Ap, Bp], tol=1e-35, maxcoeff=10000000)
    if r is not None:
        coeff_A = mp.mpf(-r[1]) / mp.mpf(r[0])
        coeff_B = mp.mpf(-r[2]) / mp.mpf(r[0])
        rational_table[(m,k)] = (coeff_A, coeff_B)
        print(f"  T_{m},{k} = {mp.nstr(T, 8)}: {mp.nstr(coeff_A, 6)}*A' + {mp.nstr(coeff_B, 6)}*B'")
    else:
        print(f"  T_{m},{k} = {mp.nstr(T, 8)}: PSLQ failed")

# Find common denominator structure
print("\n=== Checking rational structure ===")
print("Denominators appearing in coefficients:")
from fractions import Fraction
frac_table = {}
for (m, k), (coeff_A, coeff_B) in rational_table.items():
    # Round to fractions
    cA = float(coeff_A)
    cB = float(coeff_B)
    fA = Fraction(cA).limit_denominator(10000)
    fB = Fraction(cB).limit_denominator(10000)
    frac_table[(m,k)] = (fA, fB)
    print(f"  T_{m},{k}: coeff_A = {fA}, coeff_B = {fB}")

# Find the LCM of denominators and multiply to get integers
from math import gcd
def lcm(a, b):
    return a * b // gcd(a, b)

all_denoms = []
for (m,k), (fA, fB) in frac_table.items():
    all_denoms.append(fA.denominator)
    all_denoms.append(fB.denominator)

from functools import reduce
L = reduce(lcm, all_denoms, 1)
print(f"\nLCM of all denominators: {L}")

# Multiply T by L to get integer generators
print("\n=== Integer table: L * T_mk = p*A'' + q*B'' where A'' = L*A', B'' = L*B' ===")
# Actually the integrality should be: S_mk = T_mk / (some factor) ∈ Z*A' + Z*B'
# OR: find the right scaling of A', B' to be the "true" generators.
# The denominators 15 in the output suggest: the right generators are
# A'' = 15*A' = 15*U(3/2,1,1), B'' = 15*B' = 15*U(5/2,2,1)?
# Or maybe: A' and B' should be replaced by linear combinations with integer coefficients.

print("\nLet's see if there's a better basis. The rational coefficients are:")
for (m,k), (fA, fB) in frac_table.items():
    if m <= 4 and k <= 4:
        print(f"  T_{m},{k} = {fA}*A' + {fB}*B'")

# Maybe the natural generators for lambda=3/4 are:
# A'' = U(3/2,1,1) and B'' = U(5/2,2,1) but scaled:
# Note: from the GF structure, the generating series is
# sum [beta_m beta_k J_mk / pref] * (-i)^{m+k} w^m u^k = U(q)
# The Taylor coefficients of U(q) at q=1 are in Q*A' + Q*B'.
# Let me compute those directly from the ODE.
# The question is: what's the common denominator structure?

# From the calculation:
# T_00 = A'             (integer coefficients of A', B')
# T_01 = 2B' - A'       (see prediction above)
# T_02 should follow from differentiating the GF twice

# Let me compute T_01 directly
print(f"\nT_01 = (-i)^1 * beta(0)*beta(1) * J(0,1) / pref = {mp.nstr(T_table.get((0,1), 0), 12)}")
print(f"2*B' - A' = {mp.nstr(2*Bp - Ap, 12)}")

# Check T_02
T02 = T_table.get((0,2), None)
print(f"T_02 = {mp.nstr(T02, 12)}")
print(f"Prediction -3A' + 3B' = {mp.nstr(-3*Ap + 3*Bp, 12)}")
# From our earlier calculation: [w^0 u^2] U(q) = -3A' + 3B'.
# But T_02 = (-i)^2 * beta(0)*beta(2) * J(0,2) / pref.
# And [u^2] U(q) = (1/2!) d^2/du^2 U(q)|_0 ... no wait, Taylor coefficients.
# [u^2] in the GF means J_{0,2}:
# J_{0,2} = pi^2/2 * [u^2] U(q) / [(1+iu)^{3/2}]
# [u^2] f(u) = (1/2) f''(0) (second Taylor coefficient = 1/2 * second derivative)
# Wait no: [u^2] means the coefficient of u^2 in the power series, which is f''(0)/2!.
# So J_{0,2} = pi^2/2 * [coefficient of u^2 in U(q)/(1+iu)^{3/2}]
# NOT [u^2] U(q) alone.
# T_{0,2} = (-i)^2 * beta(0)*beta(2) * J_{0,2} / pref
# = 1 * 1 * beta(2) * J_{0,2} / pref
# beta(2) = 2! * Gamma(3/2)/Gamma(7/2) = 2 * (sqrt(pi)/2)/(15sqrt(pi)/8) = 2 * 4/15 = 8/15
print(f"beta(2) = {mp.nstr(beta(2), 6)} (expected 8/15 = {8/15:.6f})")

# So T_{0,2} = (8/15) * J_{0,2} / pref.
# From GF, J_{0,2} = pref * [u^2] [U(q)/(1+iu)^{3/2}] at w=0
# = pref * (coeff of u^2 in U((1-iu)/(1+iu))/(1+iu)^{3/2})
# Let g(u) = U((1-iu)/(1+iu)) / (1+iu)^{3/2}.
# [u^2] g = g''(0)/2.
# We need g''(0). This involves more computation.
# From Fourier extraction: [u^2] of GF/pref = [u^2] U(q)/(1+iu)^{3/2}
# Let's just trust the numerical PSLQ result for the correct T table.

# The correct T_01: let's just compute
print(f"\n--- Correct T_01 computation ---")
J01_val = J(0, 1)
T01_val = (-1j)**1 * beta(0) * beta(1) * J01_val / pref
print(f"T_01 = {mp.nstr(T01_val, 15)}")
r01 = mp.pslq([mp.re(T01_val), Ap, Bp], tol=1e-35, maxcoeff=1000000)
print(f"PSLQ: {r01}")

# Full correct table, using re(T_raw):
print("\n=== FINAL CORRECT TABLE ===")
for m in range(6):
    for k in range(m, 6):
        T = T_table[(m,k)]
        r = mp.pslq([T, Ap, Bp], tol=1e-35, maxcoeff=100000000)
        if r is not None:
            cA = mp.mpf(-r[1]) / mp.mpf(r[0])
            cB = mp.mpf(-r[2]) / mp.mpf(r[0])
            fA = Fraction(float(cA)).limit_denominator(10000)
            fB = Fraction(float(cB)).limit_denominator(10000)
            print(f"T_{m},{k} = {fA}*A' + {fB}*B'   (T = {mp.nstr(T,8)})")
        else:
            print(f"T_{m},{k} = PSLQ failed  (T = {mp.nstr(T,8)})")
