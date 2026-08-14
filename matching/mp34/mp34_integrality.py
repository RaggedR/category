"""Find integer table for lambda=3/4 matching matrix.

Key results from analysis:
- GF: sum J_mk w^m u^k = pi^2/2 * U(3/2, 1, q) / [(1+iw)(1+iu)]^{3/2}
  where q = (1-iw)(1-iu)/((1+iw)(1+iu))
- [w^0 u^1] U(q) = 3i*B', so U coefficients ARE in Q*A' + Q*B'.
- The right normalization is: T_mk = (-i)^{m+k} * beta_m * beta_k * J_mk / (pi^2/2)
  where beta_n = n! Gamma(3/2) / Gamma(n+3/2) = 1/(3/2)_n * n!
  so that T_mk = [w^m u^k] U(q)  and T_mk in Q*A' + Q*B'.
- A' = U(3/2,1,1), B' = U(5/2,2,1) = -2/3 * d/dz U(3/2,1,z)|_{z=1}

But we want INTEGER generators. For lambda=1/4:
  A = sqrt(2e) K_0(1/2), B = sqrt(2e) K_1(1/2) (nice closed forms)
For lambda=3/4: what are the "nice" expressions?

From the ODE z y'' + (1-z) y' - 3/2 y = 0:
- A' = U(3/2,1,1) = 0.35793...
- B' = U(5/2,2,1) = 0.21532...
- B' = 2/3 * (-y'(1)) where y'(1) = d/dz U(3/2,1,z)|_1 = -3/2 B'... wait
  Actually d/dz U(a,b,z) = -a U(a+1,b+1,z), so d/dz U(3/2,1,z) = -3/2 U(5/2,2,1) = -3/2 B'
  So B' = -2/3 y'(1). Not a separate "nice" constant -- just the derivative.

For lambda=1/4:
  U(1/2,1,2z) = e^z K_0(z)/sqrt(pi)
  U(1/2,1,1) = e^{1/2} K_0(1/2)/sqrt(pi) = sqrt(e/pi) K_0(1/2)
  A = sqrt(2e) K_0(1/2): so A = sqrt(2pi) * U(1/2,1,1).
  B = sqrt(2e) K_1(1/2): and U(3/2,2,1) = d/dz U(1/2,1,z)/(-1/2)|_1...
  d/dz U(1/2,1,z)|_1 = -1/2 U(3/2,2,1).
  U(1/2,1,z) = e^{z/2} K_0(z/2)/sqrt(pi) => d/dz = e^{z/2}/2 (K_0 - K_1)(z/2)/sqrt(pi)
  At z=1: d/dz = e^{1/2}/2 * (K_0(1/2)-K_1(1/2))/sqrt(pi).
  So U(3/2,2,1) = -2*(-1/2)*e^{1/2}*(K_0(1/2)-K_1(1/2))/(2*sqrt(pi))
  = e^{1/2}*(K_0(1/2)-K_1(1/2))/sqrt(pi)... involves both K_0 and K_1.
  The two generators for lambda=1/4 are really: A = sqrt(2pi) U(1/2,1,1) and B comes from B+A or B-A.

For lambda=3/4, the "nice" form of U(3/2,1,z):
  Check if U(3/2,1,z) = sqrt(pi) e^{z/2} * [a K_0(z/2) + b K_1(z/2)] ...?
  At z=1: U(3/2,1,1) = sqrt(pi) e^{1/2} * [a K_0(1/2) + b K_1(1/2)]
  At z=2: U(3/2,1,2) numerically and K_0(1), K_1(1)...
  This is a 2-parameter family -- can we match?
"""
import mpmath as mp
import pickle

mp.mp.dps = 50

with open("mp34_dual_matching.pkl", "rb") as f:
    nrm, C = pickle.load(f)

def c(m, k):
    return C[(min(m,k), max(m,k))]

def J(m, k):
    return c(m, k) * mp.sqrt(nrm[m] * nrm[k])

h0 = nrm[0]
J00 = J(0, 0)

# Generators
Ap = mp.hyperu(mp.mpf(3)/2, 1, 1)   # A' = U(3/2,1,1)
Bp = mp.hyperu(mp.mpf(5)/2, 2, 1)   # B' = U(5/2,2,1)

# Pochhammer: beta_n = n! / (3/2)_n
def beta(n):
    return mp.fac(n) / (mp.gamma(n + mp.mpf(3)/2) / mp.gamma(mp.mpf(3)/2))

pref = mp.pi**2 / 2   # pi^2/2

# T_mk = (-i)^{m+k} * beta_m * beta_k * J_mk / (pi^2/2) = [w^m u^k] U(q)
print("=== T_mk = (-i)^{m+k} beta_m beta_k J_mk / (pi^2/2) ===")
print("=== Testing PSLQ vs A' = U(3/2,1,1) and B' = U(5/2,2,1) ===\n")

T = {}
for m in range(7):
    for k in range(m, 7):
        Jmk = J(m, k)
        T_raw = (-1j)**(m+k) * beta(m) * beta(k) * Jmk / pref
        T_re = mp.re(T_raw)
        T_im = mp.im(T_raw)
        # Determine which component
        if (m+k) % 2 == 0:
            val = T_re
        else:
            val = T_im  # T_raw = i * (real number)
        T[(m, k)] = val
        r = mp.pslq([val, Ap, Bp], tol=1e-30, maxcoeff=100000)
        print(f"  T_{m},{k} = {mp.nstr(val, 12)},  PSLQ: {r}")

# Now show the table in terms of A', B'
print("\n=== Integer table T_mk in Z*A' + Z*B' ===")
pslq_results = {}
for (m, k), val in T.items():
    r = mp.pslq([val, Ap, Bp], tol=1e-30, maxcoeff=1000000)
    if r is not None:
        # r[0]*val + r[1]*Ap + r[2]*Bp = 0 => val = (-r[1]*Ap - r[2]*Bp)/r[0]
        alpha = -r[1] * mp.mpf(r[0])  # coefficient of Ap per unit
        beta_c = -r[2] * mp.mpf(r[0])  # coefficient of Bp per unit
        pslq_results[(m, k)] = (r[1], r[2], r[0])  # val = -(r1/r0)*A' - (r2/r0)*B'
        coeff_A = mp.mpf(-r[1]) / mp.mpf(r[0])
        coeff_B = mp.mpf(-r[2]) / mp.mpf(r[0])
        print(f"  T_{m},{k} = {mp.nstr(coeff_A,4)}*A' + {mp.nstr(coeff_B,4)}*B'")
    else:
        print(f"  T_{m},{k} = {mp.nstr(val, 8)}  (PSLQ FAILED)")

# Now try to identify what the right INTEGERS are.
# From the GF analysis: [w^0 u^1] U(q) = 3i B', so T_{0,1}/i = 3 B'.
# But T_{0,1} as defined is im((-i)^1 * beta_0 * beta_1 * J_{0,1} / pref).
# Let's check.
print(f"\nDirect check T_0,1:")
J01 = J(0, 1)
T01_raw = (-1j)**1 * beta(0) * beta(1) * J01 / pref
print(f"  T_01_raw = {mp.nstr(T01_raw, 12)}")
print(f"  Predicted = 3i * B' = {mp.nstr(3j * Bp, 12)}")
print(f"  Match? diff = {mp.nstr(abs(T01_raw - 3j*Bp), 5)}")

# Look at U(q) Taylor coefficients extracted via differentiation
# U(q) = U(3/2,1,q) as function of w (with u=0):
# q|_{u=0} = (1-iw)/(1+iw)
# dq/dw = -2i/(1+iw)^2
# At w=0: q=1, dq/dw = -2i
# d/dw U(q) = U'(q) dq/dw = -3/2 U(5/2,2,q) * (-2i) = 3i U(5/2,2,q)
# At w=0: = 3i * U(5/2,2,1) = 3i B'. ✓

# Similarly: d^2/dw^2 U(q)|_{w=0,u=0}:
# d^2/dw^2 = U''(q)(dq/dw)^2 + U'(q) d^2q/dw^2
# d^2q/dw^2 = 4i^2/(1+iw)^3 * (-i) + ... = d/dw[-2i/(1+iw)^2] = -2i*(-2)(i)/(1+iw)^3 = -4/(1+iw)^3
# At w=0: d^2q/dw^2 = -4 * (-i) = wait: d/dw[-2i/(1+iw)^2] = -2i * (-2)(1+iw)^{-3} * i = 4i^2 (1+iw)^{-3}...
# = -4(1+iw)^{-3}. At w=0: = -4.
# U''(q) = d/dz U'(z) = d/dz (-3/2 U(5/2,2,z)) = (-3/2)(-5/2) U(7/2,3,z) = 15/4 U(7/2,3,z)
# So d^2/dw^2 U(q)|_0 = 15/4 * U(7/2,3,1) * (-2i)^2 + (-3/2 B') * (-4)
# = 15/4 * U(7/2,3,1) * (-4) + 6 B'
# = -15 U(7/2,3,1) + 6 B'

# But from Fourier: [w^0 u^2] U(q) = -0.4278...
# Let me check: U(7/2,3,1) and 6B' - 15*U(7/2,3,1)
U72 = mp.hyperu(mp.mpf(7)/2, 3, 1)
val_02 = -15*U72 + 6*Bp
print(f"\n[w^0 u^2] U(q):")
print(f"  Predicted: -15*U(7/2,3,1) + 6*B' = {mp.nstr(val_02, 12)}")
# From Fourier above: -0.4278193128
print(f"  But from Fourier: -0.42782...")

# U(7/2,3,1):
print(f"  U(7/2,3,1) = {mp.nstr(U72, 12)}")
# Now: is U(7/2,3,1) in Q*A' + Q*B'?
r_U72 = mp.pslq([U72, Ap, Bp], tol=1e-30, maxcoeff=1000000)
print(f"  PSLQ U(7/2,3,1) vs A',B': {r_U72}")
# The contiguous recurrences reduce all U(a+n, b+n, 1) to A', B' via the ODE.
# U(a+1,b+1,z) = -(1/a) d/dz U(a,b,z)
# So U(5/2,2,z) = -(1/(3/2)) d/dz U(3/2,1,z) = -2/3 * (-3/2) U(5/2,2,z) ... circular.
# What reduces U(7/2,3,1)?
# Use ODE at z=1: 1*U'' + (1-1)U' - 3/2 U = 0 => U''(1) = 3/2 U(1) = 3/2 A'
# And U''(z) = U''(q)|... wait, the ODE is for U(3/2,1,z).
# U(7/2,3,z) = (1/15) U''(3/2,1,z) where prime is d/dz.
# From ODE: z U'' + (1-z) U' - 3/2 U = 0
# At z=1: U''(1) + 0 - 3/2 U(1) = 0 => U''(1) = 3/2 A'.
# So U(7/2,3,1) = (1/15) * 3/2 * A' = A'/10.
# Wait: d/dz U(3/2,1,z) = -3/2 U(5/2,2,z), so U' = -3/2 B' at z=1.
# d^2/dz^2 U(3/2,1,z) = (-3/2) d/dz U(5/2,2,z) = (-3/2)(-5/2) U(7/2,3,z) = 15/4 U(7/2,3,z)
# So U(7/2,3,1) = 4/15 * U''(1) = 4/15 * 3/2 A' = 2A'/5.
r_U72_pred = mp.mpf(2)/5 * Ap
print(f"  Predicted U(7/2,3,1) = 2/5 A' = {mp.nstr(r_U72_pred, 12)}")
# From ODE at z=1: U''(1) = 3/2 A'. And U'' = 15/4 U(7/2,3,z) => U(7/2,3,1) = 4/(15) * 3/2 A' = 2/5 A'. ✓

# So [w^0 u^2] U(q) = -15 * 2/5 A' + 6 B' = -6 A' + 6 B'.
# Predicted: -6*0.35793... + 6*0.21532... = -2.1476 + 1.2920 = -0.8556... ≠ -0.4278.
# Hmm, that's off by factor 2. Let me recompute d^2/du^2 U(q)|_{u=0,w=0}.
# Wait, I'm computing partial derivatives, not Taylor coefficients.
# [u^2] means coefficient of u^2, which is (1/2!) * d^2/du^2 U(q)|_{u=0}.
# So [w^0 u^2] = (1/2) * (d^2/du^2 U(q)|_{u=0})
# = (1/2) * (-15 U(7/2,3,1)*(-4) + (-3/2 B')*(-4))
# Wait let me redo. At w=0: q = (1-iu)/(1+iu), same as q_w=0 function.
# d/du q|_{u=0}: q = (1-iu)/(1+iu), dq/du = (-i(1+iu)-(1-iu)i)/(1+iu)^2 = -2i at u=0.
# d^2/du^2 q|_{u=0}: d/du[-2i/(1+iu)^2] = -2i * (-2)(i)/(1+iu)^3 = 4i^2/(1+iu)^3 = -4 at u=0.
# d^2/du^2 U(q)|_0 = U''(1)(dq/du)^2 + U'(1) d^2q/du^2
# = (15/4 B' ... wait, U'' here means d^2U/dq^2 not d^2U/dz^2 of the ODE).
# U'(q) = dU/dq = -3/2 B' at q=1.
# U''(q) = d^2U/dq^2 = 15/4 U(7/2,3,q), at q=1 = 15/4 * 2/5 A' = 3/2 A'.
# So d^2/du^2 U(q)|_0 = (3/2 A') * (-2i)^2 + (-3/2 B') * (-4)
# = (3/2 A') * (-4) + 6 B' = -6 A' + 6 B'
# [u^2] = (1/2) * (-6A' + 6B') = -3A' + 3B'
# Check: -3*0.35793 + 3*0.21532 = -1.0738 + 0.6460 = -0.4278. ✓ YES!

print(f"\n[w^0 u^2] prediction = -3*A' + 3*B' = {mp.nstr(-3*Ap + 3*Bp, 12)}")
print(f"From Fourier: -0.4278193128")

# Great! Now the pattern is clear.
# T_mk = [w^m u^k] U(q) = rational_combination of A' and B'.
# The integers (after removing rational factors) will appear.

# Let me compute the table systematically using differentiation.
# The coefficients [w^m u^k] U(q) can be computed via symbolic differentiation.
# Let me compute up to order 6 using mpmath.
print("\n=== Computing Taylor table of U(q) symbolically ===")
# T[m,k] = coefficient of w^m u^k in U((1-iw)(1-iu)/((1+iw)(1+iu)))
# = rational combination of A' = U(3/2,1,1) and B' = U(5/2,2,1)

# We want to write T[m,k] = p[m,k]*A' + q[m,k]*B' where p,q are rationals.
# Use: U^{(n)}(1) = d^n/dz^n U(3/2,1,z)|_1 are in Q*A' + Q*B' by the ODE.
# d^n/dz^n U = (-1)^n * (3/2)_n U(3/2+n, 1+n, z)
# = (-1)^n * Gamma(3/2+n)/Gamma(3/2) * U(3/2+n, 1+n, z)
# All reduce to A', B' via ODE: at z=1, higher U values = rational * A' or B'.
# Specifically:
# U(3/2,1,1) = A'
# U(5/2,2,1) = B'
# U(7/2,3,1) = 2/5 A'  (from ODE: y''(1) = 3/2 y(1))
# U(9/2,4,1) = ?
# From ODE: z y'' + (1-z) y' - 3/2 y = 0, differentiate:
# z y''' + y'' + (1-z) y'' - y' - 3/2 y' = 0 => z y''' + (2-z) y'' - (1+3/2) y' = 0
# At z=1: y'''(1) + y''(1) - 5/2 y'(1) = 0 => y'''(1) = 5/2 y'(1) - y''(1)
# y'(1) = -3/2 B', y''(1) = 3/2 A'
# y'''(1) = 5/2*(-3/2 B') - 3/2 A' = -15/4 B' - 3/2 A'
# U(9/2,4,z) = (1/(-1)^3 * (3/2)(5/2)(7/2)) d^3/dz^3 U(3/2,1,z) ... wait
# d^n/dz^n U(3/2,1,z) = (-1)^n (3/2)_n U(3/2+n, 1+n, z)
# So U(9/2,4,z) = d^3/dz^3 U(3/2,1,z) / ((-1)^3 * (3/2)_3)
# (3/2)_3 = 3/2 * 5/2 * 7/2 = 105/8
# d^3/dz^3 U(3/2,1,z)|_1 = y'''(1) = -15/4 B' - 3/2 A'
# U(9/2,4,1) = (-15/4 B' - 3/2 A') / (-105/8) = (15/4 B' + 3/2 A') * 8/105 = (2/7 B' + 4/35 A') ... let me compute
U94 = mp.hyperu(mp.mpf(9)/2, 4, 1)
print(f"U(9/2,4,1) = {mp.nstr(U94, 12)}")
pred_U94 = (-mp.mpf(15)/4 * Bp - mp.mpf(3)/2 * Ap) / (-mp.mpf(105)/8)
print(f"Predicted: {mp.nstr(pred_U94, 12)}")

# Now let me compute the full Taylor table numerically and verify vs PSLQ
print("\n=== Full PSLQ table ===")
for (m, k), val in T.items():
    if m <= 5 and k <= 5:
        r = mp.pslq([val, Ap, Bp], tol=1e-35, maxcoeff=10000000)
        if r is not None:
            coeff_A = mp.mpf(-r[1]) / mp.mpf(r[0])
            coeff_B = mp.mpf(-r[2]) / mp.mpf(r[0])
            print(f"  T_{m},{k} = {mp.nstr(coeff_A)} * A' + {mp.nstr(coeff_B)} * B'")
        else:
            print(f"  T_{m},{k} = {mp.nstr(val, 10)} (no PSLQ)")

# Try: express via the "integer" normalization.
# For lambda=1/4: T_mk/h_0 in ZA+ZB, where A,B = sqrt(2e)*K_{0,1}(1/2).
# The integers were "nice" because A = sqrt(2pi)*U(1/2,1,1) and B involved K_1(1/2).
# For lambda=3/4: try T_mk with a scaling of A', B' to produce integers.
# From T_{0,1} = 3*B': this is already integer * B'. Good.
# From T_{0,0} = A': integer * A'.
# From T_{0,2} = -3A' + 3B': integers.
# The pattern suggests T_mk might already be in Z*A' + Z*B'!
print("\n=== Checking if T_mk are already in Z*A' + Z*B' (no extra factors) ===")
all_ints = True
int_table = {}
for (m, k), val in T.items():
    r = mp.pslq([val, Ap, Bp], tol=1e-35, maxcoeff=10000000)
    if r is not None:
        coeff_A = mp.mpf(-r[1]) / mp.mpf(r[0])
        coeff_B = mp.mpf(-r[2]) / mp.mpf(r[0])
        # Check if coefficients are integers
        A_int = round(float(coeff_A))
        B_int = round(float(coeff_B))
        err_A = abs(coeff_A - A_int)
        err_B = abs(coeff_B - B_int)
        if err_A < 1e-15 and err_B < 1e-15:
            int_table[(m, k)] = (A_int, B_int)
            print(f"  T_{m},{k} = {A_int}*A' + {B_int}*B'  ✓")
        else:
            all_ints = False
            print(f"  T_{m},{k} = {mp.nstr(coeff_A,6)}*A' + {mp.nstr(coeff_B,6)}*B'  (not integers!)")
    else:
        all_ints = False
        print(f"  T_{m},{k} = PSLQ failed")

if all_ints:
    print("\nAll T_mk are in Z*A' + Z*B'! (A' = U(3/2,1,1), B' = U(5/2,2,1))")
