"""Find the true integer normalization.

From previous analysis:
  Poch_m * Poch_k * T_mk = p/2^a * A' + q/2^b * B'
where the denominators are powers of 2. So the "right" normalization is:
  U_mk = 2^{m+k} * (3/2)_m * (3/2)_k * T_mk should be in Z*A' + Z*B'.

This is equivalent to: U_mk = 2^{m+k} * m! * k! * (-i)^{m+k} * J_mk / pref.
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
    return mp.gamma(n + mp.mpf(3)/2) / mp.gamma(mp.mpf(3)/2)

# From the pattern: (3/2)_n = (3/2)(5/2)...(2n+1)/2 = (2n+1)!! / 2^n
# So 2^n * (3/2)_n = (2n+1)!! = 1*3*5*...*(2n+1), which is an integer.
# Therefore:
#   (2m+1)!! * (2k+1)!! * T_mk = p * A' + q * B' should have integer p, q.
# Since (2n+1)!! = 2^n * (3/2)_n, this is 2^{m+k} * Poch_m * Poch_k * T_mk.

# Let's verify: (3/2)_0 = 1 = 1!!/2^0, (3/2)_1 = 3/2 = 3!!/2 = 3/2. ✓
# (2n+1)!! = 1,3,15,105,945,...
DFactorial = [1, 3, 15, 105, 945, 10395, 135135]  # (2n+1)!! for n=0,1,...

print("(2n+1)!! sequence:", DFactorial[:6])
for n in range(6):
    poch = float(Poch(n))
    pred = DFactorial[n] / 2**n
    print(f"  n={n}: 2^n*(3/2)_n = {pred:.4f}, (2n+1)!!/{2**n} = {DFactorial[n]}/{2**n} = {DFactorial[n]/2**n:.4f}")

print("\n=== Testing U_mk = (2m+1)!! * (2k+1)!! * T_mk ===")
print("(with T_mk = (-i)^{m+k} * beta_m * beta_k * J_mk / pref, beta_n = n!/Poch_n)")
print("Equivalently: U_mk = re(2^{m+k} * Poch_m * Poch_k * (-i)^{m+k} * (beta_m*beta_k) * J_mk / pref)")
print("= re(2^{m+k} * m! * k! * (-i)^{m+k} * J_mk / pref)\n")

int_table = {}
for m in range(7):
    for k in range(m, 7):
        Jmk = J(m, k)
        # U_mk = (2m+1)!! * (2k+1)!! * T_mk
        # = (2m+1)!! * (2k+1)!! * (-i)^{m+k} * (m!/Poch_m) * (k!/Poch_k) * J_mk / pref
        # = 2^{m+k} * m! * k! * (-i)^{m+k} * J_mk / pref
        U_raw = 2**(m+k) * mp.fac(m) * mp.fac(k) * (-1j)**(m+k) * Jmk / pref
        U = mp.re(U_raw)
        err = abs(mp.im(U_raw))
        int_table[(m,k)] = U
        if m <= 5 and k <= 5:
            r = mp.pslq([U, Ap, Bp], tol=1e-35, maxcoeff=100000000)
            if r is not None:
                cA = mp.mpf(-r[1]) / mp.mpf(r[0])
                cB = mp.mpf(-r[2]) / mp.mpf(r[0])
                fA = Fraction(float(cA)).limit_denominator(100)
                fB = Fraction(float(cB)).limit_denominator(100)
                print(f"U_{m},{k} = {fA}*A' + {fB}*B'  (U={mp.nstr(U,8)})")
            else:
                print(f"U_{m},{k} = PSLQ failed  (U={mp.nstr(U,8)})")

# Also try with A'' = pi^2/2 * A' and B'' = pi^2/2 * B' (absorb the pref)
# which would give integer coefficients directly in J_mk terms.
Ap2 = pref * Ap
Bp2 = pref * Bp
print("\n=== Testing with A'' = pi^2/2 * U(3/2,1,1), B'' = pi^2/2 * U(5/2,2,1) ===")
for m in range(5):
    for k in range(m, 5):
        Jmk = J(m, k)
        U_raw = 2**(m+k) * mp.fac(m) * mp.fac(k) * (-1j)**(m+k) * Jmk
        U = mp.re(U_raw)
        r = mp.pslq([U, Ap2, Bp2], tol=1e-35, maxcoeff=100000000)
        if r is not None:
            cA = mp.mpf(-r[1]) / mp.mpf(r[0])
            cB = mp.mpf(-r[2]) / mp.mpf(r[0])
            fA = Fraction(float(cA)).limit_denominator(100)
            fB = Fraction(float(cB)).limit_denominator(100)
            print(f"  U_{m},{k} = {fA}*A'' + {fB}*B''")
        else:
            print(f"  U_{m},{k} = PSLQ failed  (U={mp.nstr(U,8)})")

# Now check if U_mk / A' or U_mk / B' are rational (find the simplest expression)
print("\n=== Checking actual integer values (after dividing by A' or B') ===")
for (m,k), U in int_table.items():
    if m > 5 or k > 5:
        continue
    # Check U/A': if rational, should have small denominator
    ratio_A = float(U / float(Ap))
    ratio_B = float(U / float(Bp))
    fA = Fraction(ratio_A).limit_denominator(10000)
    fB = Fraction(ratio_B).limit_denominator(10000)
    if abs(float(fA) - ratio_A) < 1e-8:
        print(f"U_{m},{k} = {fA} * A'")
    elif abs(float(fB) - ratio_B) < 1e-8:
        print(f"U_{m},{k} = {fB} * B'")
    else:
        print(f"U_{m},{k} = {mp.nstr(U,6)} (not simple)")

# Try: maybe rescale A', B' differently.
# From the "Poch*Poch*T" table:
# S_00 = 1*A' + 0*B'     -> A' = S_00 = T_00
# S_01 = 3/2*A' - 3*B'   -> B' = (A' - (2/3)*S_01) / 2 ... messy
# Better: define C' = A' and D' = A' - 2*B' = T_01 (a "natural" linear combination)
# Then: A' = C', B' = (C' - D')/2.
# T_00 = C'
# T_01 = D' (which is A' - 2B')
# T_11 = 11/3 A' - 20/3 B' = 11/3 C' - 10/3 (C' - D') = C'/3 + 10D'/3 ... not simpler

# Or: write in terms of E = pref*A' = J_00 and F = pref*B' = (-2/3) * d/dz (pref*U(3/2,1,z))|_1
# These are related to the ACTUAL J_mk values.

# Let me instead find by looking at the pattern of the coefficients.
# U_mk coefficients in A', B':
# U_00 = A' (coeff 1, 0)
# U_01 = ? Let's compute
U01 = int_table[(0,1)]
r01 = mp.pslq([U01, Ap, Bp], tol=1e-35, maxcoeff=100000000)
print(f"\nU_01 = {mp.nstr(U01, 12)}")
print(f"PSLQ: {r01}")
print(f"A' - 2B' = {mp.nstr(Ap - 2*Bp, 12)}")
print(f"2*(A' - 2B') = {mp.nstr(2*(Ap - 2*Bp), 12)}")
# From T_01 = A' - 2B': U_01 = 2^1 * 1! * 1! * T_01 / Poch_0 / Poch_1 * Poch_0 * Poch_1
# Wait: U_01 = 2^1 * m! * k! * re((-i)^1 * J_01 / pref)
# But T_01 = beta_0 * beta_1 * re((-i)^1 * J_01 / pref) = (1/1)*(1/Poch_1) * re(...)
# = (1 * 2/3) * re((-i) * J_01 / pref)
# U_01 = 2 * 1 * 1 * re((-i) * J_01 / pref) = 2 * (1/(beta_0*beta_1)) * T_01
# = 2 * (1/1 * 1/(2/3)) * T_01 = 2 * 3/2 * T_01 = 3 * T_01 = 3*(A' - 2B').
print(f"3*(A'-2B') = {mp.nstr(3*(Ap - 2*Bp), 12)}")

# Let me also check U_11:
U11 = int_table[(1,1)]
r11 = mp.pslq([U11, Ap, Bp], tol=1e-35, maxcoeff=100000000)
print(f"\nU_11 = {mp.nstr(U11, 12)}")
print(f"PSLQ: {r11}")
# From T_11 = 11/3 A' - 20/3 B':
# U_11 = 2^2 * 1 * 1 * (1/(beta_1*beta_1)) * T_11 = 4 * Poch_1^2 T_11 = 4 * (3/2)^2 * T_11
# = 4 * 9/4 * (11/3 A' - 20/3 B') = 9 * (11/3 A' - 20/3 B') = 33A' - 60B'.
print(f"33A'-60B' = {mp.nstr(33*Ap - 60*Bp, 12)}")

# These ARE integers times A' and B'! Let me check all:
print("\n=== Checking U_mk = integer * A' + integer * B' ===")
all_int = True
int_results = {}
for (m,k), U in int_table.items():
    if m > 5 or k > 5: continue
    r = mp.pslq([U, Ap, Bp], tol=1e-35, maxcoeff=1000000000)
    if r is not None:
        cA = mp.mpf(-r[1]) / mp.mpf(r[0])
        cB = mp.mpf(-r[2]) / mp.mpf(r[0])
        iA = round(float(cA))
        iB = round(float(cB))
        if abs(cA - iA) < 1e-8 and abs(cB - iB) < 1e-8:
            int_results[(m,k)] = (iA, iB)
            print(f"  U_{m},{k} = {iA}*A' + {iB}*B'  ✓")
        else:
            all_int = False
            print(f"  U_{m},{k} = {mp.nstr(cA,6)}*A' + {mp.nstr(cB,6)}*B'  (NOT integers!)")
    else:
        all_int = False
        print(f"  U_{m},{k} = PSLQ failed  (U = {mp.nstr(U,8)})")

if all_int:
    print("\n✓ All U_mk are in Z*A' + Z*B'!")
    print("  A' = U(3/2,1,1), B' = U(5/2,2,1)")
    print("  U_mk = 2^{m+k} m! k! * re((-i)^{m+k} J_mk) / (pi^2/2)")
