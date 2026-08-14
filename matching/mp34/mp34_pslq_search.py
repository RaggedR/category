"""Find the right normalization and generators for the lambda=3/4 matching matrix.

For lambda=1/4:
  h_n = sqrt(2)*pi^{3/2} * C(2n,n)/4^n  (from the explicit formula for MP norms)
  T_mk = 2^m m! 2^k k! / h_0 * J_mk   with J_mk = c_mk * sqrt(h_m h_k)

For lambda=3/4:
  h_n = h_0 * Gamma(n+3/2)/(n! Gamma(3/2))

The key insight: for lambda=1/4, h_n = h_0 * C(2n,n)/4^n * (something that equals 1 when n=0)
Actually h_n(lambda) = Gamma(2lambda)/(2*sqrt(pi)) * Gamma(n+2lambda)/n! * (2pi)
= pi * Gamma(n+2lambda)/n! for lambda=1/4: pi*Gamma(n+1/2)/n!
= pi * (2n-1)!!/2^n * sqrt(pi)/n! ... let me just look at what cancels.

For lambda=1/4:
  h_n / h_0 = Gamma(n+1/2)/(n! Gamma(1/2)) = C(2n,n)/4^n (since Gamma(n+1/2)/Gamma(1/2) = (2n-1)!!/2^n)
  The normalization 2^m m! * 2^k k! in T_mk: absorbed the C(2n,n) factors since
  2^n n! * sqrt(C(2n,n)) ~ 2^n n! * (2n)!/(2^n n!)^2 ...

Actually, for lambda=1/4, h_0 = sqrt(2)*pi^{3/2}, and:
  sqrt(h_m * h_k) = h_0 * sqrt(C(2m,m)*C(2k,k))/4^{(m+k)/2}
So T_mk = 2^m m! 2^k k! / h_0 * c_mk * h_0 * sqrt(C(2m,m)*C(2k,k)) / 4^{(m+k)/2}
       = m! k! * sqrt(C(2m,m)*C(2k,k)) * c_mk   -- the 2^m * 2^k / 4^{(m+k)/2} = 1!

So really: T_mk = m! k! * sqrt(C(2m,m)*C(2k,k)) * c_mk

For lambda=3/4, h_n/h_0 = Gamma(n+3/2)/(n! Gamma(3/2)).
Let's define analogously:
  T_mk^(3/4) = sqrt(h_m/h_0) * sqrt(h_k/h_0) * m! k! * c_mk / h_0 ??

Actually let's try the natural normalization from the GF perspective.
The MP GF for lambda=3/4 is:
  sum_n P_n(x) w^n = (1-iw)^{-3/4+ix} (1+iw)^{-3/4-ix}
  = exp((−3/4+ix) log(1−iw) + (−3/4−ix) log(1+iw))

The coefficient of w^n in this GF involves (3/4)_n in a certain way.
The correct denominator for integrality would absorb the poles of the Gamma function.

Let me try: for lambda=3/4, the "right" normalization divides by Gamma(n+3/4) factors.
Actually let's just do a systematic PSLQ search with various normalizations.
"""
import mpmath as mp
import pickle

mp.mp.dps = 50

with open("mp34_dual_matching.pkl", "rb") as f:
    nrm, C = pickle.load(f)

h0 = nrm[0]

def c(m, k):
    return C[(min(m,k), max(m,k))]

def J(m, k):
    return c(m, k) * mp.sqrt(nrm[m] * nrm[k])

# Identify candidate constants
# J_00 = 2pi*Gamma(3/2)^2*U(3/2,1,1)
J00 = J(0, 0)
U32 = mp.hyperu(mp.mpf(3)/2, 1, 1)
U52 = mp.hyperu(mp.mpf(5)/2, 2, 1)
U72 = mp.hyperu(mp.mpf(7)/2, 3, 1)

# U(3/2,1,1) in terms of exponential integral / Bessel?
# U(a,b,z) = 1/Gamma(a) * int_0^inf e^{-zt} t^{a-1} (1+t)^{b-a-1} dt
# U(3/2, 1, z) = 1/Gamma(3/2) * int_0^inf e^{-zt} t^{1/2} (1+t)^{-3/2} dt
# At z=1: = 2/sqrt(pi) * int_0^inf e^{-t} sqrt(t)/(1+t)^{3/2} dt
# Try: U(3/2, 1, 1) = 2/sqrt(pi) * (...)
# Actually U(1/2, 1, 2z) = sqrt(pi/(2z)) K_{1/4}(z)... no that's not right.
# Let's look at it differently.
# From DLMF 13.6.21: U(1/2, 1, z) = sqrt(pi) e^{z/2} K_0(z/2) / sqrt(pi) ...
# = e^{z/2} K_0(z/2). Wait that's what we know.
# For U(3/2, 1, z): use 13.4.24: U(a,1,z) = -1/Gamma(a)(ln(z)+psi(a)) M(a,1,z) + ...
# Or use recurrence: (b-a) U(a,b,z) = U(a-1,b,z) - z*U(a,b+1,z) -- complex
# Simpler: U(a+1,b+1,z) = -dU(a,b,z)/dz / a (diff wrt z)
# U(3/2, 1, z) = -d/dz U(1/2, 0, z) ... but b=0 is degenerate.
# Actually d/dz U(a,b,z) = -a U(a+1,b+1,z)
# So U(3/2, 2, z) = -1/(1/2) * d/dz U(1/2, 1, z) ... no wait:
# d/dz U(a,b,z)|_{a=1/2,b=1} = -(1/2) U(3/2, 2, z)
# and U(1/2,1,z) = e^{z/2} K_0(z/2), so d/dz = 1/2 e^{z/2} K_0(z/2) + e^{z/2}*1/2*(-K_1(z/2))
# = e^{z/2}/2 * (K_0(z/2) - K_1(z/2))
# So U(3/2, 2, z) = -2 * e^{z/2}/2 * (K_0(z/2) - K_1(z/2)) ...
# Hmm but we want U(3/2, 1, z), not U(3/2, 2, z).
#
# Use contiguous relation: U(a,b,z) - U(a,b-1,z) = (b-2) U(a,b-1,z) - ...
# Actually: (b-1) U(a,b,z) = U(a-1,b-1,z) + (z-b+1) U(a,b-1,z) -- Kummer's relation
# Simpler: use the integral definition directly.
# U(3/2, 1, 1): just compute numerically and find its expression.
print(f"J_00 = {mp.nstr(J00, 20)}")
print(f"h_0 = {mp.nstr(h0, 20)}")
print(f"U(3/2,1,1) = {mp.nstr(U32, 20)}")
print(f"U(5/2,2,1) = {mp.nstr(U52, 20)}")
print(f"U(7/2,3,1) = {mp.nstr(U72, 20)}")

# Key: for lambda=3/4, try to express U(3/2,1,1) and U(5/2,2,1) differently
# From the identity in the lambda=1/4 paper: U(1/2,1,2z) = e^z K_0(z)/sqrt(pi)
# For lambda=3/4, the GF will involve U(3/2,1,zeta) where zeta=phi(w)phi(u).
# Let's find the right expression for U(3/2,1,z).
#
# From the contiguous relation (DLMF 13.3.3):
# U(a,b,z) = U(a+1,b+1,z)*z + U(a,b+1,z)*(1+a-b)
# At a=1/2, b=1: U(1/2,1,z) = z U(3/2,2,z) + (1/2-1) U(1/2,2,z) + ...
# This is getting complicated. Let me try a different approach:
#
# For U(3/2, 1, z): from DLMF 13.6 "Exponential and Logarithmic Integrals":
# U(1,1,z) = e^z E_1(z) where E_1 is exponential integral
# There's also: U(n+1/2, 1, z) type formulas via parabolic cylinder functions
#
# Let me just look for A', B' numerically.
# From the structure: J_00 / h_0 = c_00 = 0.897..., J_01/h_0 = ...
# For lambda=1/4 the formula was: J_00/h_0 = A (one generator)
# But J_00/h_0 for l=1/4 = c_00 = (value from U(1/2,1,1) = e^{1/2}K_0(1/2)/sqrt(pi))
# = 2pi * pi * U(1/2,1,1) / h_0 = ...

# Let's find A', B' by looking at what T_00 and T_01 are.
# For lambda=1/4: A = sqrt(2e) K_0(1/2), B = sqrt(2e) K_1(1/2)
# These are values of sqrt(2e) * K_nu(1/2) for nu=0,1.
# For lambda=3/4: try K_nu(1/2) for various nu, or K_nu(1) for z=1...
# Or maybe we need K_0(1), K_1(1) since U(3/2,1,1) involves argument 1 not 1/2.

# From the GF for lambda=1/4: zeta = phi(w)phi(u), at w=u=0: zeta=1.
# The GF was proportional to e^{zeta/2} K_0(zeta/2) = e^{1/2} K_0(1/2).
# For lambda=3/4: the GF will be proportional to some function of zeta.
# The 3-Gamma lemma gives U(3/2, 1, zeta).
# If U(3/2, 1, z) = f(z) for some Bessel expression...

# Let me try: U(3/2, 1, z) = C * e^{z/2} * (a0*K_0(z/2) + a1*K_1(z/2)) for some C,a0,a1
# At z=1: U(3/2,1,1) = C * e^{1/2} * (a0*K_0(1/2) + a1*K_1(1/2))
z_test = 1
U_z = mp.hyperu(mp.mpf(3)/2, 1, z_test)
K0_h = mp.besselk(0, mp.mpf(1)/2)
K1_h = mp.besselk(1, mp.mpf(1)/2)
esqrt = mp.exp(mp.mpf(1)/2)
# Try PSLQ: U(3/2,1,1) = a * e^{1/2} K_0(1/2) + b * e^{1/2} K_1(1/2)
r = mp.pslq([U_z, esqrt * K0_h, esqrt * K1_h], tol=1e-25, maxcoeff=1000)
print(f"\nPSLQ U(3/2,1,1) = a*e^{{1/2}}K_0(1/2) + b*e^{{1/2}}K_1(1/2): {r}")

# Try other Bessel args
K0_1 = mp.besselk(0, 1)
K1_1 = mp.besselk(1, 1)
e1 = mp.exp(1)
r2 = mp.pslq([U_z, e1*K0_1, e1*K1_1], tol=1e-25, maxcoeff=1000)
print(f"PSLQ U(3/2,1,1) = a*e*K_0(1) + b*e*K_1(1): {r2}")

# What about U(3/2,1,1) in terms of modified Bessel I?
# Or parabolic cylinder?
# Let's try to check if U(3/2,1,z) = sqrt(pi) * something via series
# U(a,b,z) ~ z^{-a} for large z
# For a=3/2, U(3/2,1,z) ~ z^{-3/2} ... decays faster
# Derivative: d/dz U(3/2,1,z) = -3/2 * U(5/2,2,z)
# Second deriv: d^2/dz^2 U(3/2,1,z) = -3/2 * d/dz U(5/2,2,z) = 15/4 * U(7/2,3,z)

# Check if J_mk/J_00 can be expressed in U's with rational coefficients
print("\n--- Ratios J_mk / J_00 ---")
for m in range(4):
    for k in range(m, 4):
        jmk = J(m, k)
        rat = jmk / J00
        print(f"  J_{m},{k}/J_00 = {mp.nstr(mp.re(rat / (-1j)**abs(m-k)), 12)}")

# The key: find what normalization gives integers.
# Try working backwards from the GF.
# For lambda=1/4, the GF of J_mk was 2pi^{3/2} e^{zeta/2} K_0(zeta/2) / sqrt((1-iw)(1-iu))
# The generating function of P_n(x) w^n for lambda=1/4 was (1-iw)^{-1/4+ix}(1+iw)^{-1/4-ix}
# For lambda=3/4: (1-iw)^{-3/4+ix}(1+iw)^{-3/4-ix}
#
# After substituting the GF in both indices:
# sum_{m,k} J_mk w^m u^k = int Gamma(3/4+ix)^2 Gamma(3/4-ix) *
#    (1-iw)^{-3/4+ix}(1+iw)^{-3/4-ix} * (1-iu)^{-3/4+ix}(1+iu)^{-3/4-ix} dx
# = int Gamma(3/4+ix)^2 Gamma(3/4-ix) * [(1-iw)(1-iu)]^{-3/4+ix} [(1+iw)(1+iu)]^{-3/4-ix} dx
# = int Gamma(3/4+ix)^2 Gamma(3/4-ix) * [(1-iw)(1-iu)/(1+iw)(1+iu)]^{ix} *
#    [(1-iw)(1-iu)(1+iw)(1+iu)]^{-3/4} dx
#
# Let alpha = log((1-iw)(1-iu)/(1+iw)(1+iu)) = is with e^{-s} = zeta = phi(w)phi(u)
# So e^{ix*alpha} = e^{-sx} where s = -i*alpha = log(1/zeta)
# And [(1-iw)(1+iw)(1-iu)(1+iu)]^{-3/4} = [(1+w^2)(1+u^2)]^{-3/4}
#
# sum J_mk w^m u^k = [(1+w^2)(1+u^2)]^{-3/4} * int Gamma(3/4+ix)^2 Gamma(3/4-ix) e^{isx} dx
# where e^{-s} = zeta = phi(w)phi(u)
# By shifted 3-Gamma lemma (a=b=c=3/4):
# = [(1+w^2)(1+u^2)]^{-3/4} * 2pi * e^{-3s/4} * Gamma(3/2)^2 * U(3/2, 1, e^{-s})
# = [(1+w^2)(1+u^2)]^{-3/4} * 2pi * Gamma(3/2)^2 * zeta^{3/4} * U(3/2, 1, zeta)
#
# And (1+w^2)(1+u^2) = (1-iw)(1+iw)(1-iu)(1+iu)
# = [(1-iw)(1-iu)] * [(1+iw)(1+iu)]
# with zeta = (1+iw)(1+iu)/[(1-iw)(1-iu)]
# So (1-iw)(1-iu) = sqrt((1-iw)^2(1-iu)^2) = sqrt([(1+w^2)(1+u^2)/zeta])
# Hmm, let me be more careful.
# phi(w) = (1+iw)/(1-iw), so zeta = phi(w)*phi(u) = (1+iw)(1+iu)/((1-iw)(1-iu))
# (1+w^2) = (1-iw)(1+iw), (1+u^2) = (1-iu)(1+iu)
# (1+w^2)(1+u^2) = (1-iw)(1+iw)(1-iu)(1+iu)
# Let A = (1-iw)(1-iu), B = (1+iw)(1+iu), then zeta = B/A, AB = (1+w^2)(1+u^2)
# (1+w^2)(1+u^2) = AB = A * zeta * A = A^2 * zeta... no wait B = zeta*A so AB = zeta*A^2
# So [(1+w^2)(1+u^2)]^{-3/4} = (zeta*A^2)^{-3/4} = zeta^{-3/4} * A^{-3/2}
# and A = (1-iw)(1-iu)
#
# sum J_mk w^m u^k = zeta^{-3/4} * [(1-iw)(1-iu)]^{-3/2} * 2pi * Gamma(3/2)^2 * zeta^{3/4} * U(3/2, 1, zeta)
# = 2pi * Gamma(3/2)^2 * U(3/2, 1, zeta) / [(1-iw)(1-iu)]^{3/2}
# = pi^2/2 * U(3/2, 1, zeta) / [(1-iw)(1-iu)]^{3/2}
# (using Gamma(3/2) = sqrt(pi)/2, so Gamma(3/2)^2 = pi/4)

# Let's verify this GF prediction!
def GF_34(w, u):
    zeta = (1 + 1j*w) * (1 + 1j*u) / ((1 - 1j*w) * (1 - 1j*u))
    return (2 * mp.pi * mp.gamma(mp.mpf(3)/2)**2 * mp.hyperu(mp.mpf(3)/2, 1, zeta)
            / ((1 - 1j*w) * (1 - 1j*u))**mp.mpf(3)/2)

print("\n--- GF verification (compare J_mk from GF vs quadrature) ---")
# Extract Taylor coefficients by evaluating at small w, u
# J_mk = [w^m u^k] GF_34(w,u) * m! * k! (from Taylor expansion)
# Actually it's a formal power series, use numerical differentiation
# Better: use Fourier sum like in mp_dual_row_coeffs.py
Mw, Mu = 32, 32
r = mp.mpf(1) / 4
ws = [r * mp.expjpi(2 * mp.mpf(a) / Mw) for a in range(Mw)]
us = [r * mp.expjpi(2 * mp.mpf(b) / Mu) for b in range(Mu)]
vals = [[GF_34(ws[a], us[b]) for b in range(Mu)] for a in range(Mw)]

print("Comparison of J_mk (GF) vs J_mk (quadrature):")
for m in range(5):
    for k in range(m, 5):
        # Fourier extraction
        inner = [sum(vals[a][b] * mp.expjpi(-2 * mp.mpf(m*a) / Mw) for a in range(Mw)) / (Mw * r**m)
                 for b in range(Mu)]
        J_gf = sum(inner[b] * mp.expjpi(-2 * mp.mpf(k*b) / Mu) for b in range(Mu)) / (Mu * r**k)
        J_quad = J(m, k)
        diff = abs(J_gf - J_quad)
        print(f"  J_{m},{k}: GF={mp.nstr(mp.re(J_gf/(-1j)**abs(m-k)), 8)}, "
              f"quad={mp.nstr(mp.re(J_quad/(-1j)**abs(m-k)), 8)}, "
              f"diff={mp.nstr(diff, 3)}")

# Now identify the normalization that gives integers.
# GF = pi^2/2 * U(3/2, 1, zeta) / [(1-iw)(1-iu)]^{3/2}
# The MP GF for lambda=3/4: sum P_n w^n = (1-iw)^{-3/4+ix}(1+iw)^{-3/4-ix}
# Taylor coefficients involve (3/4)_n type rising factorials.
# For the T_mk normalization: look at what the Taylor coefficients of
# 1/(1-iw)^{3/2} are: [w^n] (1-iw)^{-3/2} = i^n * C(n+1/2, n) = i^n * Gamma(n+3/2)/(n! Gamma(3/2))
# So the denominator factor for each index is Gamma(n+3/2)/(n! Gamma(3/2)) times i^n.
# The right normalization should be:
# T_mk = J_mk * n_m * n_k / (J_00) where n_n = n! Gamma(3/2) / Gamma(n+3/2)
# = J_mk / (J_00 * h_m/h_0 * h_k/h_0)  -- if we normalize by h ratios?
# No. Let me think again.
#
# For lambda=1/4: the GF was e^{zeta/2}K_0(zeta/2)/sqrt((1-iw)(1-iu))
# The 1/sqrt((1-iw)(1-iu)) factor contributes [w^m] 1/sqrt(1-iw) = i^m * C(m-1/2,m) = i^m * C(2m,m)/4^m
# And 2^m m! * i^m * C(2m,m)/4^m = i^m * (2m)!/4^m / m! ... hmm 2^m m! * C(2m,m)/4^m = (2m)!/(2^m m!)
# For lambda=3/4: [w^m] 1/(1-iw)^{3/2} = i^m * Gamma(m+3/2)/(m! Gamma(3/2)) = i^m * C(m+1/2,m)
# The right normalization:
# N_m = m! / (i^m * coefficient of w^m in 1/(1-iw)^{3/2}) = m! Gamma(3/2) / Gamma(m+3/2)
# = Gamma(3/2) / (m-1+3/2)(m-2+3/2)...(3/2) = 1/((m+1/2)(m-1+1/2)...(3/2)) ...
# Actually N_m = m! Gamma(3/2) / Gamma(m+3/2) = 1 / (3/2)_m * m! where (3/2)_m = Gamma(m+3/2)/Gamma(3/2)
# N_m * N_k in the denominator would cancel out the GF prefactors.

# The correct T_mk for lambda=3/4:
# T_mk = J_mk / (J_00 * [coeff of w^m in 1/(1-iw)^{3/2}] * [coeff of u^k in 1/(1-iu)^{3/2}])
# = J_mk * [m! Gamma(3/2)/Gamma(m+3/2)] * [k! Gamma(3/2)/Gamma(k+3/2)] / J_00

print("\n--- Trying normalization T_mk = J_mk * [m! Gamma(3/2)/Gamma(m+3/2)] * [k! Gamma(3/2)/Gamma(k+3/2)] / J_00 ---")
G32 = mp.gamma(mp.mpf(3)/2)

# U(3/2,1,1) is our A' = J_00/J_00 = 1? No wait.
# The GF: sum J_mk w^m u^k = J_00 * U(3/2,1,zeta)/U(3/2,1,1) / [(1-iw)(1-iu)]^{3/2}
# So J_mk = J_00 * [w^m u^k] U(3/2,1,zeta) / [(1-iw)(1-iu)]^{3/2}
# With N_m = m! Gamma(3/2)/Gamma(m+3/2):
# N_m * N_k * J_mk / J_00 = N_m * N_k * [w^m u^k] U(3/2,1,zeta) / [(1-iw)(1-iu)]^{3/2}
# [w^m u^k] 1/(1-iw)^{3/2} = i^m * Gamma(m+3/2)/(m! Gamma(3/2))
# N_m * Gamma(m+3/2)/(m! Gamma(3/2)) = 1
# So N_m * N_k * J_mk / J_00 = [w^m u^k] U(3/2,1,zeta) (as a formal series in integers?)

# U(3/2,1,z) Taylor series: U(3/2,1,z) = sum_n u_n z^n
# Actually U(a,b,z) is NOT a power series at 0 in general (logarithmic singularity).
# U(a,1,z) ~ -1/Gamma(a) * (ln z + psi(a)) + ... as z->0+ (DLMF 13.2.11)
# So zeta -> 0 as w -> 0, u -> 0, and U has a log singularity there!
# This means the GF has a log singularity at (w,u)=(0,0), so J_mk grow...
# Actually: at w=u=0, zeta = phi(0)phi(0) = 1. So U(3/2,1,1) is finite. Good.
# As w -> something with |phi(w)| = 1, zeta has modulus 1 and the series converges...
# The issue is where the GF has singularities in w: at w = -i/1 = -i, where (1-iw) = 0.
# At w = i: zeta = (1+i*i)/(1-i*i) * phi(u) = 0 * phi(u) = 0. And U(3/2,1,0) -> infinity!
# So the singularity structure is: (1-iw) = 0 AND log singularity of U as zeta->0.

# This is analogous to lambda=1/4: K_0(zeta/2) has a log singularity as zeta->0.
# And the K_0 GF also had (1-iw)^{-1/2} prefactor.
# So the structure IS analogous, but the generators are different.

# For lambda=1/4:
# T_mk = J_mk * [m!*C(2m,m)^{1/2}/2^m] * [k!*C(2k,k)^{1/2}/2^k] * (2/h_0)
# ... actually let me re-derive from the GF.
# GF_lambda14 = 2*pi^{3/2} * e^{zeta/2} K_0(zeta/2) / sqrt((1-iw)(1-iu))
# = 2*pi^{3/2} * [function of zeta] / [(1-iw)(1-iu)]^{1/2}
#
# [w^m] 1/(1-iw)^{1/2} = i^m * C(2m,m)/4^m
# N_m^{1/4} = m! * 4^m / C(2m,m) = 4^m * (m!)^2 / (2m)! (removes the i^m * C(2m,m)/4^m factor)
# T_mk^{1/4} = N_m^{1/4} * N_k^{1/4} * J_mk / J_00
# = [w^m u^k] e^{zeta/2} K_0(zeta/2) / (e^{1/2} K_0(1/2))   -- should be rational!
# Hmm, but T_mk/h_0 was in ZA + ZB, not rational.

# Actually let's reread: T_mk = 2^m m! * 2^k k! / h_0 * J_mk
# with J_mk = c_mk * sqrt(h_m h_k).
# For lambda=1/4: h_n = sqrt(2)*pi^{3/2} * C(2n,n)/4^n
# sqrt(h_m * h_k) = sqrt(2)*pi^{3/2} * C(2m,m)^{1/2} * C(2k,k)^{1/2} / (2^{m+k})
# T_mk = 2^m m! * 2^k k! / h_0 * c_mk * sqrt(2)*pi^{3/2} * C(2m,m)^{1/2} * C(2k,k)^{1/2} / (2^{m+k})
#       = m! k! * C(2m,m)^{1/2} * C(2k,k)^{1/2} * c_mk  (since h_0 = sqrt(2)*pi^{3/2})

# So T_mk^{1/4} = m! k! sqrt(C(2m,m) C(2k,k)) c_mk.
# And this is in ZA + ZB.

# For lambda=3/4:
# h_n = h_0 * Gamma(n+3/2)/(n! Gamma(3/2))
# sqrt(h_m * h_k) = h_0 * sqrt(Gamma(m+3/2) Gamma(k+3/2) / (m! k! Gamma(3/2)^2))
# Try T_mk^{3/4} = m! k! / Gamma(3/2)^{-1} * sqrt(...) * c_mk ?
# Let alpha_n = sqrt(Gamma(n+3/2)/(n! Gamma(3/2)))  (so h_n = h_0 alpha_n^2)
# sqrt(h_m h_k) = h_0 alpha_m alpha_k
# c_mk = J_mk / (h_0 alpha_m alpha_k)
# T_mk^{3/4} = (m! Gamma(3/2)/Gamma(m+3/2)) * (k! Gamma(3/2)/Gamma(k+3/2)) * J_mk / J_00
#            = (m! Gamma(3/2)/Gamma(m+3/2)) * (k! Gamma(3/2)/Gamma(k+3/2)) * c_mk * h_0 alpha_m alpha_k / J_00
#            = (m!/Gamma(m+3/2))^{1/2} * (k!/Gamma(k+3/2))^{1/2} * c_mk * 1 / U(3/2,1,1)
# Hmm. Let me define:
# beta_n = Gamma(3/2)/Gamma(n+3/2) * n! = 1/(3/2)_n  (reciprocal of Pochhammer)
# T_mk^{3/4} = sqrt(beta_m * beta_k) * c_mk / (pi^2/2 * U(3/2,1,1) / h_0 * h_0)
# No. Let me just try different normalizations numerically and see what PSLQ finds.

# The GF is: sum J_mk w^m u^k = pi^2/2 * U(3/2,1,zeta) / [(1-iw)(1-iu)]^{3/2}
# The coefficient of w^m u^k in 1/(1-iw)^{3/2} = i^m * (3/2)_m / m! where (3/2)_m = Gamma(m+3/2)/Gamma(3/2)
# So if we define X_m = m! Gamma(3/2)/Gamma(m+3/2), then:
# i^{m+k} * J_mk * X_m * X_k / (pi^2/2) = [w^m u^k] U(3/2,1,zeta)   (up to sign from i^{m+k})
# But U(3/2,1,zeta) has a log singularity at zeta=0, so it's not a Taylor series!
# Wait: zeta = phi(w)phi(u), and at w=u=0, zeta=1. So we're expanding in powers of w,u
# around 0, which sends zeta -> phi(w)phi(u) that starts at 1 when w=u=0.
# As w,u -> 0: phi(w) ~ (1+iw)/(1-iw) ~ 1 + 2iw, so zeta ~ 1 + 2i(w+u).
# So U(3/2,1,zeta) expanded around zeta=1 is an analytic Taylor series! Good.

# Let delta = zeta - 1 = phi(w)phi(u) - 1 (which -> 0 as w,u -> 0)
# Taylor coefficients of U(3/2,1,1+delta) in delta are finite.
# And [w^m u^k] delta^j involves products of Taylor coefficients of phi(w)-1.
# phi(w)-1 = 2iw/(1-iw) ~ 2iw + 2w^2 + ...

# So the correct normalization that makes T_mk integers (or in ZA+ZB) involves
# the Taylor coefficients of U(3/2,1,1+delta) in delta times combinatorial factors.

# For A', B': U(3/2,1,zeta) and its derivative at zeta=1.
# d/dz U(3/2,1,z) = -(3/2) U(5/2,2,z)
# So A' ~ U(3/2,1,1) and B' ~ U(5/2,2,1).

# Can we relate these to Bessel functions?
# d/dz [z^{3/2} U(3/2,1,z)] ... let me try a substitution.
# Actually from DLMF: for U(a,1,z), there is a connection to the Meijer G function.
# Let's try: is there a clean expression?

# Numerical approach: try PSLQ on T_00, T_01 with many constant candidates
Ap = mp.hyperu(mp.mpf(3)/2, 1, 1)   # U(3/2,1,1)
Bp = mp.hyperu(mp.mpf(5)/2, 2, 1)   # U(5/2,2,1)
# From the contiguous relation: U(3/2,1,z) = z*U(5/2,2,z) + (1/2) U(3/2,2,z)...
# and U(3/2,2,z) can be expressed... complex. Let's just use Ap, Bp.

print(f"\nA' = U(3/2,1,1) = {mp.nstr(Ap, 20)}")
print(f"B' = U(5/2,2,1) = {mp.nstr(Bp, 20)}")

# Define normalization factors
def X(n):
    return mp.fac(n) * mp.gamma(mp.mpf(3)/2) / mp.gamma(n + mp.mpf(3)/2)

# T_mk^{normalized} = J_mk * X(m) * X(k) / (pi^2/2)
pref = mp.pi**2 / 2
print("\nT_mk = J_mk * X(m) * X(k) / (pi^2/2) -- should be in Z*Ap + Z*Bp:")
for m in range(5):
    for k in range(m, 5):
        T_raw = J(m, k) * X(m) * X(k) / pref
        # Remove phase: (-i)^{|m-k|}
        T = mp.re(T_raw / (-1j)**abs(m-k))
        r = mp.pslq([T, Ap, Bp], tol=1e-20, maxcoeff=10000)
        print(f"  T_{m},{k} = {mp.nstr(T, 10)}, PSLQ: {r}")
