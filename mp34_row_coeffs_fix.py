"""Fix the row coefficient extraction from the GF.

The GF has singularities at w = i (where 1+iw = 0), so the Taylor series
in w has radius of convergence 1. We used r=1/4 but the aliasing
might be due to the U(3/2,1,q) having additional singularities.

When w -> i, q = (1-iw)(1-iu)/((1+iw)(1+iu)) -> 0/(0*(something)) which is 0.
U(3/2,1,0) is infinite (log singularity). So the GF has a singularity at w=i.
For small r, the Taylor extraction should work.

The issue: U(3/2,1,q) for complex q. As w goes on a circle of radius r<1,
q = (1-iw)(1-iu)/((1+iw)(1+iu)) traces an arc. We need q to avoid the
negative real axis (where U has a branch cut) and 0 (where U has log singularity).

For w on circle |w|=r<1, u=0: q = (1-iw)/(1+iw). This is the Cayley transform
of w, mapping |w|<1 to Im(q)>0 (upper half plane)... wait actually:
  phi(w) = (1+iw)/(1-iw), so 1/phi(w) = (1-iw)/(1+iw) = q|_{u=0}.
For |w|<1: |phi(w)| can be > or < 1 depending on w.
For w = r e^{i theta}: |phi(w)|^2 = |1+iw|^2/|1-iw|^2.
|1+iw|^2 = (1-r sin theta)^2 + r^2 cos^2 theta = 1 - 2r sin theta + r^2.
|1-iw|^2 = 1 + 2r sin theta + r^2.
So |q|^2 = (1 + 2r sin theta + r^2) / (1 - 2r sin theta + r^2).
For theta = pi/2 (w = ir), |q|^2 = (1+2r+r^2)/(1-2r+r^2) = ((1+r)/(1-r))^2 > 1.
For theta = -pi/2 (w = -ir), |q|^2 = (1-2r+r^2)/(1+2r+r^2) < 1.
So q CAN have magnitude > 1, which is fine for U(3/2,1,q) (entire in q for Re(q)>0).
But if q has Re(q) < 0, we hit the branch cut.

Re(q) at w = r e^{i*pi} = -r (real negative): q = (1-i(-r))/(1+i(-r)) = (1+ir)/(1-ir) = phi(r).
|phi(r)| = 1 (unimodular for real r). This is on the unit circle.
For w = r (real positive): q = (1-ir)/(1+ir) = phi(-r)^{-1} = conj(phi(r)).

The key: U(3/2,1,q) is analytic for Re(q) > 0. For q on the unit circle it's still analytic.
The issue with the extraction is likely just numerical precision due to large cancellations.

Better approach: use the symbolic formula to extract Taylor coefficients exactly (as rational
multiples of A' and B'), then c_mk = T_mk / (2^{m+k} m! k! h34(m)^{1/2} h34(k)^{1/2}).

Actually wait -- the c_0k values from the GF row extension were GROWING (not decaying).
This is a sign that the GF extraction has the WRONG normalization.
Let me check: C_rows[(m,k)] = J_gf / sqrt(h34(m) * h34(k)), but the J_gf from the
Fourier sum already includes the pi^2/2 prefactor.

Actually let me trace through:
vals = GF_34(w, u) = pref * U(3/2,1,q) / ((1+iw)(1+iu))^{3/2}
GF_34(w,u) = sum_{m,k} J_mk w^m u^k.
Fourier extraction: J_mk = [w^m u^k] GF_34.
c_mk = J_mk / sqrt(h_m h_k).

So C_rows[(m,k)] = [w^m u^k] GF_34 / sqrt(h34(m) h34(k)) should equal c(m,k).
But h34(m) might be wrong! Let me check h34.

h34(n) = h0 * Poch_n / n!. From our computation: h_0 = 1.9687... and h_1 = 2.9531...
h34(1) = h0 * (3/2) / 1 = 1.9687 * 1.5 = 2.9530. ✓

So the formula should be right. The issue must be in the Fourier extraction.

Let me use larger M (more Fourier points) and smaller r.
"""
import mpmath as mp
import pickle

mp.mp.dps = 60

with open("mp34_dual_matching.pkl", "rb") as f:
    nrm, C = pickle.load(f)

h0 = nrm[0]
pref = mp.pi**2 / 2

def GF_34(w, u):
    q = (1 - 1j*w) * (1 - 1j*u) / ((1 + 1j*w) * (1 + 1j*u))
    return pref * mp.hyperu(mp.mpf(3)/2, 1, q) / ((1 + 1j*w) * (1 + 1j*u))**mp.mpf(3)/2

def h34(n):
    return h0 * mp.gamma(n + mp.mpf(3)/2) / (mp.gamma(mp.mpf(3)/2) * mp.fac(n))

def c(m, k):
    return C[(min(m,k), max(m,k))]

# Test: manually compute GF at specific (w,u) and compare to sum of known J_mk
# sum_{m,k<=11} J_mk w^m u^k for small w, u
w0, u0 = mp.mpf(1)/10, mp.mpf(1)/10
gf_test = GF_34(w0, u0)
sum_test = sum(c(m,k) * mp.sqrt(h34(m)*h34(k)) * w0**m * u0**k
               for m in range(12) for k in range(12) if m <= k) + \
           sum(c(m,k) * mp.sqrt(h34(m)*h34(k)) * w0**m * u0**k
               for m in range(12) for k in range(12) if k < m)
print(f"GF(0.1,0.1) = {mp.nstr(gf_test, 12)}")
print(f"sum J_mk w^m u^k = {mp.nstr(sum_test, 12)}")
print(f"Diff = {mp.nstr(abs(gf_test - sum_test), 5)}")

# The GF should be close to the sum if the series converges well at w=u=0.1.
# If diff is large, the GF formula is wrong.

# Also test at w=u=0: GF should equal J_00
gf_00 = GF_34(0, 0)
J00 = c(0,0) * h34(0)
print(f"\nGF(0,0) = {mp.nstr(gf_00, 12)}")
print(f"J_00 = {mp.nstr(J00, 12)}")

# The GF formula derivation was confirmed analytically and numerically (J_00 matches).
# Let me trace through the Fourier extraction more carefully.
# For r=1/4, M=80:
# [w^m u^k] GF = (1/M^2) sum_{a,b} GF(r*e^{2pi i a/M}, r*e^{2pi i b/M}) * e^{-2pi i (ma+kb)/M} / r^{m+k}
# This is exact (Cauchy formula) if |GF| series has terms r^{m+k} decaying (convergent series).
# For large m+k, J_mk/r^{m+k} grows if |J_mk| grows faster than r^{m+k}.
# From the integer table: U_mk ~ large integers. U_mk = 2^{m+k} m! k! J_mk / pref.
# So |J_mk| ~ |U_mk| * pref / (2^{m+k} m! k!).
# If |U_mk| grows fast (exponentially), then J_mk may not decay.
# From the table: U_05 = 165579 A' - 276270 B' ~ 165579 * 0.358 - 276270 * 0.215 ~ 59277 - 59399 ≈ -122.
# So U_05 is the difference of large numbers! The decay rate comes from cancellation.

# The issue with Fourier extraction: we need r such that J_mk * r^{-(m+k)} decays.
# For the GF to converge, we need |w|, |u| < 1 (singularity at w=i, u=i).
# The Taylor coefficients J_mk should decay like r_0^{m+k} for r_0 = (radius of convergence).
# But if there's strong cancellation in U_mk, J_mk itself may decay.
# From c_0k: |c_0k| = 0.897, 0.223, 0.238, 0.131, 0.123, 0.088, 0.080, 0.064, 0.058, 0.050, 0.045, 0.040.
# These are DECAYING! Good. So c(0,k) decays.
# J_0k = c(0,k) * sqrt(h0 * h_k). h_k ~ k^{1/2} (grows slowly).
# So J_0k ~ c(0,k) * sqrt(h_k) decays (times sqrt of slowly growing).
# The Fourier extraction should work with r=1/4.

# Let me check the Fourier extraction more carefully.
# The correct formula: J_mk = coeff of w^m u^k in GF(w,u) (as a power series in w,u).
# Fourier extraction via DFT: J_mk = (1/M^2) sum GF(r e^{2pi ia/M}, r e^{2pi ib/M}) e^{-2pi i(am+bk)/M} / r^{m+k}.
# Aliasing error: sum_{j: j != 0} J_{m+jM, k} r^{jM} + ...
# For r=1/4 and M=80: aliasing error ~ |J_{m+80, k}| * (1/4)^{80} = negligible.
# But for the double DFT with DIFFERENT m,k extraction:
# When we fix m and do inner sum over a, then extract k by DFT over b,
# the result is a vector of [GF(r^{m+jM}e^{...}, ...)] extracted coefficients.
# This should work.

# Let me redo the extraction with print statements to debug.
Mw, Mu = 64, 64
r = mp.mpf(1) / 4
ws = [r * mp.expjpi(2 * mp.mpf(a) / Mw) for a in range(Mw)]
us = [r * mp.expjpi(2 * mp.mpf(b) / Mu) for b in range(Mu)]

print("\nChecking GF values on the circle (w real, u=0):")
for i in [0, Mw//4, Mw//2]:
    w_test = ws[i]
    u_test = mp.mpf(0)
    gf = GF_34(w_test, u_test)
    print(f"  GF(r*e^{{2pi*{i}/{Mw}}}, 0) = {mp.nstr(gf, 8)}")

# The GF is analytic and bounded on this circle -- should be fine.
# Let me extract J_00 and J_01 manually:
vals_row0 = [GF_34(mp.mpf(0), us[b]) for b in range(Mu)]
# [u^0] GF(0,u) = J_00
J00_extracted = sum(vals_row0[b] for b in range(Mu)) / Mu / r**0
print(f"\nJ_00 extracted: {mp.nstr(J00_extracted, 12)}")
print(f"J_00 actual: {mp.nstr(J00, 12)}")

# [u^1] GF(0,u):
J01_extracted = sum(vals_row0[b] * mp.expjpi(-2*mp.mpf(b)/Mu) for b in range(Mu)) / (Mu * r**1)
J01_actual = c(0,1) * mp.sqrt(h34(0) * h34(1))
print(f"J_01 extracted: {mp.nstr(J01_extracted, 12)}")
print(f"J_01 actual: {mp.nstr(J01_actual, 12)}")

# These should match to high precision (aliasing error ~ (1/4)^{64} = 10^{-38}).
