"""Systematic analysis of lambda=3/4 matching matrix.

GF derivation:
  J_mk = int P_m(x) P_k(x) Gamma(3/4+ix)^2 Gamma(3/4-ix) dx
  Using GF: sum_n P_n(x) w^n = (1-iw)^{-3/4+ix} (1+iw)^{-3/4-ix}

  sum_{m,k} J_mk w^m u^k = int Gamma(3/4+ix)^2 Gamma(3/4-ix) *
      (1-iw)^{-3/4+ix}(1+iw)^{-3/4-ix} * (1-iu)^{-3/4+ix}(1+iu)^{-3/4-ix} dx

  = int Gamma(3/4+ix)^2 Gamma(3/4-ix) *
      [(1-iw)(1-iu)]^{-3/4+ix} [(1+iw)(1+iu)]^{-3/4-ix} dx

  Let e^{is} = (1+iw)(1+iu)/[(1-iw)(1-iu)] = phi(w)*phi(u) = zeta (Cayley product)
  So e^{isx} = [(1+iw)(1+iu)/(1-iw)(1-iu)]^{ix} = zeta^{ix}.
  And [(1-iw)(1-iu)]^{-3/4+ix}[(1+iw)(1+iu)]^{-3/4-ix}
    = [(1-iw)(1-iu)(1+iw)(1+iu)]^{-3/4} * [(1+iw)(1+iu)/(1-iw)(1-iu)]^{ix}
    = [(1+w^2)(1+u^2)]^{-3/4} * e^{isx}

  where e^{-s} = 1/zeta = phi(w)^{-1} phi(u)^{-1}.
  Wait, e^{is} = zeta means is = ln(zeta), so isx has e^{isx} = zeta^x = e^{x*ln(zeta)}.
  And e^{isx} in the 3-Gamma lemma: a=b=c=3/4, shift by e^{isx} with s from e^{-s} = zeta^{-1}?

  Careful: the shifted 3-Gamma lemma is:
    int Gamma(a+ix) Gamma(b+ix) Gamma(c-ix) e^{isx} dx
    = 2pi e^{-as} Gamma(a+c) Gamma(b+c) U(a+c, a-b+1, e^{-s})

  Here a=b=3/4, c=3/4. We need int Gamma(3/4+ix)^2 Gamma(3/4-ix) e^{isx} dx.
  The term e^{isx} in the integrand, with the phase e^{is} = zeta = (1+iw)(1+iu)/((1-iw)(1-iu)).
  But actually in our double-GF we have e^{ix ln(zeta)} = e^{isx} where is = ln(zeta),
  i.e., s = -i ln(zeta) = -i*(ln(1+iw)+ln(1+iu)-ln(1-iw)-ln(1-iu)).
  Hmm, s is complex in general. But for the 3-Gamma lemma we need real s.
  However, the formula extends by analytic continuation.

  With e^{-s} = 1/zeta = (1-iw)(1-iu)/((1+iw)(1+iu)) (the inverse Cayley):

  int Gamma(3/4+ix)^2 Gamma(3/4-ix) e^{isx} dx = 2pi e^{-3s/4} Gamma(3/2)^2 U(3/2, 1, e^{-s})
  = 2pi e^{-3s/4} Gamma(3/2)^2 U(3/2, 1, 1/zeta)

  And [(1+w^2)(1+u^2)]^{-3/4} * e^{-3s/4}:
  e^{-3s/4} = (e^{-s})^{3/4} = (1/zeta)^{3/4} = zeta^{-3/4}
  (1+w^2) = (1-iw)(1+iw), (1+u^2) = (1-iu)(1+iu)
  (1+w^2)(1+u^2) = (1-iw)(1+iw)(1-iu)(1+iu)
  Let A=(1-iw)(1-iu), B=(1+iw)(1+iu), then zeta=B/A, (1+w^2)(1+u^2)=AB.
  zeta^{-3/4} = A^{3/4}/B^{3/4}
  [(1+w^2)(1+u^2)]^{-3/4} * zeta^{-3/4} = (AB)^{-3/4} * A^{3/4}/B^{3/4}
  = A^{-3/4}/B^{3/4} * A^{3/4}/B^{3/4} = ... wait let me be careful.

  sum J_mk w^m u^k = [(1+w^2)(1+u^2)]^{-3/4} * 2pi Gamma(3/2)^2 e^{-3s/4} U(3/2, 1, e^{-s})
  = [(1+w^2)(1+u^2)]^{-3/4} * 2pi Gamma(3/2)^2 zeta^{-3/4} U(3/2, 1, 1/zeta)
  = (AB)^{-3/4} * 2pi Gamma(3/2)^2 (A/B)^{3/4} U(3/2, 1, A/B)
  = A^{-3/4} B^{-3/4} * 2pi Gamma(3/2)^2 A^{3/4} B^{-3/4} U(3/2, 1, A/B)
  = B^{-3/2} * 2pi Gamma(3/2)^2 U(3/2, 1, A/B)
  = [(1+iw)(1+iu)]^{-3/2} * 2pi Gamma(3/2)^2 U(3/2, 1, (1-iw)(1-iu)/((1+iw)(1+iu)))

So GF = 2pi * (pi/4) * U(3/2, 1, 1/zeta) / [(1+iw)(1+iu)]^{3/2}
       = pi^2/2 * U(3/2, 1, 1/zeta) / [(1+iw)(1+iu)]^{3/2}

Note: argument is 1/zeta = (1-iw)(1-iu)/((1+iw)(1+iu)), NOT zeta!
At w=u=0: 1/zeta = 1. Good, same value.

The coefficient [w^m u^k] (1+iw)^{-3/2} = (-i)^m * C(m+1/2, m) = (-i)^m * Gamma(m+3/2)/(m! Gamma(3/2))
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

J00 = J(0, 0)
h0 = nrm[0]

# The CORRECT GF:
# sum J_mk w^m u^k = pi^2/2 * U(3/2, 1, 1/zeta) / [(1+iw)(1+iu)]^{3/2}
# where zeta = (1+iw)(1+iu)/((1-iw)(1-iu)), so 1/zeta = (1-iw)(1-iu)/((1+iw)(1+iu))

def GF_34_correct(w, u):
    """Correct GF: pi^2/2 * U(3/2, 1, 1/zeta) / [(1+iw)(1+iu)]^{3/2}."""
    inv_zeta = (1 - 1j*w) * (1 - 1j*u) / ((1 + 1j*w) * (1 + 1j*u))
    denom = ((1 + 1j*w) * (1 + 1j*u))**mp.mpf(3)/2
    return mp.pi**2 / 2 * mp.hyperu(mp.mpf(3)/2, 1, inv_zeta) / denom

# Verify at w=u=0: GF should = J_00
print(f"GF(0,0) = {mp.nstr(GF_34_correct(0, 0), 15)}")
print(f"J_00 = {mp.nstr(J00, 15)}")
print(f"pi^2/2 * U(3/2,1,1) = {mp.nstr(mp.pi**2/2 * mp.hyperu(mp.mpf(3)/2,1,1), 15)}")

# Extract Taylor coefficients via Fourier sums
print("\n--- GF verification (corrected) ---")
Mw, Mu = 64, 64
r = mp.mpf(1) / 4  # radius well inside singularity at |w|=1
ws = [r * mp.expjpi(2 * mp.mpf(a) / Mw) for a in range(Mw)]
us = [r * mp.expjpi(2 * mp.mpf(b) / Mu) for b in range(Mu)]
vals = [[GF_34_correct(ws[a], us[b]) for b in range(Mu)] for a in range(Mw)]

GF_coeffs = {}
for m in range(7):
    inner = [sum(vals[a][b] * mp.expjpi(-2 * mp.mpf(m*a) / Mw) for a in range(Mw)) / (Mw * r**m)
             for b in range(Mu)]
    for k in range(m, 7):
        J_gf = sum(inner[b] * mp.expjpi(-2 * mp.mpf(k*b) / Mu) for b in range(Mu)) / (Mu * r**k)
        GF_coeffs[(m, k)] = J_gf
        J_quad = J(m, k)
        phase = (-1j)**abs(m-k)
        diff = abs(J_gf - J_quad)
        if m <= 4 and k <= 4:
            print(f"  J_{m},{k}: GF={mp.nstr(mp.re(J_gf/phase), 8)}, "
                  f"quad={mp.nstr(mp.re(J_quad/phase), 8)}, "
                  f"diff={mp.nstr(diff, 3)}")

# Now find the right normalization for integrality.
# The GF is pi^2/2 * U(3/2, 1, q) / [(1+iw)(1+iu)]^{3/2}
# where q = (1-iw)(1-iu)/((1+iw)(1+iu)) is the inverse Cayley product.
#
# [w^m] 1/(1+iw)^{3/2} = (-i)^m * C(m+1/2, m) = (-i)^m * Gamma(m+3/2)/(m! Gamma(3/2))
# So the natural normalization to extract U coefficients is:
# T_mk = J_mk * m! Gamma(3/2)/Gamma(m+3/2) * k! Gamma(3/2)/Gamma(k+3/2) * (-i)^{m+k} / (pi^2/2)
# = [w^m u^k] U(3/2, 1, q)   ... if U(3/2,1,q) has integer Taylor coefficients around q=1.

# Actually: U(3/2,1,q) around q=1 is NOT a standard power series (it's analytic there),
# so write q = 1 - epsilon where epsilon = 1 - (1-iw)(1-iu)/((1+iw)(1+iu)).
# As w->0, u->0: epsilon ~ 2i(w+u) + ... and U(3/2,1,1-epsilon) = U(3/2,1,1) + higher.
# The Taylor coefficients of U(3/2,1,1-epsilon) in epsilon involve U and its derivatives,
# which are rational combinations of U(3/2,1,1) and U(5/2,2,1) (via the ODE).

# The key: generators are A' = U(3/2,1,1) and B' = U(5/2,2,1).
# Let's verify that c_mk entries lie in QA' + QB' by using a different normalization.

# For lambda=1/4 with h_n = h_0 * C(2n,n)/4^n:
# c_mk = J_mk/sqrt(h_m h_k) = J_mk/(h_0 sqrt(C(2m,m)C(2k,k)/16^{(m+k)/2}))
# T_mk^{1/4} = m! k! sqrt(C(2m,m)C(2k,k)) c_mk
#            = m! k! sqrt(C(2m,m)C(2k,k)) J_mk / (h_0 sqrt(C(2m,m)C(2k,k)/16^{...}))

# For lambda=3/4 with h_n = h_0 * (3/2)_n/n! where (3/2)_n = Gamma(n+3/2)/Gamma(3/2):
# c_mk = J_mk/(h_0 sqrt((3/2)_m (3/2)_k / (m! k!)))
# = J_mk * sqrt(m! k!) / (h_0 * sqrt((3/2)_m (3/2)_k))

# Try: T_mk^{3/4} = m! k! sqrt((3/2)_m (3/2)_k) * (-i)^{m+k} * J_mk / (pi^2/2)
# This makes T_mk^{3/4} "extracted from U(3/2,1,q)", possibly rational * (A', B') combination.

def Poch(a, n):
    """Pochhammer (rising factorial): (a)_n = Gamma(a+n)/Gamma(a)."""
    return mp.gamma(a + n) / mp.gamma(a)

Ap = mp.hyperu(mp.mpf(3)/2, 1, 1)
Bp = mp.hyperu(mp.mpf(5)/2, 2, 1)
print(f"\nA' = U(3/2,1,1) = {mp.nstr(Ap, 20)}")
print(f"B' = U(5/2,2,1) = {mp.nstr(Bp, 20)}")
print(f"B'/A' = {mp.nstr(Bp/Ap, 10)}")

pref = mp.pi**2 / 2
print("\n--- PSLQ with T_mk = m! k! Poch(3/2,m) Poch(3/2,k) * J_mk * (-i)^{m+k} / (pi^2/2) ---")
for m in range(6):
    for k in range(m, 6):
        Jmk = J(m, k)
        phase = (-1j)**(m+k)  # note: m+k parity
        T_raw = mp.fac(m) * mp.fac(k) * Poch(mp.mpf(3)/2, m) * Poch(mp.mpf(3)/2, k) * Jmk * phase / pref
        # Should be real
        T = mp.re(T_raw)
        if abs(mp.im(T_raw)) > 1e-20 * abs(T):
            print(f"  T_{m},{k} has Im = {mp.nstr(mp.im(T_raw), 5)} -- removing phase differently")
            T_raw2 = T_raw / (-1j)**abs(m-k) / (-1j)**(2*min(m,k))
            T = mp.re(T_raw2)
        r = mp.pslq([T, Ap, Bp], tol=1e-25, maxcoeff=100000)
        print(f"  T_{m},{k} = {mp.nstr(T, 10)}, PSLQ: {r}")

# Try another normalization: without the (3/2)_n factors
print("\n--- PSLQ without Pochhammer factors: J_mk * (-i)^{|m-k|} * (-1)^{?} / J_00 ---")
for m in range(6):
    for k in range(m, 6):
        Jmk = J(m, k)
        phase = (-1j)**abs(m-k)
        T = mp.re(Jmk / phase / J00)
        r = mp.pslq([T, Ap, Bp], tol=1e-25, maxcoeff=100000)
        print(f"  J_{m},{k}/J_00 = {mp.nstr(T, 10)}, PSLQ vs A',B': {r}")

# Also: try to find what the ODE for U(3/2,1,z) gives as generators
# Kummer's ODE: z U'' + (b-z) U' - a U = 0, here z U'' + (1-z) U' - 3/2 U = 0
# So U' = 3/2 U' ... wait: z y'' + (1-z) y' - 3/2 y = 0
# At z=1: y''(1) + (1-1)y'(1) - 3/2 y(1) = 0 => y''(1) = 3/2 y(1)
# With y'(1) = d/dz U(3/2,1,z)|_1 = -(3/2) U(5/2,2,1) = -(3/2) Bp
# So y(1) = Ap, y'(1) = -(3/2) Bp, y''(1) = 3/2 * Ap
# y'''(1): differentiate ODE once: z y''' + (1-z) y'' - (3/2+1) y' - y'' = 0
# z y''' + (1-z-1) y'' - 5/2 y' = 0 => z y''' - z y'' - 5/2 y' = 0
# At z=1: y'''(1) - y''(1) - 5/2 y'(1) = 0
# y'''(1) = y''(1) + 5/2 y'(1) = 3/2 Ap + 5/2*(-3/2 Bp) = 3/2 Ap - 15/4 Bp

# So all higher derivatives at z=1 are in Z*Ap + Z*Bp (rationally).
# Good! The Taylor series of U(3/2,1,1+delta) around delta=0 has coefficients in Q*Ap + Q*Bp.

# Now: q = (1-iw)(1-iu)/((1+iw)(1+iu))
# At w=u=0: q=1.
# delta = q - 1 = [(1-iw)(1-iu) - (1+iw)(1+iu)] / [(1+iw)(1+iu)]
# numerator = -2iw - 2iu = -2i(w+u)
# So delta = -2i(w+u) / [(1+iw)(1+iu)] + higher order
#
# The Taylor coefficients of delta in w,u involve (1+iw)^{-1}(1+iu)^{-1} factors.
# [w^m u^k] delta^j involves (1+iw)^{-j} etc.
# The whole thing is: [w^m u^k] U(3/2,1,q) = [w^m u^k] sum_j U^{(j)}(1)/j! * (q-1)^j
#
# The overall picture: T_mk are rational multiples of A' and B'.
# Let me look for the right denominators via the GF.

# The GF: sum J_mk w^m u^k = pi^2/2 * U(3/2, 1, q) / [(1+iw)(1+iu)]^{3/2}
# where q = (1-iw)(1-iu)/((1+iw)(1+iu))
#
# Change variables: let w = it/(1+t) (so that 1+iw = 1/(1-t)... no)
# Alternative: use generating variable for P_n.
# The key combinatorial step at lambda=1/4 was that the denominators 2^m m! absorbed
# the Taylor coefficients of 1/sqrt(1-iw) = sum_m (-1/2 choose m)(-iw)^m.
# The binomial coefficient (-1/2 choose m) = (-1)^m C(2m,m)/4^m.
# So [w^m] (1/(1-iw))^{1/2} = (i)^m C(2m,m)/4^m, and 4^m/(C(2m,m) m!) * (i)^m * C(2m,m)/4^m = i^m/m!.
# Thus (multiplying by m! 2^m = 4^m/C(2m,m) ... hmm.
# Actually for lambda=1/4: T_mk = m! k! sqrt(C(2m,m)C(2k,k)) c_mk.
# sqrt(C(2m,m)) = sqrt((2m)!)/(m!) ... this absorbs from h_n.
#
# For lambda=3/4: [w^m] (1/(1+iw))^{3/2} = (-i)^m (3/2)_m / m!
# The natural removal: multiply by m! / (3/2)_m = m! Gamma(3/2)/Gamma(m+3/2).
# Let beta_m = m! Gamma(3/2)/Gamma(m+3/2).
# T_mk^{3/4} = beta_m * beta_k * J_mk * (-i)^{m+k} * 2/(pi^2) = [w^m u^k] U(3/2,1,q)
#
# This is a Taylor expansion of U around q=1 (extracted from series).
# Let's verify this makes sense by extracting numerically.

print("\n--- Extract [w^m u^k] U(3/2,1,q) via Fourier ---")
def U_of_q(w, u):
    q = (1 - 1j*w) * (1 - 1j*u) / ((1 + 1j*w) * (1 + 1j*u))
    return mp.hyperu(mp.mpf(3)/2, 1, q)

# The Taylor coefficients of U_of_q should be in Q*A' + Q*B'.
Mw2, Mu2 = 32, 32
r2 = mp.mpf(1) / 8  # small radius
ws2 = [r2 * mp.expjpi(2 * mp.mpf(a) / Mw2) for a in range(Mw2)]
us2 = [r2 * mp.expjpi(2 * mp.mpf(b) / Mu2) for b in range(Mu2)]
vals2 = [[U_of_q(ws2[a], us2[b]) for b in range(Mu2)] for a in range(Mw2)]

print("Checking [w^m u^k] U(3/2,1,q) vs PSLQ A',B':")
for m in range(5):
    inner2 = [sum(vals2[a][b] * mp.expjpi(-2*mp.mpf(m*a)/Mw2) for a in range(Mw2)) / (Mw2 * r2**m)
              for b in range(Mu2)]
    for k in range(m, 5):
        coeff = sum(inner2[b] * mp.expjpi(-2*mp.mpf(k*b)/Mu2) for b in range(Mu2)) / (Mu2 * r2**k)
        # Remove phase: [w^m u^k] of U(q) where q has odd powers of iw,iu
        # q = 1 - 2i(w+u)/(1+iw)(1+iu) + ..., leading term ~ 1 - 2iw - 2iu
        # So [w^0 u^0] = U(3/2,1,1) = A', real
        # [w^1 u^0]: d/dw U(3/2,1,q)|_{w=u=0} = U'(3/2,1,1) * dq/dw|_0
        # dq/dw = -2i/((1+iw)(1+iu)) - (1-iw)(1-iu)*i(1+iu)/((1+iw)(1+iu))^2 ...
        # At w=u=0: dq/dw = -2i*1 + ... Actually let me compute numerically.
        T = mp.re(coeff) if (m+k)%2==0 else mp.im(coeff)
        rr = mp.pslq([coeff.real, Ap, Bp, coeff.imag], tol=1e-15, maxcoeff=1000)
        print(f"  [w^{m} u^{k}] U(q) = {mp.nstr(coeff, 8)}, PSLQ: {rr}")

# Try with rational scaling by known denominators
print("\n--- Try: is [w^m u^k] U(q) * (integer) in Z*A' + Z*B'? ---")
# Look at m=0 row:
# [u^0] U(q)|_{w=0} = U(3/2,1,1) = A'
# [u^1] U(q)|_{w=0} = dU/du|_0 = U'(3/2,1,1) * dq/du|_0
# dq/du at w=u=0: q = (1-iu)/(1+iu) at w=0, dq/du = (-i(1+iu)-(1-iu)*i)/(1+iu)^2 = -2i/(1+iu)^2|_0 = -2i
# So [u^1] U(q)|_{w=0} = -3/2 Bp * (-2i) = 3i Bp
# Thus [w^0 u^1] = 3i B' (purely imaginary, consistent with m+k=1 odd)
inner_m0 = [sum(vals2[a][b] * mp.expjpi(-2*mp.mpf(0*a)/Mw2) for a in range(Mw2)) / (Mw2 * r2**0)
            for b in range(Mu2)]
for k in range(6):
    coeff = sum(inner_m0[b] * mp.expjpi(-2*mp.mpf(k*b)/Mu2) for b in range(Mu2)) / (Mu2 * r2**k)
    print(f"  [w^0 u^{k}] U(q) = {mp.nstr(coeff, 10)}")
    print(f"    3i*B' = {mp.nstr(3j*Bp, 10)}" if k==1 else "")
