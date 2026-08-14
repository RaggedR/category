"""
Discrete Hermite matching matrix — the Krawtchouk analogue.

Continuous case (Langer–Claude 2026):
  I_{mk} = ∫ H_m(t) H_k(t) e^{-3t²/2} dt
  Generating function: √(2π/3) · exp(-(s² - 4su + u²)/3)
  Rescaled D_{mk} are integers counting signed bipartite matchings.

Discrete analogue:
  I(m,k,n,q) = Σ_x K_m(x;n) K_k(x;n) · C(n,x) q^x (1-q)^{n-x}
             = E[K_m(X;n) K_k(X;n)]  for X ~ Bin(n,q)

  Generating function: [1 + α(s+u) + su]^n   where α = 1-2q

  Coefficient extraction:
    I(m,k,n,q) = Σ_c n!/((n-m-k+c)!(m-c)!(k-c)!c!) · α^{m+k-2c}

  When α ∈ ℤ, every entry is automatically an integer.

This script explores the matrix for small n and various α.
"""
from fractions import Fraction
from math import comb, factorial


def I_entry(m, k, n, alpha):
    """
    Compute I(m,k,n,q) via coefficient extraction from [1 + α(s+u) + su]^n.

    I(m,k,n,q) = Σ_c n!/((n-m-k+c)!(m-c)!(k-c)!c!) · α^{m+k-2c}

    where c ranges over non-negative integers satisfying:
      n-m-k+c ≥ 0, m-c ≥ 0, k-c ≥ 0, c ≥ 0
    """
    if m > n or k > n:
        return Fraction(0)
    total = Fraction(0)
    c_min = max(0, m + k - n)
    c_max = min(m, k)
    for c in range(c_min, c_max + 1):
        a = n - m - k + c  # copies of "1"
        b = m + k - 2 * c   # copies of "α(s+u)"
        # multinomial: n! / (a! (m-c)! (k-c)! c!)
        # Note: b = m+k-2c and the binomial C(b, m-c) distributes s vs u
        # Combined: n! / (a! b! c!) * C(b, m-c) = n! / (a! (m-c)! (k-c)! c!)
        coeff = Fraction(factorial(n),
                         factorial(a) * factorial(m - c) * factorial(k - c) * factorial(c))
        total += coeff * Fraction(alpha) ** b
    return total


def print_matrix(n, alpha, label=""):
    """Print the full (n+1)×(n+1) matrix I(m,k,n,α)."""
    size = n + 1
    print(f"\n{'='*60}")
    print(f"n = {n}, α = {alpha}  {label}")
    print(f"Generating function: [1 + {alpha}(s+u) + su]^{n}")
    print(f"{'='*60}")

    # Compute matrix
    mat = [[I_entry(m, k, n, alpha) for k in range(size)] for m in range(size)]

    # Check if all entries are integers
    all_int = all(mat[m][k].denominator == 1 for m in range(size) for k in range(size))

    # Find max width for formatting
    strs = [[str(int(mat[m][k])) if mat[m][k].denominator == 1 else str(mat[m][k])
             for k in range(size)] for m in range(size)]
    max_w = max(len(s) for row in strs for s in row)

    # Header
    hdr = "m\\k" + "".join(f"{k:>{max_w+1}}" for k in range(size))
    print(hdr)
    print("-" * len(hdr))
    for m in range(size):
        row = f"{m:>3}" + "".join(f"{strs[m][k]:>{max_w+1}}" for k in range(size))
        print(row)

    print(f"\nAll entries integer: {all_int}")

    # Check symmetry
    sym = all(mat[m][k] == mat[k][m] for m in range(size) for k in range(size))
    print(f"Symmetric: {sym}")

    # Check parity
    zero_odd = all(mat[m][k] == 0 for m in range(size) for k in range(size) if (m+k) % 2 == 1)
    print(f"Vanishes when m+k odd: {zero_odd}")

    return mat


def krawtchouk_direct(j, x, n):
    """Krawtchouk polynomial K_j(x; n, 1/2) via explicit sum."""
    total = 0
    for r in range(j + 1):
        if r > x or j - r > n - x:
            continue
        total += (-1)**r * comb(x, r) * comb(n - x, j - r)
    return total


def verify_generating_function(n, alpha, num_q_vals=3):
    """Verify I(m,k,n,α) against direct Krawtchouk computation."""
    q = Fraction(1 - alpha, 2)
    if q < 0 or q > 1:
        print(f"  (q = {q} outside [0,1], skipping direct verification)")
        return True

    print(f"\nVerifying against direct Krawtchouk sum for n={n}, α={alpha} (q={q})...")
    ok = True
    for m in range(n + 1):
        for k in range(n + 1):
            # Direct computation
            direct = Fraction(0)
            for x in range(n + 1):
                km = krawtchouk_direct(m, x, n)
                kk = krawtchouk_direct(k, x, n)
                w = Fraction(comb(n, x)) * q**x * (1 - q)**(n - x)
                direct += km * kk * w
            # Generating function computation
            gf = I_entry(m, k, n, alpha)
            if direct != gf:
                print(f"  MISMATCH at (m,k)=({m},{k}): direct={direct}, gf={gf}")
                ok = False
    if ok:
        print(f"  All {(n+1)**2} entries match. ✓")
    return ok


# ─── Main exploration ───

print("DISCRETE HERMITE MATCHING MATRIX")
print("================================")
print()
print("Generating function: F(s,u) = [1 + α(s+u) + su]^n")
print("where α = 1-2q, q = parameter of Bin(n,q)")
print()
print("Coefficient extraction:")
print("  I(m,k,n,α) = Σ_c n!/((n-m-k+c)!(m-c)!(k-c)!c!) · α^{m+k-2c}")

# Verify the generating function formula
verify_generating_function(6, Fraction(1, 3))  # q = 1/3, α = 1/3

# Small examples with integer α
for n in [4, 6, 8]:
    for alpha in [1, -1, 2, -2, 3]:
        q = Fraction(1 - alpha, 2)
        label = f"(q = {q})"
        print_matrix(n, alpha, label)

# The orthogonality case
print("\n" + "="*60)
print("ORTHOGONALITY CHECK: α = 0 (q = 1/2)")
print("="*60)
mat = print_matrix(6, 0, "(q = 1/2)")
print("\nDiagonal entries should be C(n,k):")
for k in range(7):
    print(f"  I({k},{k}) = {mat[k][k]}, C(6,{k}) = {comb(6,k)}")

# Rational α — look for structure
print("\n" + "="*60)
print("RATIONAL α: looking for non-integer α with integer rescaling")
print("="*60)

for n in [6]:
    for alpha_num, alpha_den in [(1,3), (2,3), (1,2), (3,2)]:
        alpha = Fraction(alpha_num, alpha_den)
        mat = print_matrix(n, alpha, f"(q = {Fraction(1-alpha, 2)})")

        # Try rescaling by α^{-(m+k)} to see if structure emerges
        print(f"\n  Rescaled by {alpha_den}^(m+k):")
        for m in range(n + 1):
            vals = []
            for k in range(n + 1):
                rescaled = mat[m][k] * Fraction(alpha_den) ** (m + k)
                vals.append(str(int(rescaled)) if rescaled.denominator == 1 else str(rescaled))
            print(f"    m={m}: {vals}")


# Diagonal analysis
print("\n" + "="*60)
print("DIAGONAL ANALYSIS: I(k,k,n,α) for integer α")
print("="*60)
print("\nα=2:")
for n in range(1, 16):
    diag = [I_entry(k, k, n, 2) for k in range(n + 1)]
    print(f"  n={n:2d}: {[int(d) for d in diag]}")

print("\nα=1:")
for n in range(1, 16):
    diag = [I_entry(k, k, n, 1) for k in range(n + 1)]
    print(f"  n={n:2d}: {[int(d) for d in diag]}")

print("\nDiag I(0,0,n,α) = [1 + 0 + 0]^n = 1 for all n (trivially).")
print("Diag I(1,1,n,α) = α² · n + n = n(α²+1)  (verify below):")
for alpha in [1, 2, 3]:
    for n in [4, 8, 12]:
        val = int(I_entry(1, 1, n, alpha))
        expected = n * (alpha**2 + 1)
        print(f"  α={alpha}, n={n}: I(1,1) = {val}, n(α²+1) = {expected}, match={val==expected}")


# Edge entries
print("\n" + "="*60)
print("EDGE ENTRIES: I(0,k,n,α)")
print("="*60)
print("These should be [s^0 u^k] [1+α(s+u)+su]^n = C(n,k) · α^k")
for alpha in [1, 2, -1, 3]:
    print(f"\nα = {alpha}:")
    for n in [6, 8, 10]:
        edges = [int(I_entry(0, k, n, alpha)) for k in range(n + 1)]
        expected = [comb(n, k) * alpha**k for k in range(n + 1)]
        match = edges == expected
        print(f"  n={n}: {edges}")
        print(f"       = C(n,k)·α^k: {match}")


# Compare continuous vs discrete structure
print("\n" + "="*60)
print("STRUCTURAL COMPARISON: continuous vs discrete")
print("="*60)
print("""
Continuous (Hermite):
  GF = exp(-(s² - 4su + u²)/3)
  = exp of quadratic form
  Rescaled entries D_{mk} count signed bipartite matchings
  Parity: vanishes when m+k odd
  Edge: D_{0,2n} = (-1)^n (2n)!/n!
  Diagonal: A081671, OGF = 1/√(1-8x+12x²)

Discrete (Krawtchouk):
  GF = [1 + α(s+u) + su]^n
  = n-th power of (1 + linear + bilinear)
  Entries I(m,k) = multinomial sums × α^{m+k-2c}
  Integer when α ∈ ℤ (automatic!)
  Parity: does NOT vanish when m+k odd (for α ≠ 0)
  Edge: I(0,k) = C(n,k) · α^k
  Diagonal: I(k,k,n,α) = ? (to investigate)
""")
