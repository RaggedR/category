In 1955, Green startled the mathematical community by publishing the complete, general formula for the complex character tables of all finite general linear groups GL_n(q).

GL_n(q) is the group of invertible n×n matrices over the finite field 𝔽_q. Its size depends on q, its structure depends on q, and nobody expected a uniform answer across all q. Green found one anyway.

The polynomials Q^λ_μ(q) are indexed by two partitions of n: λ specifies a unipotent conjugacy class, μ specifies a maximal torus. They don't give the full character table directly — the irreps of GL_n(𝔽_q) are indexed by a much richer set than single partitions. But they give the essential building block: the unipotent slice, from which the full table can be assembled.

The Change of Basis

Green's polynomials are transition coefficients between two bases of the symmetric function ring Λ: p_μ = Σ_λ Q^λ_μ(q) · P_λ(x; q⁻¹), where p_μ are power-sum symmetric functions and P_λ(x; t) are Hall–Littlewood symmetric functions — a one-parameter family interpolating between Schur functions at t = 0 and monomial symmetric functions at t = 1.

Three Miracles

First: Q^λ_μ(q) ∈ ℤ[q]. Always. You're doing linear algebra over ℚ(q) and the denominators miraculously cancel. This integrality is what lets you plug in any prime power and get honest integers.

Second: Green's polynomials factor. The characters of the symmetric group are themselves a change of basis: p_μ = Σ_λ χ^λ_μ s_λ. The Kostka–Foulkes polynomials provide a second: s_λ = Σ_ν K_λν(q) P_ν(x; q). Composing the two gives Q = χᵀ · K, the character table of Sₙ times the Kostka–Foulkes matrix. The Kostka–Foulkes polynomials have coefficients in ℤ_{≥0} — strictly non-negative. For twenty years nobody could explain why. Then Lascoux and Schützenberger (1978) showed each coefficient counts semistandard Young tableaux with a given value of a statistic called "charge."

Third: set q = 1 and the finite field vanishes. Hall–Littlewood functions collapse to monomial symmetric functions, and Green's polynomials become ordinary character values of the symmetric group: Q^λ_μ(1) = χ^λ_μ. Green built a q-deformation where Sₙ is the ground state and GL_n(𝔽_q) is the excited state. One polynomial, infinitely many groups.

This is more than analogy. Tits proposed that Sₙ literally *is* GL_n(𝔽₁) — the general linear group over "the field with one element." 𝔽₁ doesn't exist as a field in the classical sense, but set q = 1 everywhere and linear algebra becomes combinatorics: a k-dimensional subspace of 𝔽_q^n — q^k vectors closed under linear combination — becomes a k-element subset, and the Gaussian binomial [n choose k]_q becomes C(n,k). The q = 1 specialisation isn't a trick but evidence of a geometry over a field that doesn't exist yet.

This is magic as far as I am concerned

#Mathematics
