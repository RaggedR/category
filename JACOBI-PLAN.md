# Plan: General Jacobi Matching Matrices

*Handover for the next agent. Written 2026-08-10 after commits e95fe23 (phase-error fix) and 04f2812 (U + asymmetric integrality).*

First: read `~/.claude/AGENT.md`, then the memory index at
`~/.claude/projects/-Users-robin-git-category/memory/MEMORY.md`
(especially `project_matching_synthesis.md`). The papers are
`matching-coefficients.tex` (THE synthesis paper, §2.2 Chebyshev) and
`jacobi-matching.tex` (standalone, corrected).

## The construction (one sentence each)

Kernel W → orthonormal polynomials P̂ₙ → functions Φₙ = P̂ₙ√W.
Forward matching coefficients bₘₖ = ∫P̂ₘP̂ₖW^{3/2}dx (these *approximate*
Φₘ by polynomials — say "approximation", never "expansion", for the forward
direction). Reverse coefficients use W^{1/2}.

## Established results (all quadrature-verified, do not re-derive from scratch)

Quartic factorials: P(n)=1·5···(4n−3), Q(n)=3·7···(4n−1),
P̃(n)=5·9···(4n+1), Q⁺(n)=7·11···(4n+3). Write s=(m+k)/2, d=|m−k|/2.

| matrix | kernel | rescale by | entry | edge | diagonal |
|---|---|---|---|---|---|
| forward T | (1−x²)^{−1/2} | Q(s)/J(0) | ½[P(d)Q(s)/Q(d)+P(s)] ∈ ℤ₊ | P(n) = A007696 | ½[Q+P] = 1,2,13,138 ∉ OEIS |
| reverse T | same, W^{1/2} | P̃(s)/J̃(0) | ∈ ℤ, signs alternate | −Q(n−1) | ½[P̃−Q(n−1)] = 1,2,21,282 ∉ OEIS |
| forward U | (1−x²)^{+1/2} | 3Q⁺(s+1)/J(0) | ½[Q⁺(s+1)ρ(d)+3P(s)] ∈ ℤ, ρ(n)=−3P(n−1)/Q⁺(n) | −15P(n−1) | ½[Q⁺(n+1)+3P(n)] = 5,40,585 ∉ OEIS |
| asymmetric (½,−½) | (1−x)^{1/2}(1+x)^{−1/2} | 2^k√2 | ℤ[½]·√2, whole matrix dyadic | row 0: 2C(2n,n) even, −(4n+3)Cat(n) odd | — |

Generator constant J(0) = √π Γ(¼)/Γ(¾) = Γ(¼)²/√(2π) — the lemniscatic
period (CM curve y² = 1−x⁴). All odd-parity entries of symmetric kernels
vanish. Odd entries of the asymmetric kernel do NOT vanish.

Key mechanism: x=cosθ collapses TₘTₖ to two cosines; the Wallis integral is
𝓘(p,q) = ∫₀^π sinᵖθ cos(qθ)dθ = **cos(πq/2)** · πΓ(p+1) / (2ᵖ Γ((p+q+2)/2) Γ((p−q+2)/2)).

## ⚠ Pitfalls (each one burned a previous session)

1. **The phase trap.** The cos(πq/2) factor in the Wallis formula is
   essential. Dropping it creates phantom odd-parity entries and sign errors
   that *look like structure* (a fake second generator γ_o lived in
   jacobi-matching.tex for days). Always verify closed forms against
   `mp.quad` of the *defining integral*, never against a re-derivation of the
   same combinatorics.
2. **Quadrature at endpoint singularities.** x-space quadrature of
   (1−x²)^{−3/4} is unreliable near ±1. Substitute θ: integrate
   `mp.quad(lambda t: ..., [0, mp.pi])` at `mp.dps = 40`, splitting at
   π/4, π/2, 3π/4 for stubborn cases.
3. **Rational recognition.** Use `Fraction(float(v)).limit_denominator(10**6)`
   and then assert agreement to ≥ 10 digits; for two-generator lattices use
   mpmath PSLQ with `maxcoeff ≥ 1e10`.
4. **Robin's conventions.** Real theorem environments with his statements
   verbatim; "triangular" matrix (not "symmetric") if he says so; forward =
   approximation, reverse = expansion; no colons mid-sentence in prose;
   unicode math in conversation, LaTeX only in .tex.

## The central conjecture (why this matters)

**Matching arithmetic = period arithmetic.** The generator ring of the
matching matrix appears to be the period lattice of the weight:
Gaussian → ℤ; Chebyshev → quartic Pochhammers (¼)ₙ, (¾)ₙ over the
lemniscatic period Γ(¼)²; asymmetric half-integer → dyadic Catalan/central
binomial (ballot paths); MP → Bessel-period lattice ℤA⊕ℤB.

## ✅ Session results 2026-08-10 — sextic test DONE, conjecture graduates (refined)

Scripts + logs in repo: `jacobi_sextic.py` (exact Gram–Schmidt, 40-digit
verification of every entry against the defining integral),
`jacobi_sextic_analysis.py` (primitive-ℤ frame, integrality),
`jacobi_sextic_connection.py` (the structural theorem). α = β = −5/9.

1. **Sextic generator confirmed.** ν₀ = B(½,⅙) = 2^{2/3}πΓ(⅓)/Γ(⅔)²
   (plan's evaluation checked, 40 digits). General ratio, a=(p+2)/2:
   𝓘(p,2n)/𝓘(p,0) = (1−a)ₙ/(a)ₙ; p=−⅔ gives (⅓)ₙ/(⅔)ₙ — thirds.
2. **Two-lattice refinement of the central conjecture.** The entry
   integral (W^{3/2}) contributes the sextic/thirds lattice; the
   *polynomials themselves* (weight W) contribute a second lattice at
   eighteenths: monic Gegenbauer λ=−1/18 recurrence
   βₙ = 9n(9n−10)/((18n−1)(18n−19)) (verified n≤8), so norms live on
   ∏(9j−10), ∏(18j±1) classes. Chebyshev's single quartic lattice is the
   λ=0 degeneration (Tₙ∈ℤ[x], β=¼). Matrix = (eighteenths)ᵀ(thirds)(eighteenths).
3. **Integrality.** In the primitive-ℤ polynomial frame (PZₙ, contents =
   sliding windows of (18j−1) primes), den(Eₘₖ) = 2^a × primes from
   {5,11,17,23,…} = **6j−1 — a genuine sextic factorial family**, depending
   essentially on s=(m+k)/2 only. So an s-dependent rescale gives ℤ entries
   (same shape as Chebyshev's Q(s)); there is no clean pure-row factor.
4. **Structural theorem (the real prize; verified ALL m,k≤8 exactly).**
   With λ = α+½ (polys), μ = 3α/2+½ (i.e. W^{3/2} = (1−x²)^{μ−½}), the
   classical Gegenbauer connection formula gives every entry as a finite
   Pochhammer sum — the matching matrix is the **Gram matrix of connection
   coefficients** C^λ → C^μ in the μ-norm:
   ⟨πₘ,πₖ⟩_{W^{3/2}} = mfₘ mfₖ Σⱼ A_{(m−j)/2}(m) A_{(k−j)/2}(k) hⱼ^μ,
   A_k(n) = (λ−μ)ₖ(λ)_{n−k}(n−2k+μ)/(k!(μ)_{n−k+1}),
   hⱼ^μ/h₀^μ = (2μ)ⱼ μ/(j!(j+μ)), mfₙ = n!/(2ⁿ(λ)ₙ), and λ−μ = −α/2
   for every symmetric Jacobi weight. Row-0 closed form:
   q_{0,2n} = (2n)!(λ−μ)ₙ(λ)ₙμ/(2^{2n}(λ)_{2n} n!(μ)ₙ₊₁) — the (5/18)ₙ
   factor explains the observed primes 5,23,41,59 (≡5 mod 18).
   The Chebyshev two-term closed form ½[P(d)Q(s)/Q(d)+P(s)] should be the
   λ→0 limit. This is a *proof mechanism*: general (α,β) follows from
   Jacobi connection coefficients (also classical, single-sum).

**New pitfall (add to the list):** sympy `nsimplify` numerically re-searches
even exact Rationals and can return surd garbage (e.g. 621/14840 →
225·2^{78/145}·…). Never call nsimplify on exact quantities.

### Next steps for the next agent (revised order)

1. ✅ **DONE 2026-08-10 — Write it up.** New §7 "The Sextic Test" in
   `jacobi-matching.tex` (Wallis-ratio prop, ν₀ generator, two lattices,
   E-table, connection-Gram theorem + row-0 corollary, λ→0 corollary);
   abstract item 5, Discussion + open-problem updates. The λ→0 limit is
   NOT analytic — every monic entry is a *rational function of α* (the
   1/(λ)ₙ of the monic factor cancels the (λ)_{n−ℓ} zeros), so Chebyshev
   is evaluation at α=−½. Verified all m,k≤8 exactly:
   `jacobi_lambda0_limit.py` / `.log` (also checks the integer two-term
   form and the degenerate Tₙ-in-C^{−1/4} connection coefficients).
   Synthesis paper §2.2 NOT touched (Robin's choice).
2. **Prove general symmetric α symbolically** (sympy, symbolic α): the
   connection-Gram identity is classical, so the content is the resulting
   generator-ring statement: entries ∈ ℚ·ν₀ in any rational-normalized
   frame, ν₀ = B(½, (3α+2)/2·…) — write the beta value once and check the
   task-5 prediction for which α give CM periods (thirds ↔ α ∈ −⅓ℤ-ish
   families; the Pochhammer steps are (−α/2)ₙ, (λ)ₙ, (μ)ₙ, (2μ)ₙ).
3. Then original Tasks 2–4 below (Wallis-power-zero family, reverse U /
   asymmetric, ballot-path model). For Task 3 note the reverse kernel has
   its own μ′ = α/2+½, λ−μ′ = +α/2·(−1)… — same mechanism, sign flipped.
4. Task 5 (general (α,β)) now has a mechanism: use Jacobi connection
   coefficients instead of PSLQ guessing.

## Tasks, in order (original — Task 1 done, see above)

### 1. The sextic test (decisive, do first)

Prediction: a symmetric Jacobi weight whose Wallis power p = 3α+1 has
(p+2)/2 ∈ ⅓+ℤ should produce **sextic** structure — Γ(⅓)-type generator
(period of y² = 1−x⁶, CM by ℤ[ω]) and Pochhammer ratios stepping by
thirds. Candidate: α = β = −5/9 gives p = −⅔, so
𝓘(p,0) = πΓ(⅓)/(2^{−2/3}Γ(⅔)²) — check this evaluation yourself.
- The polynomials are Gegenbauer C^{(λ)}, λ = α+½ = −1/18; they do NOT
  collapse to single cosines. Expand them in the Chebyshev basis (exact
  rational coefficients via recurrence) so every entry becomes a finite
  ℚ-combination of 𝓘(−⅔, 2n).
- First derive 𝓘(p,2n)/𝓘(p,0) as a ratio of shifted Pochhammers for
  general p (reflection formula; watch the phase-sign cancellation).
- Compute rows 0–4 exactly. Question: is there a single row factor making
  entries integral? Does a "sextic factorial" family appear?
- If yes: the conjecture graduates. If no: report *how* it fails — a
  rational-but-not-integral answer is still informative.

### 2. The Wallis-power-zero family (cheap, surprising if true)

(α,β) = (a,−a) always has p = 3(a−a)/2·… — check: does every member of
this one-parameter family have Γ-free (rational or algebraic) entries, or
is a = ½ special because the polynomials are trigonometric (3rd/4th-kind
Chebyshev)? Test a = ¼ numerically. Also do the mirror case
(−½, ½) (third-kind Vₙ, cos((n+½)θ)/cos(θ/2)) — expect the same
sequences up to signs; confirm.

### 3. Reverse matrices for U and asymmetric (bookkeeping + one honest check)

Reverse-U uses W^{1/2} = (1−x²)^{1/4}, effective p = ½, whose Pochhammer
ratios are (−¼)ₙ/(5/4)ₙ — this REUSES classes P and P̃. That would break
the "four classes each exactly once" sentence at the end of §2.2 of
matching-coefficients.tex. Compute it, and if the reuse is real, soften
that sentence (it currently refers only to the four matrices computed
so far — make the scope explicit or restate the poetry honestly).
Reverse asymmetric: expect dyadic again; find the row-0 combinatorics.

### 4. A ballot-path model for the Catalan row (the deep one)

Row 0 of the asymmetric kernel is 2C(2n,n) / −(4n+3)Cat(n). Catalan
numbers demand a bijection. Find a weighted-lattice-path or matching model
whose transfer matrix reproduces the full dyadic matrix, in the style of
§4 (Automata) of the synthesis paper — states, transitions, weights,
start/end rule. Search first (Robin's rule): fourth-kind Chebyshev moment
combinatorics, Viennot/Flajolet weighted Motzkin paths for W = tan(θ/2)
weights, "central binomial + Catalan interleaved" in OEIS/literature.

### 5. Only then: general (α,β)

With data points at exponents {0, ±½, −5/9}, guess the general law for
the generator ring as a function of (α,β) (expect: generated by the beta
value B(3α/2+…, 3β/2+…), i.e. a period of y^N = xᵃ(1−x)ᵇ). State it as a
conjecture in the synthesis paper only when at least three independent
exponents confirm.

## Success criteria

- Every claimed closed form checked against 40-digit quadrature of the
  defining integral for m,k ≤ 8 before it enters any .tex file.
- New sequences checked against OEIS (JSON API: oeis.org/search?q=...&fmt=json;
  a bare `null` response means no match).
- tectonic compiles clean; commit to ~/git/category (repo exists, main branch).
- Update `project_matching_synthesis.md` memory + MEMORY.md index.
- If an approach fails three times, stop and write down why (three-strike rule).
