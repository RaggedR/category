# Plan: language.tex — The Dictionary Paper

*Handover for the next agent. Written 2026-08-11 after commit 64cb017 (MP
ℤA⊕ℤB integrality proved via Cayley factorization + Kummer recurrence).*

First: read `~/.claude/AGENT.md`, then the memory index at
`~/.claude/projects/-Users-robin-git-category/memory/MEMORY.md`
(especially `project_matching_synthesis.md` and `project_mp_dual_pair.md`).

## What this paper is

NOT a survey of the seven matching-matrix computations — those live in
`matching-coefficients.tex` and the standalone pdfs. This paper's object
is the **dictionary itself**: a three-column language (analysis /
combinatorics / arithmetic) in which the family results are entries.
The claim, agreed with Robin, is that the individual theorems keep
deflating to classical mechanisms, but the dictionary that makes the
questions askable is new. Model: the Askey scheme contains almost no
theorems and is the most-used object in the field. Write the paper so
the *introduction claims the dictionary, not the entries*.

## The dictionary (this is the paper's centrepiece)

| analysis | combinatorics | arithmetic |
|---|---|---|
| GF kernel | automaton | generator lattice |
| ODE order of kernel | number of states | lattice rank |
| Taylor coeffs uₙ at base point | walk weights | all signs + generators |
| substitution (Cayley/Laguerre σ) | species chassis: defects, chains, matchings | universal, sign-free; denominators clear by Legendre |
| kernel singularity type | — | obstruction level |
| D-module reducibility | — | algebraic vs transcendental generators |

Sentences only this language can say (the paper's §"claims"):
1. Gaussian rank 1; MP rank 2; Hahn rank 2 with one algebraic generator
   (reducible D-module, Jordan–Hölder visible in the lattice); Wilson no
   finite rank. No classical reference makes this classification.
2. The combinatorics is universal; the family only chooses uₙ
   (proved for MP¼ fwd+rev; the mechanism of commit 64cb017).
3. Degeneration is evaluation, not limit (Jacobi rationality trick).
4. Obstruction hierarchy: entire weight → geometric decay + ℤ;
   branch point at finite distance → power-law decay; Barnes Γ(2ix) →
   no lattice at all.

## Family data (do not re-derive; sources in brackets)

| family | GF kernel | order | lattice / generators | status |
|---|---|---|---|---|
| Hermite | exp(quadratic) | 1 | ℤ; signed bipartite matchings; diag A081671 | proved [hermite-matching-matrix] |
| Laguerre | 2/(1+w+u−3wu) | 1 | ℤ; Delannoy | proved [laguerre-matching] |
| Charlier | Charlier kernel | 1 | ℤ; automaton thm | proved [charlier-matching] |
| Chebyshev/Jacobi | Wallis/Beta mass | — | quartic classes P,Q,P̃,Q⁺ over J(0)=Γ(¼)²/√(2π); sextic ν₀=B(⅙,½); connection-Gram | proved [jacobi-matching §7] |
| MP¼ fwd | e^{ζ/2}K₀(ζ/2)=√πU(½,1,ζ) | 2 | ℤA⊕ℤB, A,B=√(2e)K₀,K₁(½) | **proved 2026-08-11** [meixner-pollaczek-simple Thm 3.1+4.3, mp_symbolic_lattice.py] |
| MP¼ rev | e^{−1/ζ} | 1 | (2π/e)ℤ; defect/chain/Charlier species | proved [ibid §6] |
| MP¾ | Tricomi | 2 | A′=U(3/2,1,1), B′=U(5/2,2,1) | empirical; **predicted provable by the same Kummer mechanism** [mp34-matching] |
| unbalanced Γ¼Γ¾ | first-order | 1 | single A″=eE₁(1); one-state automaton | proved [unbalanced-dual] |
| Hahn | √ζ′·U(1,1,·)-type | 2 | ℚ𝒜⊕ℚℬ, 𝒜=√πe²E₁(2), ℬ=√π; reducible module ⇒ ℬ algebraic | empirical (PSLQ) [continuous-hahn §4.3] |
| Wilson | Barnes Γ(2ix) | ≥3? | NO compression to 1e18; moments ℚ-rank ≥16 | open [wilson_rank_test_hp] |

Also available if wanted: Weil-form results (arch side rank 2 over
√(2π), γ+log8π; prime side rank ≥16 — "the explicit formula's two
halves are algebraically incommensurable in this basis")
[weil-form.tex]. Probably a remark, not a section — ask Robin.

## Suggested structure (Robin will cut; draft minimal)

1. **The dictionary** — the table + the four sentences, immediately.
2. **The chassis** — one page: W → P̂ → Φ; forward approximation /
   reverse expansion; the universal species (defects (2d−1)!!, Lah
   chains, matchings) and why Legendre's formula clears denominators.
   Cite the MP¼ proof as the worked instance of "combinatorics
   universal, arithmetic in uₙ".
3. **The entries** — one short paragraph per family, data table only,
   pointing to standalone pdfs. No proofs.
4. **The obstruction hierarchy** — decay/lattice axis with the three
   levels and which family sits where.
5. **The test** — Wilson as the falsifiable prediction of the
   language; state the D-module-rank question precisely. If the
   language can't make Wilson well-posed, say so honestly.
6. **Status table** — proved / empirical / predicted, per claim.
   This paper must not blur the three.

## ⚠ Pitfalls

1. **Honesty layering.** Hahn integrality is PSLQ-empirical (and the
   original "PSLQ failed" verdict was a maxcoeff=1000 artifact —
   always maxcoeff ≥ 1e10, confirm at 2× dps). MP¾ is *predicted*, not
   proved. Wilson "rank ≥16" comes from the moment ceiling argument.
   Do not upgrade any of these in prose.
2. **Robin's conventions.** Write HALF the length you think is right —
   he edits by deletion. No colons mid-sentence in prose. Real theorem
   environments with existing statements verbatim. Forward =
   "approximation", reverse = "expansion". Unicode in conversation,
   LaTeX only in .tex. tectonic to compile (not pdflatex).
3. **Don't recompute.** Every number in the family table is
   quadrature-verified in the standalone papers/scripts. If you must
   verify, theta-substitution for endpoint singularities, mp.dps ≥ 40.
   Never sympy nsimplify on exact rationals.
4. **The Askey-scheme framing is Robin-approved in conversation** (this
   session) but he has not seen it in print — flag it as a framing
   choice in your first draft message, don't bury it.
5. **OEIS** — JSON API oeis.org/search?q=...&fmt=json; bare `null`
   means no match. New sequences from 64cb017 (uₙ parts, MP diagonal)
   are all unmatched; citable as such.

## Optional pre-work with upside (do BEFORE drafting if energy allows)

- **MP¾ via the same mechanism** (~an afternoon): Tricomi kernel is
  also Kummer-type; if u′ₙ obeys an integer recurrence the "predicted"
  row upgrades to "proved" and claim 2 gets a second instance. This
  strictly strengthens the paper.
- **Wilson kernel ODE order** (harder, unbounded): even a lower bound
  ≥3, or evidence of irregularity from Γ(2ix), converts §5 from
  question to partial result. Do not sink more than a session; the
  paper stands without it.

## Success criteria

- language.tex compiles clean with tectonic → language.pdf.
- Every claim tagged proved/empirical/predicted, consistent with the
  status table; nothing upgraded.
- No proofs duplicated from standalone papers; pointers only.
- Draft is SHORT (target ≤ 8pp before Robin's cut).
- Commit to ~/git/category main; update
  `project_matching_synthesis.md` + MEMORY.md index.
