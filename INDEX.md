# Category — Index

Research on orthogonal polynomial matching coefficients, the Riemann zeta function, quantum cohomology, and category theory.

## Directory Structure

### [`matching/`](matching/) — Matching Coefficients

The core research programme. Given a weight function w(x) and its orthogonal polynomials {P_n}, the *matching coefficients* b_{N,k} approximate w(x) as a finite sum ∑ b_{N,k} P̂_k(x), and the *matching matrix* M_{n,k} encodes how the first N polynomial sections converge to the function.

| Subfolder | Paper | Weight / family |
|-----------|-------|-----------------|
| — | **matching-coefficients.tex** | Synthesis paper (7 families) |
| — | **language.tex** | Dictionary of matching-coefficient concepts |
| `hermite/` | hermite-matching-matrix.tex | e^{−x²} — Hermite polynomials |
| `laguerre/` | laguerre-matching.tex | x^α e^{−x} — Laguerre polynomials |
| `charlier/` | charlier-matching.tex | e^{−μ} μ^x / x! — Charlier (discrete) |
| `mp/` | meixner-pollaczek.tex (+ 2 variants) | Meixner–Pollaczek at λ = ½ |
| `mp34/` | mp34-matching.tex | Meixner–Pollaczek at λ = ¾ |
| `joint-sech/` | joint-sech.tex | Joint ¼ ⊗ ¾ sech weight |
| `unbalanced/` | unbalanced-dual.tex | Unbalanced Γ(¼)Γ(¾) dual |
| `hahn/` | continuous-hahn.tex | Continuous Hahn polynomials |
| `jacobi/` | jacobi-matching.tex | Jacobi polynomials (endpoint singularities) |
| `wilson/` | wilson-matching.tex | Wilson polynomials |
| `plans/` | orthogonal-plan.tex | Roadmap: kernels and dual bases |

Supporting files: each subfolder contains `.py` scripts (numerical computation), `.pkl` caches, and `.png` plots.

### [`zeta/`](zeta/) — Riemann ξ and ζ

Papers on polynomial approximation of the Riemann Ξ function and related analytic number theory.

| Paper | Topic |
|-------|-------|
| zero.tex | Zeros of polynomial sections of Hermite functions |
| zeta.tex | Polynomial approximation of Ξ(t) |
| zeta-basis.tex | A family of zeta functions and their inner product |
| mellin.tex | Log-Hermite functions, Mellin transform, local RH |
| simple.tex | The zeta family, simply |

### [`cohomology/`](cohomology/) — Quantum Cohomology & Weil

The matching matrix reinterpreted through mirror symmetry: ℤA ⊕ ℤB as periods of the mirror of ℙ¹.

| Paper | Topic |
|-------|-------|
| cohomology.tex | Geometry behind the Bessel lattice / QH(ℙ¹) |
| weil-form.tex | Weil explicit formula in the MP basis |
| negative-result.tex | The prime side does not compress |

### [`category-theory/`](category-theory/) — Pure Category Theory

| Item | Topic |
|------|-------|
| two-colour-trees.tex | Two-coloured trees, container composition |
| limits-colimits-rep.tex | Limits and colimits in Rep(G) |
| coloured-graphs-as-directed-containers.pdf | Coloured graphs as directed containers |
| snoc-trees.pdf | Snoc trees |
| hopf-torus.html | Interactive Hopf–torus visualisation |
| `conway/` | Game of Life via store comonad |
| `knowledge-graph/` | Coloured graphs paper (extended) + Alastair collaboration |

### [`schur-weyl/`](schur-weyl/) — Schur–Weyl Duality

Four papers on explicit Schur–Weyl duality, Iwahori–Hecke algebras, Mackey theory, and local rules for Young tableaux.

### [`topology-debate/`](topology-debate/) — LLM Topology Debate

Experiment: Claude instances debating point-set topology, with evaluation framework and findings.

### [`linkedin/`](linkedin/) — LinkedIn Posts

Category theory posts: Connes–Kreimer, phantom types, Heegaard splitting, F₁, Hopf monads, Lean verification, Haskell–Hopf algebras.

### [`blog/`](blog/) — Blog Posts & Essays

- sn-is-gln1-blog.md — S_n is GL_n(𝔽₁)
- greens-polynomials.md — Green's polynomials
- s-expression-problem.md — The expression problem
- starting-simple-ai-tax.pdf — Starting simple (AI tax)
