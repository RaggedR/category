# Star Topology Results: 6-Agent Analysis of Purbhoo Migration Bug

## Summary

The star topology produced HIGHER diversity (0.450) than the ring (0.318), **reversing** the GECCO prediction. But this is not a refutation — it reveals that the GECCO ordering measures RAW diversity, while for bug diagnosis, convergence quality matters more than diversity quantity.

## Agent Outputs

| Agent | Role | Words | Position |
|-------|------|-------|----------|
| Vanilla | Generalist center (hub) | 1,717 | Center |
| Exegete | Paper parsing | 2,174 | Peripheral |
| Coords | Coordinate translation | 1,798 | Peripheral |
| Orient | Orientation analysis | 1,938 | Peripheral |
| Boundary | Domain & termination | 2,062 | Peripheral |
| Tracer | Figure 2 empirical | 2,310 | Peripheral |

## Key Findings

### 1. The direction rule is NOT the (remaining) bug (NEW finding)

The Orient agent discovered that `_pivot_vertex` only produces spokes {2, 3, 4}. Spokes 0, 1, 5 are unreachable. For these three reachable spokes, the existing direction logic is correct. The collision bug must be elsewhere.

### 2. `ensure_hex_cell` causes overshooting (NEW finding)

The Boundary agent identified that dynamically extending the domain during migration prevents proper termination. The fix: pre-compute the nest-B domain with correct finite extent, fill with all-1 edges, then DON'T extend during migration.

### 3. Proposition 3.1 constrains reachable orientations (NEW finding)

The Exegete noted that during correct migration, a rhombus only uses 2 of 3 orientations, ruling out the nest-C orientation. The NW ambiguity is an artifact of considering unreachable states.

### 4. Figure 2 Tracer reached a WRONG conclusion (star isolation effect)

The Tracer, isolated from other peripherals, concluded that ALL rotations are uniformly CCW. This contradicts the "avoid NW" rule (CCW from NE spokes gives NW). The vanilla center caught and corrected this error in its synthesis — but only because it read all peripherals. Had the center not been available, this error would have persisted.

### 5. Nest-B domain extents are computable (NEW finding)

The Coords agent provided a formula: nest B extends to `j = N + max(beta_i) - 1` where beta is the complement of mu.

## Convergence Analysis

| Conclusion | Ring (6/6) | Star (6/6) |
|---|---|---|
| `_MIGRATION_SHIFT` unused → use it | 6/6 | 4/6 |
| `direction_score` is a bug | 6/6 | 5/6 |
| NW is the forbidden orientation | 5/6 | 4/6 |
| NW → H via CCW | 5/6 | 3/6 |
| Domain/boundary is an issue | 6/6 | 5/6 |
| Direction rule is NOT the only bug | 1/6 (Agent 4) | 3/6 (Orient, Boundary, Vanilla) |

The star produced LESS convergence on the primary diagnosis (4-5/6 vs 6/6) but MORE insight into secondary issues (domain setup, termination).

## Diversity Comparison

### Headline Numbers

| Metric | Ring | Star | Delta |
|--------|------|------|-------|
| **Mean pairwise distance** | **0.318** | **0.450** | **+0.132** |
| Std deviation | 0.083 | 0.115 | |
| Min (most similar pair) | 0.183 | 0.237 | |
| Max (most diverse pair) | 0.455 | 0.741 | |

**GECCO prediction: ring > star. Result: REVERSED (star > ring by 41%).**

### Why the Reversal?

The GECCO paper measures raw semantic diversity — how different the texts are from each other. For this bug-diagnosis task:

1. **Ring diversity was LOW because agents converged on the CORRECT answer.** The ring's sequential information flow (each agent reads predecessor) created a convergence cascade. By Agent 5, four prior agents' findings were incorporated. The result: high agreement, low diversity.

2. **Star diversity was HIGH because isolated peripherals reached DIVERGENT conclusions.** The Tracer concluded "uniform CCW" (wrong). The Orient concluded "direction rule is fine, bug is elsewhere" (partially right but divergent from group). Without peer correction, these conclusions persisted as outliers.

3. **The vanilla center partially homogenized** — it's closest to Exegete (0.237) and Orient (0.318), absorbing their framing. But it couldn't fully correct the Boundary agent (0.468) or Coords (0.493), which developed independent analyses.

### The Most Diverse Pair: Orient vs Boundary (0.741)

This pair had a cosine distance of 0.741 — nearly orthogonal! In the ring, the same role-pair (Agent 3 vs Agent 4) had distance 0.386. The star's isolation doubled their divergence. The Orient agent focused purely on the CW/CCW algebra; the Boundary agent focused on domain geometry and termination. Without peer visibility, they developed completely different vocabularies and framings.

### Vanilla Center Effect

| Peripheral | Distance to Vanilla |
|---|---|
| Exegete | 0.237 (closest) |
| Orient | 0.318 |
| Tracer | 0.333 |
| Boundary | 0.468 |
| Coords | 0.493 (furthest) |

The center absorbed the Exegete's framing most strongly. The Coords and Boundary agents stayed furthest from the center — their technical vocabulary (grid coordinates, domain extents) was too specialized for the generalist to fully synthesize.

## Cross-Experiment Comparison

| Experiment | Topology | Task | Diversity | Quality |
|---|---|---|---|---|
| Ring debate (R3) | 4-ring | Open question | 0.437 | High diversity, speculative |
| Star debate (R3) | 5-star | Open question | 0.362 | Converged, fewer ideas |
| **Ring analysis** | **6-ring** | **Bug diagnosis** | **0.318** | **Correct diagnosis, high convergence** |
| **Star analysis** | **6-star** | **Bug diagnosis** | **0.450** | **More insights, less convergence, 1 error** |

### The GECCO Ordering and Task Type

The GECCO paper's ordering (ring > star for diversity) holds for **open-ended exploration** (debate experiments: 0.437 > 0.362). It is REVERSED for **diagnosis tasks** (analysis experiments: 0.318 < 0.450).

The mechanism:
- **Exploration tasks**: Ring's sequential propagation preserves local disagreements. Star's center homogenizes too aggressively, killing productive diversity.
- **Diagnosis tasks**: Ring's sequential propagation creates a convergence cascade — each agent reads and agrees with predecessors, compressing diversity. Star's isolation preserves peripheral independence, maintaining diversity but also preserving errors.

### Implications

1. **For bug diagnosis**: Ring topology produces better convergence on the correct answer. Star topology produces more independent analyses but includes uncorrected errors.

2. **For mathematical research**: Use ring for convergent synthesis (what IS the answer?), star for divergent exploration (what COULD the answer be?). This aligns with the debate experiment's recommendation.

3. **The GECCO ordering is task-dependent**: Ring > star for diversity in exploration. Star > ring for diversity in convergence tasks. The topology effect on diversity is real but its DIRECTION depends on whether the task rewards agreement or disagreement.
