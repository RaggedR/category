# Testable Hypotheses

Extracted from 3-round Opus ring-topology debate (2026-05-04, negative coupling).
All hypotheses concern the open question: does the KTW puzzle transfer matrix define a growth diagram whose face labels biject puzzles with LR tableaux?

Minimal test case: c^(3,2,1)_{(2,1),(2,1)} = 2, with 12 free ∇ choices and 3 tableau entries.
Second test case: c^(4,4,3,1)_{(3,2,1),(3,2,1)} = 3, with 22 free ∇ choices and 6 tableau entries.

## H1: Active/Passive Decomposition (Computational Agent, Round 2)

**Claim**: Not all ∇ tiles contribute to the bijection. "Active" ∇ tiles are those whose branching choice affects subsequent forced △ tiles. "Passive" ∇ tiles sit in frozen regions where both choices lead to identical downstream configurations.

**Prediction**: The number of active ∇ tiles equals |ν/λ| (the number of boxes in the skew shape). For the minimal case: 3 active out of 12 total. For the second case: 6 active out of 22 total.

**Status**: REFUTED (2026-05-04)

**Result**: ALL free ∇ tiles are locally active. The branching table guarantees this: `0 → {(0,0), (10,1)}` and `1 → {(0,10), (1,1)}` — the two options always produce different top─ values. There is no "passive" free ∇ in the local sense.

The H1 definition of "active" was too weak. A stronger definition would be: a ∇ is "globally active" if flipping its choice (while propagating downstream to maintain validity) changes the final puzzle. This remains untested because the path-dependence of free positions (see H2 result) complicates perturbation analysis — flipping one ∇ may change which subsequent ∇ tiles are even free.

## H2: Top Faces of the Growth Diagram (Growth Diagram Agent, Round 3)

**Claim**: The active ∇ tiles are precisely those whose position in the triangular grid borders the skew shape ν/λ. These are the "top faces" of the growth diagram when restricted to the skew region. The remaining ∇ tiles are in the "bulk," contributing to intermediate partition evolution but not to the final recording tableau.

**Prediction**: There is a geometric criterion (expressible in terms of (row, col) coordinates in the puzzle triangle and the partitions λ, μ, ν) that identifies which ∇ positions are active, without needing perturbation analysis.

**Status**: PARTIALLY REFUTED (2026-05-04)

**Result**: The set of free ∇ positions is **not the same across puzzles** within a triple. Different valid puzzles have different sets of 12 (or 22) free positions. This demolishes the "top faces" claim as stated — there is no single fixed set of geometrically-determined ∇ tiles that encode the bijection. The identity of which ∇ tiles are free depends on the upstream choices, not just the boundary partitions.

However, the **count** of free ∇ tiles is constant across all puzzles in a triple (12 for c=2, 22 for c=3). This is a non-trivial conservation law that the "top faces" framing failed to explain.

## H3: Bulk Invariance / Gauge Symmetry (Yang-Baxter Agent, Round 3)

**Claim**: Changing "bulk" ∇ choices (those not bordering the skew shape) while preserving "boundary" choices always produces an invalid puzzle — OR produces a valid puzzle with the same LR tableau. In gauge theory language: bulk degrees of freedom are pure gauge and don't affect the observable (the tableau).

**Prediction**: For c^(3,2,1)_{(2,1),(2,1)} = 2, there exist ∇ positions where flipping the choice bit either (a) breaks puzzle validity or (b) leaves the corresponding LR tableau unchanged.

**Status**: CONFIRMED (structurally) (2026-05-04)

**Result**: The choice space is extremely sparse:

| Case | Free ∇ (F) | 2^F | Valid (c) | Ratio |
|------|-----------|-----|-----------|-------|
| c=2  | 12 | 4,096 | 2 | 0.049% |
| c=3  | 22 | 4,194,304 | 3 | 0.00007% |

Every puzzle has the same number of free ∇ choices (F is constant within a triple). The enumeration tree is path-dependent (free positions shift between puzzles), so "2^F" overstates the actual reachable space. Nevertheless, the extreme sparsity confirms that puzzle validity imposes massive constraints on the choice vectors.

**Key discovery — Frozen Positions**: Some choice-vector positions are **identical across all valid puzzles** within a triple:

| Case | Free ∇ | Frozen positions | Frozen values | Variable positions |
|------|--------|-----------------|---------------|-------------------|
| c=2 | 12 | 6 (positions 3,4,6,8,9,11) | `(0,1,0,1,0,1)` | 6 |
| c=3 | 22 | 7 (positions 1,6,7,9,18,19,21) | `(0,1,0,0,1,0,1)` | 15 |

This is a form of bulk invariance, but not the one the Yang-Baxter Agent predicted. The frozen positions aren't "gauge" — they're **rigidly determined by the boundary conditions** (λ,μ,ν). The variable positions carry all the puzzle-distinguishing information.

## H4: Error-Correcting Code Structure (Computational Agent, Round 3)

**Claim**: The free ∇ choices implement a (12, 3) binary code where 3 information bits encode the LR tableau and 9 parity bits ensure valid puzzle structure. Valid puzzles = codewords satisfying all parity checks from △ determinism.

**Prediction**:
1. Exactly c of the 2^12 vectors are valid (matches the LR coefficient)
2. The minimum Hamming distance between valid codewords equals the partition lattice distance between corresponding LR tableaux
3. "Syndrome decoding" of invalid vectors reveals which 3 positions are information bits

**Status**: REFUTED as stated, but reveals deeper structure (2026-05-04)

**Result**: The code is NOT (F, |ν/λ|). The number of distinguishing positions (where any two codewords differ) exceeds |ν/λ|:

| Case | F | |ν/λ| | Distinguishing positions | Frozen positions |
|------|---|------|-------------------------|-----------------|
| c=2 | 12 | 3 | 6 | 6 |
| c=3 | 22 | 6 | 15 | 7 |

The Hamming distances between valid codewords:

**c=2 case**: Puzzle 1 vs 2: distance 6, differing at positions {0,1,2,5,7,10}.

**c=3 case**:
- Puzzle 1 vs 2: distance 6, positions {10,11,12,15,17,20}
- Puzzle 1 vs 3: distance 14, positions {0,2,3,4,5,8,11,12,13,14,15,16,17,20}
- Puzzle 2 vs 3: distance 10, positions {0,2,3,4,5,8,10,13,14,16}

The code has more "information bits" than |ν/λ|. The path-dependence of free positions (H2 result) explains this: since different puzzles have different free ∇ positions, the choice vectors live in different (overlapping) subspaces. A simple linear code over a fixed set of positions cannot capture this structure.

**What IS true**: The frozen positions form a rigid "skeleton" — identical values across all valid puzzles. This is consistent with a **nonlinear code** or a **constraint satisfaction problem** where the frozen variables are forced by unit propagation from the boundary conditions.

## Status

- [x] H1: REFUTED — all free ∇ are locally active (the definition was too weak)
- [x] H2: PARTIALLY REFUTED — free positions are puzzle-dependent, not geometrically fixed
- [x] H3: CONFIRMED — extreme sparsity + frozen positions discovered
- [x] H4: REFUTED as (F, |ν/λ|) code — but nonlinear constraint structure exists

## Unexpected Discoveries

### From Experiment 1 (text-only debate, `test_hypotheses.py`)

#### 1. Path-dependent free positions
Different valid puzzles have different sets of free ∇ positions, even though the count F is constant. This means the choice vectors for different puzzles live in overlapping but distinct subspaces of {0,1}^(total ∇ tiles). Any theory of the bijection must account for this: the "recording data" is not written to a fixed set of face positions.

#### 2. Frozen skeleton
Some choice-vector positions are rigidly determined by (λ,μ,ν) — identical across all valid puzzles. These frozen values likely propagate from the boundary conditions via △ determinism (unit propagation in the constraint graph). The frozen positions could be the "structural overhead" mentioned in OPEN_QUESTION.md.

| Case | Free ∇ | Frozen positions | Frozen values | Variable positions |
|------|--------|-----------------|---------------|-------------------|
| c=2 | 12 | 6 (positions 3,4,6,8,9,11) | `(0,1,0,1,0,1)` | 6 |
| c=3 | 22 | 7 (positions 1,6,7,9,18,19,21) | `(0,1,0,0,1,0,1)` | 15 |

#### 3. The constraint graph perspective
The results suggest viewing puzzles as solutions to a constraint satisfaction problem rather than codewords in a linear code. The ∇ branching choices are variables, the △ determinism rules are constraints, and the boundary conditions (λ,μ,ν) are fixed. Valid puzzles = satisfying assignments. The frozen positions are variables forced by constraint propagation. This connects to the statistical physics of random CSPs (phase transitions, freezing thresholds) more naturally than to coding theory.

### From Experiment 2 (tool-augmented debate, `sonnet_tools_v2.json`)

#### 4. Puzzles differ at exactly |ν/λ| levels (VERIFIED)

The Computational Agent discovered — and we independently verified — that for c^(3,2,1)_{(2,1),(2,1)} = 2, the two puzzles have identical `levels` tuples everywhere except at exactly 3 levels, matching |ν/λ| = |{(3,2,1)}| - |{(2,1)}| = 3.

```
Level 0: (1, 0, 1, 0, 1, 0) == (1, 0, 1, 0, 1, 0)   ← ν boundary (identical)
Level 1: (1, 0, 1, 10, 0)   != (1, 10, 0, 1, 0)      ← 3 element diffs
Level 2: (0, 1, 1, 0)       != (1, 0, 0, 1)           ← 4 element diffs
Level 3: (10, 1, 0)         != (1, 0, 10)              ← 3 element diffs
Level 4: (0, 1)             == (0, 1)                   ← identical
Level 5: (10,)              == (10,)                    ← identical
Level 6: ()                 == ()                       ← apex
```

The differing levels (1, 2, 3) form the "middle band" of the triangle — the region between the ν boundary (bottom) and the λ/μ convergence zone (top). This middle band is exactly where the skew shape ν/λ lives.

**Significance**: This is a different decomposition from the choice-vector analysis (Experiment 1). It operates on the *output representation* (puzzle levels), not the *enumeration process* (∇ choices). The fact that both decompositions yield |ν/λ| — one as the number of differing levels, the other as the predicted number of "information bits" — suggests a deep structural connection.

**Verified on c=3 case** (c^(4,4,3,1)_{(3,2,1),(3,2,1)}, |ν/λ| = 6, N=8 so 9 levels):

```
Levels identical across all 3 puzzles: [0, 6, 7, 8]  (count=4)
Levels that differ (in at least one pair): [1, 2, 3, 4, 5]  (count=5)

Pairwise:
  Puzzles 1 vs 2: 3 differing levels at [3, 4, 5]
  Puzzles 1 vs 3: 5 differing levels at [1, 2, 3, 4, 5]
  Puzzles 2 vs 3: 3 differing levels at [1, 2, 3]
```

**Result**: 5 differing levels, not 6. The |ν/λ|-levels conjecture does NOT hold exactly for c=3. However, the "middle band" pattern persists: levels 0 (ν boundary) and 6-8 (top) are always identical. The differing levels form a contiguous band [1..5] in the interior.

The c=2 match (3 differing levels = |ν/λ| = 3) may have been coincidental, or the relationship is more nuanced (e.g., the number of differing levels equals N - 1 - (length of longest part of λ), or similar).

#### 5. Debate self-correction via tool use

In Round 1 of Experiment 2, the Bijection Agent claimed "no branching choices exist in the puzzle representation" — confusing the final puzzle object with the enumeration process. The Computational Agent caught this error in the same round: "The branching IS present — it's in the transfer matrix, not the static output!" By Round 2, the debate had converged on the correct interpretation: branching occurs during `_row_transitions`, but the final puzzles are determinate configurations that differ at specific levels.

This self-correction did not happen in Experiment 1 (text-only), where agents fabricated computational claims without verification. Tool use made the debate empirically grounded and self-correcting.

## Experiment Comparison

| | Experiment 1 (text-only) | Experiment 2 (tools) |
|---|---|---|
| Model | Opus | Sonnet |
| Rounds | 3 | 2 |
| Diversity (mean) | 0.43 → 0.53 → 0.44 | 0.35 → 0.30 |
| Hallucinated data | Yes — fabricated choice vectors | No — all claims verified by code |
| Key finding | Frozen skeleton + path-dependent free positions | Puzzles differ at exactly |ν/λ| levels |
| Self-correction | No | Yes (Bijection Agent corrected in Round 1) |
| Cost | ~$0.45 | ~$0.80 |

**Diversity was lower with tools** (0.35 vs 0.43). Agents converge when they can see the same data. Text-only agents diverge more freely but make unfounded claims.

## Bijection Attempt (2026-05-04, manual investigation)

After the two debate experiments, we attempted to construct the actual bijection: puzzle → LR tableau. This is the core question in OPEN_QUESTION.md.

### Setup

Minimal test case: c^(3,2,1)_{(2,1),(2,1)} = 2. Skew shape ν/λ has 3 boxes: (0,2), (1,1), (2,0). Content μ = (2,1) means entries are two 1's and one 2.

The two LR tableaux:
```
T1: . . 1    T2: . . 1
    . 1          . 2
    2            1
```

### What we found

**1. Intermediate levels with 10 are NOT partition boundaries.** The puzzle levels at odd positions contain the label 10 (rhombus), which is not in {0,1}. Applying `boundary_to_partition` to these levels gives nonsensical results. Only even levels are pure {0,1} and decode to valid partitions.

**2. Each R (=10) in a level "carries" exactly one box.** Expanding 10 → (1,0) gives an "upper" partition; expanding 10 → (0,1) gives a "lower" partition. The skew between them is a single box — the box the R is transporting through the triangle.

```
Puzzle 1:                         Puzzle 2:
  Level 1: R carries box (2,1)     Level 1: R carries box (1,2)
  Level 3: R carries box (0,1)     Level 3: R carries box (1,0)
  Level 5: R carries box (0,0)     Level 5: R carries box (0,0)
```

**3. R appears only at odd levels — exactly where the left boundary has a 1-bit.** For λ = (2,1), the boundary is (0,1,0,1,0,1). The 1-bits at positions 1,3,5 correspond to "row boundaries" of λ. Each R at level 2k+1 corresponds to processing row k of both λ and μ simultaneously.

**4. The partition chain at even levels differs between puzzles.**

```
Puzzle 1: (3,2,1) → (1,1) → () → ()
Puzzle 2: (3,2,1) → (2,)  → () → ()
```

The intermediate partition at level 2 is (1,1) vs (2,) — both have 2 boxes but different shapes. This determines how the skew is being decomposed.

### What we could NOT find

**The map from R-carried boxes to tableau entries.** We tried:

| Attempt | P1 → entries | P2 → entries | Correct? |
|---------|-------------|-------------|----------|
| Carried box row | [2,0,0] | [1,1,0] | No — T1=[1,1,2], T2=[1,2,1] |
| n − carried row | [1,3,3] | [2,2,3] | No |
| R column in level | [3,0,0] | [1,2,0] | No |
| R column from right | [1,2,0] | [3,0,0] | No |

No simple function of the R position or carried-box coordinates gives the tableau entry. The bijection appears to require tracking the full interaction pattern between left (λ) and right (μ) boundary paths as they propagate through the triangle, not just the R positions in isolation.

### Why it's hard

The core difficulty: the bijection between puzzles and LR tableaux classically goes through **Purbhoo's migration** (puzzle → mosaic → jeu de taquin → tableau), which is inherently non-local. A single jdt slide can cascade across the entire tableau. The OPEN_QUESTION asks whether the growth diagram framework can provide a *local* alternative — reading ∇ face labels directly — but our investigation suggests the encoding is more entangled than a per-face-label readoff.

### What a solution would look like

Based on what we've found, a successful bijection would need to:

1. **Identify a reading order** on the ∇ faces — likely related to the migration sequence, not the left-to-right row-by-row order
2. **Project** the full choice sequence onto the |ν/λ| = 3 skew boxes — the projection must use the boundary conditions (λ,μ,ν) to determine which faces contribute
3. **Map** each projected face label to a tableau entry from {1,...,len(μ)} — this step probably uses the μ boundary (right side) to assign entry values

The R-at-odd-levels observation narrows the search: the bijection data lives in the odd-level faces, and the even-level partitions provide the structural scaffold. The reading order probably follows the R-paths through the triangle rather than the row-by-row scan.

### Partial pattern at level 3

At R-level 3 (λ-row 1 crossing), the carried box row directly gives the tableau entry:
- P1 carries (0,1) → row 0 → entry 1 for skew box (1,1) — matches T1
- P2 carries (1,0) → row 1 → entry 2 for skew box (1,1) — matches T2

But this formula (entry = carried_row + 1) fails at levels 1 and 5. The bijection appears to depend on the full migration context, not just individual R-carried boxes in isolation. The non-locality of jeu de taquin manifests here: the entry at one skew box depends on what happened at other skew boxes during earlier steps of the migration.

### Status

**Open.** The session produced structural observations that narrow the search:
1. R appears only at odd levels = λ row boundaries
2. Each R carries exactly one box (upper/lower partition split)
3. The carried box's row correlates with the tableau entry at some levels
4. The full bijection requires tracking the cumulative effect of all R's — it's not a per-R formula

The next step would be to implement Purbhoo's migration algorithm directly on the mosaic representation, which requires reconstructing the corner rhombi from the puzzle boundary data and implementing the local hexagonal rotations from Section 3.1.

## Mosaic Structure — Key Findings (2026-05-04/05, from reading Purbhoo arXiv:0705.1184)

### The mosaic is NOT the puzzle re-labeled — it tiles a DIFFERENT region

| | Puzzle | Mosaic |
|---|---|---|
| **Domain** | Equilateral triangle, side n | Hexagon A'AB'BC'C (150° at A,B,C; 90° at A',B',C') |
| **Tile types** | 2 triangles + 1 rhombus (the "10" piece) | Triangle + **square** + rhombus |
| **Rhombi location** | Distributed throughout interior | **Only in the 3 nests** (corners A, B, C) |
| **Interior** | Triangles + rhombi | Triangles + **squares** (no rhombi) |

### The puzzle↔mosaic bijection (Section 2.4)

Physically: take a mosaic, remove all nest rhombi, then imagine the squares have corner hinges. Pull corners A', B', C' outward — the hexagon straightens to a triangle, and **squares flex into the puzzle's 10-labeled rhombi** (30°/60° rhombi). The 0/1 edge labels come from which way the line segments rotated during straightening.

**This is why KTW puzzles are rigid under hexagonal flips** (verified computationally: 0 available flips for both c=2 puzzles). The squares have already become rigid rhombi — there's no slack for local rearrangement.

### Rhombus counts confirm the structure

| Case | Puzzle rhombi | n×m | Shared across puzzles | |ν/λ| |
|------|--------------|------|----------------------|-------|
| c=2 (N=6) | 9 | 3×3=9 | 3 | 3 |
| c=3 (N=8) | 16 | 4×4=16 | — | 6 |

Each puzzle rhombus = one square in the mosaic = one box of the n×m bounding rectangle. The 3 shared rhombi (invariant across puzzles) are at the boundary — forced by boundary conditions alone.

### Migration algorithm (Section 3.1)

**Input**: mosaic, accessible flock F = (ζ, (E,N), f), compass direction.

**Process** (for north/east migration): take the rhombus in ζ which is **maximal** in the standard order. Proceeding largest-to-smallest:

For each rhombus ◇ in the flock:
1. Look to the right of ◇; find the **smallest hexagon** containing ◇ that is entirely weakly right of ◇'s leftmost point
2. **2-fold symmetric** hexagon → rotate tiling 180°
3. **3-fold symmetric** hexagon → rotate tiling 120° (direction chosen so some edge of ◇ becomes horizontal or vertical)
4. Repeat until ◇ reaches the target nest

**Key insight**: the 2-fold rotation is a rhombus↔square swap. The 3-fold rotation cycles triangle+square+rhombus. These are the "local hexagonal rotations" — they operate on the mosaic's mixed tile types, not on the rigid KTW tiles.

### The bijection (Corollary 3.2)

```
Mosaic(α,β,γ) → give α orientation (-E_A, -N_A)
              → unique flock F on α with that orientation (Lemma 2.5)
              → migrate F from nest A to nest B (south)
              → F' at nest B is identified with an LR-tableau
```

Inverse: LR-tableau on λ/μ → form mosaic (∅, γ∨, γ) where γ∨ is identified with λ → flock G in γ∨ → migrate from B to A → recover original mosaic.

### Migration = jeu de taquin (Proposition 5.2)

Each single-rhombus migration from B to C (or vice versa) corresponds to one Schützenberger slide. This means the full migration is a sequence of jdt slides.

### Implementation path

To implement the migration computationally:
1. **Build the hexagonal mosaic** from a puzzle by reversing the straightening: add nests at corners, convert puzzle rhombi (10-edges) back to squares
2. **Place nest rhombi** at nest A according to boundary conditions
3. **Label the flock** using Lemma 2.5 (unique labeling given orientation)
4. **Implement 2-fold and 3-fold rotations** on hexagonal sub-regions
5. **Run migration** largest-to-smallest in standard order
6. **Read LR-tableau** from the flock at nest B

Alternative (from Section 5.1): since migration = jdt, could implement directly as Schützenberger slides on tableaux. But this requires translating the mosaic geometry into tableau coordinates, which is itself non-trivial.

### Code produced

- `~/git/puzzles/mosaic.py` — Mosaic class with edge-label storage, face enumeration, adjacency maps, rhombus detection. Confirmed both c=2 puzzles are valid and have 9 rhombi each (= n×m). Flip graph BFS confirms 0 reachable mosaics from either puzzle (rigidity).

## Open Questions for Next Investigation

1. ~~Does the |ν/λ|-levels result hold for c=3?~~ **Answered**: No — 5 differing levels for |ν/λ| = 6. But the "contiguous middle band" pattern holds. What determines the band width?
2. **What determines the frozen positions?** Can they be predicted from (λ,μ,ν) without enumerating puzzles?
3. **Connection between frozen choice positions and identical levels?** The 6 frozen choice positions (Experiment 1) and the 4 identical levels (Experiment 2) must be related — what is the map?
4. **How does the free-position set shift between puzzles?** The positions that are free in puzzle 1 but forced in puzzle 2 might reveal the growth diagram structure.
5. **Does the constraint graph have the structure of a growth diagram?** The DAG of constraint propagation (which ∇ forces which △ forces which ∇) might itself be the growth diagram.
6. **Can we implement the hexagonal mosaic directly?** The key challenge is representing the hexagonal domain with its nests, implementing the square↔rhombus conversion, and coding the 2-fold/3-fold rotations.
7. **Is the jdt shortcut viable?** Since migration = Schützenberger slides, can we bypass the full mosaic and compute jdt directly from the puzzle structure? This would require understanding which jdt slides correspond to which puzzle data.

## Test Scripts

- `test_hypotheses.py` — H1-H4 tests on Experiment 1 output, using `lr_coefficients.py`
- `sonnet_tools_v2.json` — Full Experiment 2 output including all tool calls and agent code
- Verification of Discovery 4 ran inline (see above)
