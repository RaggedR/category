# Debate Findings: Tested and Untested Ideas

## Summary

Two multi-agent debate experiments (ring topology, star topology) produced a total of ~15 distinct mathematical ideas about the open question: does the KTW puzzle transfer matrix define a growth diagram whose face labels biject puzzles with LR tableaux?

Of these, 4 were rigorously tested (from the ring debate). Most were refuted, but the refutations revealed unexpected structure. The star debate produced 5 additional ideas, several overlapping with already-refuted ring hypotheses (confirming that the star converges toward previously-explored territory).

---

## Part 1: Ring (Cycle) Debate — Tested Results

### Experimental Setup
- **Topology**: 4-node ring (Growth Diagram → Yang-Baxter → Bijection → Computational)
- **Model**: Claude Opus, text-only, negative coupling, 3 rounds
- **Test case**: c^(3,2,1)_{(2,1),(2,1)} = 2 (12 free ∇, 3 tableau entries)
- **Second case**: c^(4,4,3,1)_{(3,2,1),(3,2,1)} = 3 (22 free ∇, 6 tableau entries)

### H1: Active/Passive ∇ Decomposition

**Source**: Computational Agent (Round 2)

**Claim**: Not all ∇ tiles contribute to the bijection. "Active" ∇ tiles affect subsequent forced △ tiles; "passive" ∇ tiles sit in frozen regions where both choices lead to identical downstream configurations. Number of active ∇ = |ν/λ|.

**Prediction**: 3 active out of 12 total (c=2 case); 6 active out of 22 (c=3 case).

**VERDICT: REFUTED**

ALL free ∇ tiles are locally active. The branching table guarantees this: `0 → {(0,0), (10,1)}` and `1 → {(0,10), (1,1)}` — the two options always produce different top values. There is no "passive" free ∇ in the local sense.

**What we learned**: The concept of "active" needs strengthening. A ∇ is locally active (trivially) but may not be *globally* active — flipping it might force a cascade of downstream changes that ultimately produces the same puzzle. However, testing global activity requires solving a path-dependent constraint propagation problem.

---

### H2: Top Faces of the Growth Diagram

**Source**: Growth Diagram Agent (Round 3)

**Claim**: The active ∇ tiles are precisely those whose position in the triangular grid borders the skew shape ν/λ. These are the "top faces" of the growth diagram. There is a geometric criterion — expressible in terms of (row, col) coordinates — that identifies which ∇ positions are active.

**Prediction**: A fixed set of ∇ positions (determined by λ, μ, ν alone) encodes the tableau.

**VERDICT: REFUTED**

The set of free ∇ positions is **not the same across puzzles within a triple**. Different valid puzzles have different sets of 12 (or 22) free positions. This demolishes the "top faces" claim: there is no single fixed set of geometrically-determined ∇ tiles that encode the bijection.

**What we learned**: The **count** of free ∇ tiles is constant across all puzzles in a triple (12 for c=2, 22 for c=3). This conservation law is non-trivial and unexplained by the "top faces" framing. Any theory must account for path-dependent free positions with a conserved total.

---

### H3: Gauge Symmetry / Bulk Invariance

**Source**: Yang-Baxter Agent (Round 3)

**Claim**: Bulk ∇ choices (those not bordering the skew shape) are "pure gauge" — changing them either breaks validity or leaves the tableau unchanged. The 3 "boundary" ∇ choices carry all observable information.

**Prediction**: There exist positions where flipping the choice bit either breaks the puzzle or preserves the tableau.

**VERDICT: CONFIRMED (structurally, with unexpected form)**

The choice space is extremely sparse:

| Case | Free ∇ (F) | 2^F | Valid puzzles (c) | Valid fraction |
|------|-----------|-----|-------------------|---------------|
| c=2  | 12 | 4,096 | 2 | 0.049% |
| c=3  | 22 | 4,194,304 | 3 | 0.00007% |

**Key discovery — Frozen Positions**: Some choice-vector positions are identical across ALL valid puzzles:

| Case | Free ∇ | Frozen positions | Frozen values | Variable positions |
|------|--------|-----------------|---------------|-------------------|
| c=2 | 12 | 6 (at indices 3,4,6,8,9,11) | `(0,1,0,1,0,1)` | 6 |
| c=3 | 22 | 7 (at indices 1,6,7,9,18,19,21) | `(0,1,0,0,1,0,1)` | 15 |

This is not the "gauge symmetry" predicted (which would be continuous freedom in bulk), but rather **rigid constraints**. The frozen positions are forced by boundary conditions via unit propagation from △ determinism.

---

### H4: Error-Correcting Code (12,3)

**Source**: Computational Agent (Round 3)

**Claim**: The free ∇ choices implement a (12,3) binary code where 3 information bits encode the LR tableau and 9 parity bits ensure valid puzzle structure.

**Predictions**:
1. Exactly c of the 2^12 vectors are valid (matches LR coefficient) ✓
2. Minimum Hamming distance between codewords = partition lattice distance ✗
3. "Syndrome decoding" reveals which 3 positions are information bits ✗

**VERDICT: REFUTED as linear code**

The code is NOT (F, |ν/λ|). The Hamming distances:

| Pair | Distance | Distinguishing positions |
|------|----------|------------------------|
| c=2: P1 vs P2 | 6 | {0,1,2,5,7,10} |
| c=3: P1 vs P2 | 6 | {10,11,12,15,17,20} |
| c=3: P1 vs P3 | 14 | {0,2,3,4,5,8,11,12,13,14,15,16,17,20} |
| c=3: P2 vs P3 | 10 | {0,2,3,4,5,8,10,13,14,16} |

The number of distinguishing positions (6 for c=2, 15 for c=3) exceeds |ν/λ| (3 and 6 respectively). The code has more "information bits" than a (F, |ν/λ|) code allows.

**What we learned**: The path-dependence of free positions (from H2) explains this — since different puzzles have different free ∇ positions, the choice vectors live in overlapping but distinct subspaces. A linear code over a fixed set of positions cannot capture this. The structure is a **nonlinear constraint satisfaction problem**, not a linear code.

---

### Unexpected Discoveries (from testing)

1. **Path-dependent free positions**: Different puzzles have different free ∇ sets, but the count F is constant. The choice vectors live in overlapping subspaces of {0,1}^(total ∇ tiles).

2. **Frozen skeleton**: 6 of 12 positions (c=2) and 7 of 22 (c=3) are rigidly forced by (λ,μ,ν). These are structural overhead — determined by boundary conditions alone.

3. **Constraint graph perspective**: Puzzles are satisfying assignments of a CSP. The ∇ choices are variables, △ determinism rules are constraints, boundary conditions are fixed. Frozen positions = unit propagation from boundaries.

4. **Puzzles differ at contiguous middle-band levels** (from tool-augmented experiment): For c=2, the two puzzles have identical level tuples everywhere except at exactly 3 levels (matching |ν/λ|). For c=3, 5 levels differ (not 6). The differing levels form a contiguous band in the triangle interior.

---

## Part 2: Star Debate — Untested Ideas

### Experimental Setup
- **Topology**: 5-node star (Vanilla center ↔ 4 specialist peripherals)
- **Model**: Claude Opus, text-only, negative coupling, 3 rounds
- **Same test cases as ring**

### S1: Migration Zone = Partial Growth Diagram

**Source**: Bijection Agent (Rounds 2-3)

**Claim**: Only ∇ triangles in positions corresponding to boxes of ν/λ participate in Purbhoo's migration. The growth diagram exists but only on a **subgraph** of the full puzzle — specifically, the ∇ tiles that have mosaic squares which actually move during migration. The other ∇s contribute to puzzle multiplicity but not tableau content.

**Prediction**: For c^(3,2,1)_{(2,1),(2,1)}, exactly 3 of the 12 ∇ tiles correspond to migrating mosaic squares.

**Status**: UNTESTED

**Relationship to H1**: Overlaps with H1 (active/passive decomposition), but the criterion is different. H1 defined "active" via downstream influence on △ tiles (refuted — all are locally active). S1 defines "relevant" via participation in Purbhoo's migration (a global, non-local property). These are distinct claims.

**Relationship to H2**: H2's refutation (free positions are puzzle-dependent) may not apply here, because S1 asks about *migrating squares* in the mosaic representation, not about free/forced status in the enumeration tree.

**Implementation path**: Requires implementing Purbhoo's migration algorithm on the mosaic representation. Infrastructure exists (`~/git/puzzles/mosaic.py`).

---

### S2: Transfer Matrix Kernel at x=1

**Source**: Vanilla Agent (Round 2)

**Claim**: The 12→3 projection is the kernel of the transfer matrix at the free-fermionic specialisation x=1. At generic spectral parameter x, the transfer matrix has full rank. At x=1, the kernel dimension equals the redundant bits: dim(ker) = F − |ν/λ| = 12 − 3 = 9.

**Prediction**: Compute the transfer matrix symbolically; its kernel at x=1 has dimension 9 (c=2) or 16 (c=3).

**Status**: UNTESTED

**Assessment**: Conceptually interesting but mathematically vague. The "transfer matrix" in Zinn-Justin's framework is an operator on a Hilbert space (tensor product of local state spaces), not a 12×12 matrix acting on choice vectors. The claim would need to be made precise: which matrix exactly, acting on which space?

---

### S3: Yang-Baxter Orbits on Choice Vectors

**Source**: Computational Agent (Round 2)

**Claim**: The Yang-Baxter equation defines an equivalence relation on the 4096 choice vectors (for c=2). Apply R-matrix moves as local transformations on adjacent ∇ choices. The number of orbits under this action should equal c (the LR coefficient). Each orbit corresponds to one LR tableau.

**Prediction**: Generate all 2^12 vectors, apply YBE moves, count orbits. Should get exactly 2 orbits for c=2.

**Status**: UNTESTED

**Assessment**: The most theoretically exciting idea. If true, it provides a *new proof* of the LR rule via quantum group symmetry. The key question: how exactly does R act on adjacent ∇ choices? Zinn-Justin's R-matrix acts on pairs of row states, not directly on ∇ choice bits. Translating between these representations is non-trivial but potentially doable.

**Concern**: The 4096 vectors include invalid puzzles. The YBE action preserves validity (it permutes solutions of the YBE). So the 2 valid puzzles should form 2 singleton orbits among the valid vectors — which is trivially true and uninteresting. The interesting version would be: YBE acts on some larger set and the orbits landing in "valid puzzle space" correspond to tableaux.

---

### S4: Spectral Parameter Criticality

**Source**: Yang-Baxter Agent (Rounds 2-3)

**Claim**: ∇ triangles contribute to the tableau precisely when their spectral parameters satisfy a criticality condition (u−v = 0). The 3 "relevant" ∇s aren't at fixed positions — they're at positions where the spectral parameters of adjacent transfer matrices coincide.

**Prediction**: For each ∇ at position (i,j), compute spectral parameters u_i, v_j from the boundary data. The ∇ contributes to the tableau iff u_i = v_j.

**Status**: UNTESTED

**Assessment**: Requires understanding how spectral parameters are assigned in the free-fermionic puzzle model. In Zinn-Justin's formulation, spectral parameters are attached to *rows and columns* of the lattice, not individual vertices. The condition u_i = v_j has a clear meaning (resonance point), but whether it selects exactly |ν/λ| positions needs verification.

---

### S5: Non-crossing Migration Paths as Poset Structure

**Source**: Bijection Agent (Round 3)

**Claim**: Migration paths in Purbhoo's algorithm are non-crossing (Lemma 8 of arXiv:0705.1184). This non-crossing constraint defines a partial order on the migrating tiles. This partial order is the poset structure that a growth diagram needs.

**Prediction**: The non-crossing paths define a lattice structure on the |ν/λ| migrating squares, and the growth diagram's local rule operates on this lattice (not on the full triangular grid).

**Status**: UNTESTED

**Assessment**: This is the most structurally coherent proposal for *how* a partial growth diagram could work. Non-crossing paths → Tamari lattice or similar → natural poset for growth diagram faces. Would need to verify that the local rule (Burge-type) is well-defined on this restricted poset.

---

## Part 3: Ideas Shared Between Both Debates (Convergent)

These ideas appeared in both the ring and star debates, suggesting they are natural attractors for this problem:

1. **12:3 ratio is a projection/compression** — appeared in all agents across both topologies. The precise mechanism is disputed (linear code, kernel, migration zone, gauge quotient) but the framing is universal.

2. **Cylindric growth diagrams (Langer §4.3)** — Growth Diagram Agent proposed this in both debates. The claim: each puzzle row is one application of L_i, and the ∇ choices are Burge parameters. Neither confirmed nor refuted — remains the leading structural framework.

3. **Yang-Baxter ⇒ growth diagram coherence** — appeared in both. The claim that YBE ⇒ diamond lemma ⇒ well-defined growth diagram is correct as a *necessary* condition but not proven sufficient for KTW puzzles specifically.

4. **Purbhoo's migration is fundamentally non-local** — Bijection Agent stressed in both debates that jdt cascades across the full tableau. Any "local" growth diagram interpretation must explain how non-local migration emerges from local face labels.

---

## Part 4: What to Implement Next

### Priority 1: Purbhoo's Migration Algorithm

**Why**: It's the only existing bijection (puzzle → LR tableau). Understanding it concretely — by running it — will reveal:
- Which mosaic squares actually migrate (testing S1)
- Whether migration paths are non-crossing (testing S5)
- The exact non-locality structure (informing all other hypotheses)

**Infrastructure**: `~/git/puzzles/mosaic.py` and `~/git/puzzles/hex_mosaic.py` exist with edge-label storage, face enumeration, nest geometry, flock construction, and hexagonal rotations.

**Current status**: Partially implemented. Steps 1-4 work correctly:
- `HexMosaic.from_puzzle()` builds the hexagonal mosaic with nests
- `build_flock_left()` labels the flock in standard order
- `rotate_180()` and `rotate_120()` correctly permute spoke labels (10-count conserved when domain is complete)
- `migrate_single()` correctly migrates the FIRST rhombus through the interior

**Blocking issue**: The hexagonal domain boundary is incompletely defined. After the first rhombus migrates through the interior into nest B, it encounters the PRE-PLACED nest-B rhombi and keeps swapping with them, overshooting its target position. The termination condition (when to stop) is wrong.

**Root cause**: Purbhoo's algorithm operates on a hexagonal domain where nests have specific tile structures (rhombi packed in partition shapes, separated by gap triangles). The current code either:
1. Leaves nest edges undefined (None) → rotations break at boundary
2. Fills them with default 1-edges → creates fake rhombi that cause overshooting

**What's needed to fix**:
1. **Remove nest-B rhombi before migration** (nest B starts empty in Purbhoo's setup — the migrating flock CREATES the nest-B content). Keep nest-A rhombi (the flock source) and nest-C rhombi.
2. **Properly define the nest-B domain**: the empty nest B should have all edges = 1 (all-1 triangles) in the region where migrating rhombi will settle. No pre-placed 10-edges.
3. **Fix termination**: a rhombus is "settled" when it's in nest B AND the pivot vertex's opposite spoke has a 0-label (boundary of the packed region), not just when it crosses the boundary.

Alternative approach: **bypass hexagonal migration entirely** and implement the bijection via the BZ-pattern / strand-diagram route, which operates on the combinatorial puzzle structure without geometric domain issues.

### Priority 2: Yang-Baxter Orbit Computation (S3)

**Why**: If orbits under R-matrix action correspond to LR tableaux, this is a new result connecting quantum groups to the LR rule — potentially publishable independently.

**What's needed**:
1. Read Zinn-Justin equation (7) carefully to define R-matrix on local state pairs
2. Translate R-action from row states to ∇ choice bits
3. Generate orbit structure on the valid choice vectors
4. Check if orbit count = c

**Risk**: High — the translation step (2) may reveal that R doesn't act cleanly on choice vectors at all.

### Priority 3: Spectral Parameter Assignment (S4)

**Why**: Quick to test if we understand the parameter assignment. Would immediately confirm or refute whether "criticality" selects |ν/λ| positions.

**What's needed**: Read Zinn-Justin's paper to determine how u_i and v_j are assigned from (λ,μ,ν). Check if |{(i,j) : u_i = v_j}| = |ν/λ|.

---

## Part 5: Lessons for Multi-Agent Debate

### What the ring produced that the star didn't:
- More diverse and **persistent** disagreements (diversity stayed above 0.43)
- Ideas that were wrong in interesting ways (H4's error-correcting code framing led to discovering the frozen skeleton)
- Genuine conceptual tension between the bijection and algebraic perspectives

### What the star produced that the ring didn't:
- A clear synthesis (Vanilla Agent's "migration zone = partial growth diagram")
- Faster convergence to implementable proposals
- Less hallucinated computational detail (fewer unfounded claims about specific bit vectors)

### Topology recommendation for mathematical research:
- **Ring** for divergent exploration (generating hypotheses)
- **Star** for convergent synthesis (finding the best path forward)
- **Sequence**: Ring first (explore), then Star (synthesize), then implement
