# Agent 4 Findings: Boundary & Termination Analysis

## 1. The Initial State of Nest B Before Migration

### What Purbhoo's Paper Tells Us

The paper constructs a mosaic with boundary data (alpha, beta, gamma) where alpha, beta, gamma are the collections of rhombi packed in nests A, B, C respectively. The mosaic is a COMPLETE tiling of the hexagonal region -- every tile is placed, every edge is labelled. The three nests are FULL from the start.

Lemma 2.6 is critical here: "If alpha = emptyset, then for any beta, there is a unique mosaic (alpha, beta, gamma). In this mosaic gamma is the complement of beta." This tells us what happens when a nest is EMPTY of rhombi -- but this is a special case. In the general case all three nests contain rhombi.

For the bijection (Corollary 3.2), the mosaic starts as a COMPLETE tiling with ALL three nests populated. Nest B contains beta (determined by mu). The migration operation then moves alpha's rhombi from nest A toward nest B. The key insight is: **nest B is NOT empty before migration. It is full of beta-rhombi.** The migrating flock must displace these existing tiles as it moves through the interior and into the nest region.

### What This Means for the Code

The `from_puzzle` constructor in `hex_mosaic.py` (line 195-234) correctly builds the full mosaic with rhombi in ALL three nests. The right nest gets H-oriented rhombi via `_place_h_rhombus`. This is correct -- the mosaic should be complete.

However, the HANDOFF's reproduction script does something strange: it REMOVES nest-B edges before migration:
```python
for ek in [ek for ek in hm.edges if any(rv >= 0 and jv > N-rv for rv,jv in ek)]:
    del hm.edges[ek]
```
This strips out the entire nest-B region. It then fills with 1-edges in the left and bottom nests. The HANDOFF even warns "Don't fill nest B with 1-edges" and "Don't fill nest B with pre-placed rhombi" -- indicating previous attempts to manipulate nest B's initial state all failed.

### My Analysis: The Fundamental Misconception

I believe the problem is a misunderstanding of what migration DOES. In Purbhoo's framework:

1. The mosaic is a complete tiling. Nest B has rhombi beta.
2. Migration does NOT require nest B to be "empty" or "cleared". Migration is a sequence of LOCAL rotations (180-degree and 120-degree) that propagate through the tiling.
3. Each rotation displaces tiles locally. The migrating rhombus pushes through the existing tiling, rearranging tiles along its path (creating the "wake").
4. The beta-rhombi in nest B get rearranged by these rotations. They don't collide with the migrating rhombus -- they get MOVED OUT OF THE WAY by the rotation operations.

**The correct initial state for nest B is: whatever beta-rhombi the mosaic naturally contains.** Don't delete them. Don't add 1-edges. The mosaic should be left as `from_puzzle` constructs it.

The reason previous approaches failed:
- Deleting nest-B edges creates an incomplete tiling. The rotation operations at the boundary then encounter missing edges and behave incorrectly.
- Filling with 1-edges creates a tiling where every edge is 1 (no rhombi), which is the wrong mosaic -- it's the unique mosaic with alpha = emptyset (Lemma 2.6).
- Pre-placing rhombi incorrectly causes symmetry-type mismatches, leading to wrong rotation choices.

## 2. The Termination Condition: What Does "Packed Into a New Nest" Mean?

### Purbhoo's Statement

Page 14: "Repeat this process, until the rhombus is packed into a new nest." Proposition 3.1: "All rhombi in zeta' are in the same nest. Consequently, zeta' is a transformed diagram."

The proof of Proposition 3.1 is elegant: "the rhombi can only ever be in two of three orientations during the migration process, thereby ruling out one of the nests." Since each nest has a characteristic orientation pair, a rhombus that has only been in two orientations can only end up in one of two nests -- and since it started at A and is moving toward B, it ends at B.

### What "Packed" Actually Means

From Section 2.4 (page 9): rhombi are "packed" in the nest if they "lie inside the cone, have edges parallel to the cone, and are a linear transformation of a straight diagram whose corner is the corner of the cone."

This means a rhombus is packed when:
1. It is geometrically inside the nest-B cone region
2. Its edges are aligned with the nest-B orientation (E_B, N_B)
3. Together with other packed rhombi, it forms a transformed Young diagram shape

### The Current `_is_in_nest_b` Function

The current code (line 833-851) checks if a rhombus is GEOMETRICALLY past the right boundary:
```python
if orient == 'nw': return j + 1 > N - r
elif orient == 'h': return j >= N - r - 1
elif orient == 'ne': return j + 1 >= N - r
```

But this is only a NECESSARY condition, not a SUFFICIENT one. Being past the boundary doesn't mean the rhombus is "packed". The current code then checks `_is_stuck` (line 854-883) -- whether the rhombus has no more valid moves.

### The Real Termination Condition

The correct termination condition should be: **the rhombus has no valid move, period.** Not "the rhombus is in nest B AND stuck," but simply "stuck." Here is why:

Purbhoo's algorithm says: repeat the rotation process until the rhombus is packed. The algorithm is DETERMINISTIC -- at each step, there is exactly one valid rotation (either the 180-degree or one of the 120-degree options). When no valid rotation exists, the rhombus has necessarily arrived at its final position.

The `_is_in_nest_b` check is actually unnecessary and potentially harmful. If the direction rule is correct, the rhombus WILL end up in nest B (by Proposition 3.1). The algorithm should simply iterate until no move is possible.

However, there is a subtlety: at intermediate positions inside the interior, a rhombus should ALWAYS have a valid move. If it doesn't, something has gone wrong (wrong rotation direction, or missing edges). So `_is_stuck` being true in the interior is a BUG SYMPTOM, not a termination condition.

**Proposed termination logic:**
1. At each step, attempt to find a valid rotation.
2. If no rotation is available AND the rhombus is in the nest-B region (geometrically), it has terminated successfully.
3. If no rotation is available AND the rhombus is NOT in nest-B, the algorithm has a bug.

The current two-phase approach (enter nest B, then wait for stuck) is actually the right structure, but `_is_stuck` has a flaw: it only checks if the opposite spoke has label 10 (for 180-degree) or if the 120-degree spokes exist. It doesn't account for the case where the spoke exists but leads to an impossible configuration.

## 3. Interaction Between Subsequent Rhombi and Already-Settled Rhombi

### The Wake Structure (Section 4)

Purbhoo's Section 4 introduces the concept of the WAKE: the set of tiles displaced during a single rhombus's journey. The wake has upper and lower components, each one tile wide, forming a path from the initial position to the final position.

The Wake Crossing Lemma is the key result: if rhombus X travels before Y, and Y begins above X's midline, then Y's journey does not cross the midline. In particular, Y does not disturb X's lower wake, and X's upper wake is intact for Y.

### What This Means for Multi-Rhombus Migration

The processing order matters critically. For south or west migration, the paper says to process from "minimal" to "maximal" in the standard order (page 13). The current code uses `flock_standard_order` which sorts by (label ascending, j descending, r descending).

The Wake Crossing Lemma guarantees that if rhombi are processed in this order, later rhombi do NOT disturb the positions of earlier ones. Specifically:

- Rhombus i's final position is STABLE -- no later rhombus j (j > i in processing order) will displace it.
- Rhombus j may travel through tiles that were rearranged by rhombus i's wake, but it will not cross into rhombus i's settled position.

The collisions reported in the HANDOFF (rhombi 3-6 landing at the same position as earlier rhombi) therefore indicate that either:
1. The direction rule is wrong, causing rhombi to take incorrect paths, OR
2. The nest-B initial state is wrong, causing rhombi to settle at wrong positions, OR
3. Both.

### The Already-Settled Rhombus Problem

When rhombus 1 migrates and settles at position (1,4), it becomes part of the tiling. Its edges are now 10-labels at its settled position. When rhombus 2 migrates, it encounters a tiling that includes rhombus 1 at (1,4). This is CORRECT -- rhombus 1's presence modifies the tiling, and rhombus 2's path should naturally route around it (the Wake Crossing Lemma guarantees this).

If instead rhombi 3-6 collide with rhombus 1, the Wake Crossing Lemma is being violated. This means either the lemma's hypotheses are not satisfied (wrong processing order, or the wake is not intact), or the rotation direction is wrong (so the path is wrong and crosses where it shouldn't).

## 4. Domain Boundary Effects on the Rotation Algorithm

### The Problem at the Boundary

The rotation algorithm requires examining spokes at a vertex. At interior vertices, all 6 spokes are well-defined. But at boundary vertices (near the edge of the hexagonal domain), some spokes may lie OUTSIDE the domain. The code's `_is_stuck` function checks `spokes[cw_spoke] in hm.edges` -- if a spoke's edge is not in the edge dictionary, the spoke "doesn't exist."

At the domain boundary, three things can happen:
1. **Missing spokes toward the exterior**: These correctly indicate no valid move in that direction.
2. **Missing spokes due to incomplete nest filling**: If the nest-B region was deleted (as the reproduction script does), spokes that SHOULD exist are missing, causing false "stuck" signals.
3. **Boundary spokes of the hexagonal region**: The hexagon A'AB'BC'C has specific edge orientations at each boundary segment. These determine which rotations are available near the boundary.

### The _complete_hex_domain Function

The code has a function `_complete_hex_domain` (line 144) that fills missing spokes with label 1 at vertices carrying rhombus edges. This is meant to handle boundary cases. But if nest B has been deleted, this function cannot correctly reconstruct what should be there.

### My Recommendation

The domain should be the COMPLETE mosaic as constructed by `from_puzzle`. No edges should be deleted from nest B. The boundary of the hexagonal region naturally constrains which rotations are available -- spokes that point outside the hexagon simply don't exist in the edge dictionary, and the algorithm correctly treats them as unavailable.

The issue is NOT that boundary spokes are missing. The issue is that the ARTIFICIAL deletion of nest-B edges creates an incorrect tiling that violates the assumptions of Purbhoo's algorithm.

## 5. Synthesis: The Root Cause

The collision bug has TWO interacting causes:

**Cause 1 (Primary): Wrong direction rule.** The 120-degree rotation direction uses a position-based heuristic (`direction_score`) instead of Purbhoo's orientation-based rule ("in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical"). This is the direction rule that Agent 3 should be analyzing.

**Cause 2 (Secondary): Wrong nest-B initialization.** The reproduction script deletes nest-B edges before migration, creating an incomplete tiling. Even if the direction rule were correct, migrating through an incomplete tiling would produce wrong results because:
- Rotation operations at vertices near the deleted region encounter missing edges
- The `_is_stuck` check triggers prematurely or at wrong positions
- The wake structure is corrupted because tiles that should be displaced don't exist

**The fix should be:**
1. Do NOT delete nest-B edges. Start with the complete mosaic from `from_puzzle`.
2. Fix the direction rule (Agent 3's domain).
3. Trust that the Wake Crossing Lemma guarantees non-collision when the direction rule and initial state are both correct.
4. The termination condition (`_is_stuck` while `_is_in_nest_b`) is structurally correct but should only apply when the mosaic is complete.

## 6. Response to Agent 3's Findings: The Geometry Question and the Table

Agent 3 makes several observations that interact with my boundary analysis. Let me address them.

**Agent 3's key discovery (Section 7):** Every rhombus in the equilateral grid has edges of ALL three orientations (H, NE, NW). This means Purbhoo's rule "some edge ends up horizontal or vertical" cannot be about the grid-level orientation -- it must be about the MOSAIC-level geometry. This is a crucial finding that I agree with entirely.

**Connection to the boundary problem:** This geometric subtlety is SEPARATE from the boundary initialization issue I identified, but the two bugs interact destructively. Even if the direction rule were perfectly correct, operating on an incomplete tiling (with nest B deleted) would produce wrong results. Conversely, even with a complete tiling, the wrong direction rule produces collisions. The test case likely needs BOTH fixes simultaneously.

**Agent 3's recommendation to use the `_MIGRATION_SHIFT` table:** I agree this is the immediate first step. The table at least makes the algorithm deterministic and enables systematic testing. But I add a caveat: the table should be tested with a COMPLETE mosaic (no nest-B deletion), not with the HANDOFF's reproduction script which deletes nest-B edges. The reproduction script's results are unreliable because they compound two distinct bugs.

**Agent 3's 60% confidence on NW -> H (CCW):** From the boundary perspective, I note that the NW -> H transition produces a rhombus whose shared edge is horizontal. At the boundary of nest B, horizontal edges are the CHARACTERISTIC orientation of nest-B rhombi (the `_place_h_rhombus` function places H-oriented rhombi in nest B). This is a structural argument for NW -> H: as the rhombus approaches nest B, transitioning to H aligns it with the nest's native orientation, facilitating "packing." This raises my confidence to about 70%.

However, Agent 3's observation about mosaic-level vs grid-level geometry means this argument from the nest orientation could be misleading. The mosaic rhombus is a 30/150-degree shape, not a 60/120-degree shape, so "alignment with the nest" is more nuanced than matching grid orientations.

## 7. What Lemma 2.6 Does NOT Say

One might think Lemma 2.6 (alpha = emptyset => unique mosaic with gamma = complement of beta) means that nest B should start empty during migration. This is wrong. Lemma 2.6 describes what the mosaic looks like when NO rhombi are at nest A -- it's about the structure of the tiling in the ABSENCE of a flock, not about the initial condition for migration.

Before migration, the mosaic has all three nests full: (alpha, beta, gamma). After migration, the flock has moved from A to B, producing a new mosaic where nest B contains both the original beta-rhombi (rearranged) and the migrated alpha-rhombi (which form the new flock zeta'). The paper proves that zeta' is a valid flock and can be identified with an LR-tableau.

The key mathematical point: migration does not "empty" nest A and "fill" nest B. It rearranges the ENTIRE tiling through a sequence of local rotations. The result is a new tiling of the same hexagonal region, with the flock now located at a different nest.
