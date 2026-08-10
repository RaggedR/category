# Vanilla Center: Final Synthesis of the Migration Direction Rule Bug

*Incorporating findings from all five specialists: Paper Exegete, Coordinate Translator, Orientation Analyst, Figure 2 Tracer, and Boundary/Termination Analyst.*

---

## Executive Summary

The migration direction rule (CW vs CCW at 120-degree rotations) is **almost certainly NOT the bug**. The current code's direction logic is correct for all reachable spoke indices. The collision bug resides in the **domain setup and termination logic** -- specifically, how nest B is prepared before migration and how the rhombus knows when to stop. The direction rule analysis produced a genuine mathematical insight (the NW-spoke ambiguity and its resolution via the two-orientation constraint), but this insight points away from the direction rule as the bug's cause.

---

## Part I: The Direction Rule Is Correct

### 1.1 Textual Analysis (All Specialists Agree)

Purbhoo's direction rule for 3-fold symmetric hexagons (page 13):

> "Rotate 120 degrees in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical."

All five specialists agree on the mapping for A-to-B migration in the hex_mosaic.py coordinate system:
- **H** = "horizontal" (parallel to migration direction, rightward) -- ALLOWED
- **NE** = "vertical" (perpendicular to migration) -- ALLOWED
- **NW** = oblique (neither horizontal nor vertical) -- FORBIDDEN

The Coordinate Translator provides the structural explanation: each nest uses one orientation (A=NW, B=H, C=NE), and the forbidden orientation NW is the source-nest orientation. The Exegete strengthens this via Proposition 3.1: during migration, rhombi only ever occupy two of three orientations, ruling out one nest.

### 1.2 The NW Spoke Ambiguity and Its Resolution

The central mathematical finding of this analysis is that Purbhoo's stated rule appears to be ambiguous for NW-oriented spokes (2 and 5): both CW and CCW rotations produce allowed orientations (NE and H respectively). This was identified independently by all specialists who analyzed the spoke table.

**The Orientation Analyst's resolution is definitive**: the `_pivot_vertex` function only returns spokes 2, 3, and 4. Spokes 0, 1, and 5 are never the starting spoke at any migration step. The ambiguity at spoke 5 is therefore moot -- it never occurs. For spoke 2 (the only reachable NW spoke), the code chooses CCW (giving H), which is the geometrically correct choice: it transitions the rhombus from source-nest orientation (NW=A) toward target-nest orientation (H=B).

The resulting direction rule for A-to-B migration:

| Reachable spoke | Orient | Direction | Result orient |
|---|---|---|---|
| 2 | NW | CCW | H |
| 3 | H | CCW | NE |
| 4 | NE | CW | H |

The code's logic `direction = 1 if my_spoke in (1, 4) else -1` is correct. Including spoke 1 is harmless (unreachable). The `_MIGRATION_SHIFT` table agrees for all relevant entries.

### 1.3 The Tracer's Error

The Figure 2 Tracer initially found results consistent with the above, but then revised to claim uniform CCW rotation for all spokes. This is incorrect: for spoke 1 (NE), CCW gives spoke 5 (NW), which is the FORBIDDEN orientation. The tracer's revised analysis assumes a continuous NW->H->NE->NW cycle, but Purbhoo's rule prevents the NE->NW transition. The cycle is not NW->H->NE->NW; rather, after transitioning out of NW (the initial nest-A orientation), the rhombus oscillates between H and NE, never returning to NW. The tracer acknowledges 70% confidence and the difficulty of reading Figure 2 at PDF resolution.

### 1.4 Why Purbhoo's Claim of Unambiguity Is Consistent

Purbhoo says "this restriction rules out one of the two possible 120 degree rotations." The Exegete's Proposition 3.1 analysis explains why: during a correctly proceeding migration, the rhombus is never in the NW orientation at a 3-fold hexagon (except possibly at the very first step when leaving the source nest). Therefore, in practice, the rule always disambiguates: one rotation produces an allowed orientation, the other produces NW. The NW-spoke case, where both rotations produce allowed orientations, simply does not arise during normal operation because it would require the rhombus to be in the forbidden orientation already. If a rhombus IS at NW after the first step, something has already gone wrong.

---

## Part II: The Actual Bug -- Domain Setup and Termination

### 2.1 The Boundary Analyst's Diagnosis

The Boundary/Termination Analyst provides the most operationally significant findings. The collision bug is caused by incorrect preparation of the nest-B domain:

**Failed approach 1 (from HANDOFF): Leave nest-B edges undefined.** The `ensure_hex_cell` function dynamically fills missing edges with label 1, causing the domain to grow incrementally as the rhombus approaches the boundary. This creates an unstable boundary that keeps retreating, allowing the rhombus to overshoot indefinitely.

**Failed approach 2 (from HANDOFF): Pre-place nest-B rhombi.** If `right_nest_positions` places 10-labelled edges in nest B, migrating rhombi encounter these as 2-fold symmetric hexagons and trigger 180-degree swaps. This cascades through the pre-placed rhombi, displacing them from their intended positions.

**Correct approach: Pre-fill nest B with finite all-1 domain.** Before migration, nest B should contain only 1-labelled edges (empty triangles), bounded by a correctly computed extent. The extent is determined by the complement of mu: beta = complement(mu, n, m) defines how many rhombus positions are available at each level.

### 2.2 Termination

The current termination logic uses `_is_in_nest_b` AND `_is_stuck`. The Boundary Analyst identifies the critical issue: if the domain extends too far (filled with 1s indefinitely), `_is_stuck` never returns True because there are always more all-1 cells to rotate through. If the domain is correctly bounded, `_is_stuck` works because spokes at the boundary are outside the domain (return None), blocking further rotation.

**Refined termination**: A rhombus has arrived when it is in nest B AND the pivot vertex sits at the domain boundary (at least one required spoke is absent from `hm.edges`). The Boundary Analyst also suggests a purely geometric criterion: the pivot vertex (cr, cj) satisfies cj >= N - cr + width_at_level(cr), where width is derived from the beta partition shape.

### 2.3 Dynamic Domain Extension Is the Root Cause

The `ensure_hex_cell` call at line 962 of `migrate_single` is identified by the Boundary Analyst as the most dangerous operation. It fills missing edges with label 1 on every migration step, effectively extending the domain in front of the rhombus. This prevents termination in nest B. The fix: either remove this call for boundary vertices or pre-compute the complete domain and guard against extension beyond it.

### 2.4 Processing Order

The Boundary Analyst raises the processing-order question. Purbhoo specifies: for south/west migration, process from minimal to maximal; for north/east migration, from maximal to minimal. The code uses `flock_standard_order`. For A-to-B migration in the code's coordinate system (leftward to rightward = roughly east), the maximal rhombus should be processed first. If `flock_standard_order` gives minimal-first, the order is wrong and could cause collisions.

---

## Part III: Additional Issues

### 3.1 The `rotate_120` Function and Partial Domains

The `rotate_120` function (line 649) only permutes spokes that are in the domain, silently skipping out-of-domain spokes. If 4 of 6 spokes are in the domain, the rotation permutes those 4 and leaves the other 2 unchanged. This produces an invalid local state -- the three edge labels around each face must satisfy the DELTA/NABLA pattern, and a partial rotation breaks this. This could corrupt the mosaic state near boundaries, causing downstream collisions.

### 3.2 The `hex_rotate` Function as Alternative

There is a separate function `hex_rotate` (line 673) that permutes all 12 edges (6 spokes + 6 outers) and treats out-of-domain edges as label 1. This is more robust than `rotate_120` for boundary cases because it does not skip any positions. However, `migrate_single` uses `rotate_120`, not `hex_rotate`. Switching to `hex_rotate` (or a variant that also handles outers) might fix boundary rotation issues.

### 3.3 Pivot Vertex for NE Orientation

The `_pivot_vertex` function returns spoke 4 for NE-oriented rhombi, at center (r+1, j+1). But the NE edge connects (r, j+1) and (r+1, j+1). The "rightward leading" vertex should be the one with higher j, which is (r+1, j+1) -- but both vertices have the same j-value (j+1)! The tiebreak uses higher r, giving (r+1, j+1). This seems correct for rightward migration (the upper vertex is further along the migration path because rows narrow as r increases), but it should be verified against Figure 2.

---

## Part IV: The Correct Fix (Ranked)

1. **Pre-compute the nest-B domain** with correct extent based on complement(mu). Fill all edges within this domain with label 1. Do not place any 10-labelled edges in nest B. (Boundary Analyst, Priority: CRITICAL)

2. **Remove or guard `ensure_hex_cell`** at line 962. It must not extend the domain beyond the pre-computed boundary. Replace with a domain membership check. (Boundary Analyst, Priority: CRITICAL)

3. **Verify processing order** in `flock_standard_order`. For A-to-B migration (east-like), process maximal rhombus first. (Boundary Analyst, Priority: HIGH)

4. **Add assertions**: After each 120-degree rotation, assert the resulting orientation is H or NE, not NW. Assert that `_pivot_vertex` returns spoke 2, 3, or 4. (Orientation Analyst, Priority: MEDIUM -- for debugging)

5. **Keep the current direction rule**. The logic `direction = 1 if my_spoke in (1, 4) else -1` is correct. (Orientation Analyst, Coordinate Translator, Exegete, Priority: DO NOT CHANGE)

6. **Consider switching from `rotate_120` to `hex_rotate`** for boundary robustness, or fixing `rotate_120` to handle partial domains correctly. (Priority: MEDIUM)

---

## Part V: Confidence Assessment

- The direction rule is correct for all reachable spokes: **95% confidence**. Four of five specialists agree (the fifth, the Tracer, self-reported 70% confidence in a contradictory conclusion that contains a demonstrable logical error).

- The bug is in domain setup / termination, not in the direction rule: **85% confidence**. The Boundary Analyst's analysis is structurally sound and explains both the "overshooting" and "cascade" failure modes described in the HANDOFF. However, I cannot be certain without testing.

- The specific fix (pre-computed finite all-1 nest-B domain + guarded ensure_hex_cell): **75% confidence**. This is the most natural approach, but the exact extent computation needs verification, and there may be subtleties in how the boundary interacts with the rotation mechanics.

---

*This synthesis resolves the star topology experiment by unifying all five specialist perspectives. The key mathematical insight -- the two-orientation constraint from Proposition 3.1 resolving the NW spoke ambiguity -- is a genuine contribution to understanding Purbhoo's algorithm. The practical conclusion is that the direction rule code is correct and the bug lies in domain preparation.*
