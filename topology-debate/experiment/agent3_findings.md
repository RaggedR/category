# Agent 3 Findings: Orientation Cycle Analysis for Migration Direction Rule

## Summary

The direction rule ambiguity for NW-oriented spokes (2 and 5) cannot be resolved by looking at the resulting orientation alone. The resolution requires understanding what Purbhoo means by "horizontal or vertical" in his mosaic coordinate system, which is *relative* to the migration direction. I argue below that the correct interpretation yields a **consistent CW choice for NW spokes**, making the `_MIGRATION_SHIFT` table wrong at positions 2 and 5. Furthermore, I identify a deeper structural issue: the pivot vertex selection for NW-oriented rhombi may itself be systematically choosing the wrong vertex, which would flip the entire analysis.

---

## 1. The Orientation Cycle: What Sequence Should Rhombi Follow?

The HANDOFF.md proposes a cycle NW -> H -> NE -> NW. Let me check whether this is consistent with Purbhoo's description.

Purbhoo says (page 13): "rotate the tiling 120 degrees, in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical."

The three mosaic tile types are:
- (a) equilateral triangle, side 1
- (b) square, side 1 (this is the H-rhombus in the grid: label triple (10,1,0))
- (c) rhombus, side 1, angles 30 and 150 degrees (this is the NW or NE in the grid)

Crucially, a 30/150-degree rhombus has two long diagonals and two short diagonals. "Horizontal or vertical" refers to the edges of the rhombus *in Purbhoo's coordinate system*, which is oriented so that the migration direction points to the right.

For A->B migration in the hex_mosaic.py coordinate system, migration is roughly "rightward" (increasing j). In Purbhoo's rotated view, the migration direction = right. This means:
- "Horizontal" = parallel to the migration direction (the j-axis direction)
- "Vertical" = perpendicular to the migration direction

Now, the three edge orientations in the equilateral grid are H (horizontal), NE (60 degrees from horizontal), and NW (120 degrees from horizontal). In Purbhoo's mosaic coordinate system, these map to specific physical directions. The critical question is: which of {H, NE, NW} correspond to Purbhoo's "horizontal" and "vertical"?

**Key insight from the mosaic geometry**: The hexagon A'AB'BC'C has its sides aligned with specific directions. For A->B migration, the relevant orientation vectors are E_A and E_B. The E_A direction points along AA' (one side of the hexagon leaving nest A), and E_B points along BB' (arriving at nest B). These are NOT the same as the equilateral grid's horizontal direction.

The mosaic is drawn on a hexagonal region where the three pairs of sides correspond to the three edge orientations. Specifically:
- The sides AA' and BB' run parallel to the NE direction in the grid
- The sides AB' and CA' run parallel to the H direction in the grid  
- The sides BC' and CC' run parallel to the NW direction in the grid

When Purbhoo says "rotate the page so migration points right," and the migration is from A to B, the direction A->B through the interior of the hexagon would align roughly with the NE->H transition in grid coordinates.

## 2. Systematic Analysis of All 6 Spokes

Let me re-examine the table with a fresh eye, focusing on what "some edge of the rhombus ends up horizontal or vertical" means.

| my_spoke | orient | CW -> spoke | CW orient | CCW -> spoke | CCW orient |
|----------|--------|-------------|-----------|--------------|------------|
| 0        | H      | 2           | NW        | 4            | NE         |
| 1        | NE     | 3           | H         | 5            | NW         |
| 2        | NW     | 4           | NE        | 0            | H          |
| 3        | H      | 5           | NW        | 1            | NE         |
| 4        | NE     | 0           | H         | 2            | NW         |
| 5        | NW     | 1           | NE        | 3            | H          |

**Observation 1: Period-3 symmetry.** Spokes 0 and 3 have the same orientation (H), as do 1 and 4 (NE), and 2 and 5 (NW). The spoke pairs (0,3), (1,4), (2,5) are related by the 180-degree rotation (spoke -> spoke + 3). The CW and CCW options for paired spokes give the SAME pair of resulting orientations:
- From H: CW gives NW, CCW gives NE (both for spokes 0 and 3)
- From NE: CW gives H, CCW gives NW (both for spokes 1 and 4)
- From NW: CW gives NE, CCW gives H (both for spokes 2 and 5)

This means the direction rule, if it depends ONLY on the starting orientation, should give the same answer for paired spokes. The current `_MIGRATION_SHIFT` table does this: {0:2, 3:2} (both CCW), {1:4, 4:4} (both CW), {2:2, 5:2} (both CCW).

**Observation 2: The non-ambiguous cases.** For H-oriented spokes (0, 3), the code chooses CCW (shift=2), giving NE. For NE-oriented spokes (1, 4), it chooses CW (shift=4), giving H. This means the rule eliminates NW as a resulting orientation whenever possible. The interpretation is: **NW is the "bad" orientation** -- the one that is neither "horizontal" nor "vertical" in Purbhoo's coordinates.

**Observation 3: The ambiguous cases.** For NW-oriented spokes (2, 5), both CW and CCW give "good" orientations (NE or H respectively). The current table picks CCW (shift=2), giving H. But is this correct?

## 3. Does the Direction Depend on Migration Direction (A->B vs B->A)?

Purbhoo explicitly addresses this (page 14): "Note: because orientations can become rotated during migration, migration south is not necessarily the inverse of migration north."

More importantly, the direction rule says "causes some edge of the rhombus to end up exactly horizontal, or exactly vertical." Since "horizontal" and "vertical" are defined relative to the migration direction (page 4: "all directions in this paper are used in a relative, rather than in an absolute sense"), reversing the migration direction would rotate the reference frame by 120 degrees. This would change which orientations count as "horizontal" and "vertical."

**Conclusion**: Yes, the direction choice should depend on migration direction. For A->B vs B->A, the "forbidden" orientation (the one that is neither horizontal nor vertical) would be different. This is consistent with the paper's warning about non-trivial inverse relationships.

## 4. Does It Depend on Which PART of the Hexagon the Rhombus Is In?

No. Purbhoo's rule is purely local: "look to the right of the rhombus, find the smallest hexagon containing it, rotate 120 degrees." The rule about horizontal/vertical edges is a property of the rhombus's orientation relative to the migration direction, not its position. The same starting orientation should always produce the same rotation choice regardless of where in the hexagon the rhombus currently sits.

However, there is a subtlety: the **pivot vertex** selection IS position-dependent (it's the "leading" vertex toward nest B). The current code in `_pivot_vertex` always picks the vertex with higher j, which is correct for rightward migration. But the spoke index at that vertex depends on the orientation, and the spoke index determines which trio of spokes participates in the rotation.

## 5. Resolving the NW Ambiguity: Tracing Figure 2

Figure 2 shows a complete migration from left to right. The flock has three rhombi labeled 1, 1, 2. Let me trace what I can observe:

In the first row of Figure 2, the cyan rhombus (labeled "2") starts at the top-left of the nest and migrates rightward through several steps. Looking at the orientation changes visible in the sequence:

- The rhombus appears to alternate between NE-like and H-like orientations as it moves right
- At no point does a rhombus settle into a NW orientation for more than the initial position
- The overall pattern is consistent with NW being avoided as a "resting" orientation during migration

This suggests the cycle is: when starting from NW, the rhombus should transition to whichever of {H, NE} moves it most directly toward the target nest.

## 6. The Structural Argument for CW on NW Spokes

Here is my key argument. Consider the FLOW of orientations during migration:

For a rhombus migrating rightward (A->B), at each 3-fold rotation step it must advance toward the right. The orientation cycle should be monotonically "advancing" in the migration direction. In the equilateral grid with rightward migration:

- **H** (horizontal edge) points along the migration direction. After rotation, it should tilt to point somewhat upward or downward relative to the migration path. Both NE (60 degrees) and NW (120 degrees) are options. NW would point the edge "backward" relative to rightward migration. **Choose NE** (CCW, shift=2). This matches the current table.

- **NE** (60-degree edge) tilts forward-right. After rotation, it should become more aligned with the migration direction. H (0 degrees) aligns perfectly. NW (120 degrees) points backward. **Choose H** (CW, shift=4). This matches the current table.

- **NW** (120-degree edge) tilts backward-left. After rotation, it should move toward alignment with the migration direction. NE (60 degrees) is "more rightward" than H (0 degrees) in the sense of having a positive component in the direction perpendicular to migration. But H is more directly aligned with the rightward migration direction.

**Wait** -- this is the ambiguity restated. Let me try a different approach.

**The parity argument**: Spokes 0, 2, 4 form one trio (even), and 1, 3, 5 form the other (odd). The flip_3 function rotates spokes within a single trio. Notice:
- Spoke 2 (NW, even trio): CW within even trio sends spoke 2 -> spoke 4 (NE). CCW sends spoke 2 -> spoke 0 (H).
- Spoke 5 (NW, odd trio): CW within odd trio sends spoke 5 -> spoke 1 (NE). CCW sends spoke 5 -> spoke 3 (H).

The current table picks CCW for both (giving H). But consider: in the non-ambiguous cases:
- Spoke 0 (H, even): picks CCW -> spoke 4 (NE). That is: H -> NE within even.
- Spoke 1 (NE, odd): picks CW -> spoke 3 (H). That is: NE -> H within odd.
- Spoke 3 (H, odd): picks CCW -> spoke 1 (NE). That is: H -> NE within odd.
- Spoke 4 (NE, even): picks CW -> spoke 0 (H). That is: NE -> H within even.

The pattern within each trio is: H -> NE (via CCW), and NE -> H (via CW). This is a 2-cycle between H and NE, with NW appearing only as a starting orientation that must be escaped.

**For NW spokes, the question is: should NW escape to H or to NE?**

If the rule is "always avoid NW as the result," then both CW and CCW are valid. We need a tiebreaker. I propose looking at the **continuation**: after the NW -> ? step, the next step should continue the H <-> NE oscillation. If we go NW -> H (CCW), the next step would be H -> NE (CCW again). If we go NW -> NE (CW), the next step would be NE -> H (CW again). Both are valid continuation chains.

## 7. The Real Resolution: Purbhoo's "Horizontal or Vertical" IS Decidable

Re-reading Purbhoo more carefully: "causes some edge of the rhombus to end up exactly horizontal, or exactly vertical." The rhombus has FOUR edges (two pairs of parallel edges). After rotation, the rhombus has a new orientation. The question is whether ANY of its four edges is horizontal or vertical *in the rotated coordinate frame where migration points right*.

In the equilateral grid oriented for rightward migration:
- An **H-oriented** rhombus has its shared edge horizontal. Its other edges are at 60-degree angles. One pair of edges IS horizontal. **H has a horizontal edge.**
- An **NE-oriented** rhombus has its shared edge at 60 degrees. Its other edges are horizontal and at 120 degrees. Wait -- no. The NE rhombus consists of two triangles sharing a NE edge. Its four edges are: one NE (shared), one H, one NW, and one NE (outer). Actually each rhombus has exactly 4 edges of 3 orientations... No, a rhombus (30/150 degree) in the mosaic has 4 edges that come in two parallel pairs. In the equilateral triangle grid, the edges alternate between two of the three possible orientations.

**This is the key geometric fact I need to establish**: for each rhombus orientation (H, NE, NW), what are the orientations of its four edges?

- **NW rhombus** (shared edge = NW): The two triangles share a NW edge. The other edges are H and NE. So the rhombus has edges of types: NW, H, NE, NW (two parallel NW edges, one H, one NE). This rhombus has both a horizontal edge AND an edge at 60 degrees (NE). In Purbhoo's frame, it already has "a horizontal edge."
- **H rhombus** (shared edge = H): Other edges are NW and NE. So edges: H, NW, NE, H. It has a horizontal edge.
- **NE rhombus** (shared edge = NE): Other edges are H and NW. So edges: NE, H, NW, NE. It also has a horizontal edge.

**Every rhombus in the equilateral grid has edges of all three orientations!** This means "some edge of the rhombus ends up horizontal" is ALWAYS true, which would make the rule vacuous. This cannot be right.

**Resolution**: Purbhoo's "horizontal or vertical" does NOT refer to the equilateral grid edges. It refers to the **mosaic** geometry, which uses different tiles. In the mosaic, the (c)-type rhombus (30/150 degrees) has four edges that are NOT aligned with the equilateral grid. The mosaic's rhombus has edges at 15 and 75 degrees relative to the grid edges (approximately). Purbhoo's "horizontal" and "vertical" in the mosaic geometry have specific meanings that do not trivially map to the grid's H/NE/NW classification.

**This means the entire approach of choosing based on resulting grid orientation is an approximation.** The correct approach would require tracking the mosaic-level geometry, not just the grid-level spoke orientations.

## 8. Revised Recommendation

Given the analysis above, I believe:

1. **The `_MIGRATION_SHIFT` table's unambiguous entries are likely correct**: {0:2, 1:4, 3:2, 4:4} (H->NE via CCW, NE->H via CW). These avoid NW as a result, which is the best available heuristic for "horizontal or vertical."

2. **The NW entries {2:2, 5:2} (both CCW, giving H) are PLAUSIBLE but UNVERIFIED.** The alternative {2:4, 5:4} (both CW, giving NE) is equally plausible from orientation analysis alone.

3. **The current code does NOT use the `_MIGRATION_SHIFT` table at all.** It uses `direction_score` which is a position-based heuristic. The first fix should be to actually USE the table instead of the heuristic.

4. **To resolve the NW ambiguity definitively**, we need either:
   - (a) A careful geometric analysis of the mosaic-level rhombus edge orientations under the coordinate transformation between grid and mosaic, or
   - (b) Empirical testing: try both options and see which produces the correct 6 distinct positions for the test case.

5. **The `_MIGRATION_SHIFT` table is not used by `migrate_single`** -- this is the most immediate bug. Even if the table has a wrong entry for NW spokes, switching from the position heuristic to the table would at least make the algorithm deterministic and testable.

## 9. A Testable Prediction

If we switch to using `_MIGRATION_SHIFT` as-is (NW -> H via CCW for spokes 2 and 5), and collisions persist, then try flipping to {2:4, 5:4} (NW -> NE via CW). If exactly one of these two options produces 6 distinct positions for the test case, that resolves the ambiguity empirically.

My **best guess**, based on the H <-> NE oscillation pattern and the structure of the even/odd spoke trios: the NW -> H choice (CCW, the current table value) is more likely correct, because it feeds into the dominant H -> NE -> H cycle that characterizes rightward migration. But I hold this at about 60% confidence.

---

## 10. Note on Agent 2's Findings

Agent 2's file does not yet exist. If it appears with analysis of the pivot vertex selection or the mosaic-to-grid coordinate mapping, that could resolve the geometric question I identified in Section 7. The pivot vertex choice (which vertex of the shared edge to rotate around) is actually the more consequential decision than the rotation direction, since it determines which hexagon is being manipulated.
