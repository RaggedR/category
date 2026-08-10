# Figure 2 Trace: Empirical Direction Choices in Purbhoo's Migration

## Preamble

This document traces every visible rotation step in Figure 2 of Purbhoo (arXiv:0705.1184, page 5), which shows a flock with content {1, 1, 2} migrating from the left nest to the right nest through 14 subfigures arranged in four rows. The goal is to extract the empirical 120-degree direction choices (CW vs CCW) at each step, which will resolve the ambiguity described in CONTEXT.md for NW-oriented spokes (spokes 2 and 5).

## Coordinate Frame Warning

Purbhoo states on page 13: "the direction of migration is north, and thus the page should be turned 60 degrees clockwise" for Figure 2. This means:
- As drawn on the page, the migration direction is roughly upper-right (northeast on the page).
- Purbhoo's "right" (the migration direction) corresponds to north in his coordinate system.
- When he says "horizontal or vertical", these are relative to the rotated frame where migration points right.

So "horizontal" in the direction rule means edges parallel to the migration direction (northeast as drawn), and "vertical" means edges perpendicular to it (northwest-southeast as drawn). This is critical for interpreting what the rule does.

## The Flock

The flock has three rhombi labelled 1, 1, 2. In the starting configuration (subfigure 1), they are packed in the left nest. The colors are:
- **Cyan (turquoise)**: the "1"-labelled rhombi (two of them)
- **Pink/magenta**: the "2"-labelled rhombus

The standard order for north/east migration says to process the maximal rhombus first. Since the flock has content (1, 1, 2), the maximal rhombus in standard order is the "2" rhombus, followed by the two "1" rhombi.

However, examining Figure 2 carefully: the migration processes ALL three rhombi together as a flock, showing the step-by-step rotation operations. The rhombi are migrated one at a time in order, and each row of subfigures shows several steps of migration.

## Tile Geometry Reminder

The mosaic tiles are:
- (a) equilateral triangle, side 1 -- has 3-fold symmetry about its center
- (b) square, side 1 -- has 2-fold rotational symmetry about each edge-center with a rhombus
- (c) rhombus, side 1, angles 30 and 150 degrees -- this is the tile that migrates

At each step, the algorithm looks to the right of the rhombus, finds the smallest hexagon containing it, and:
- If the hexagon has 2-fold symmetry (case (a) on page 13): rotate 180 degrees
- If 3-fold symmetry (case (b) on page 13): rotate 120 degrees, choosing the direction that makes some edge of the rhombus end up "exactly horizontal or exactly vertical"

## Row-by-Row Trace

### Row 1: Subfigures 1-4 (First rhombus migration)

**Subfigure 1 (starting position):** The initial mosaic. Three rhombi are packed in the left nest. The cyan "1" rhombus is at the top-left of the nest, with a second cyan "1" below it, and the pink "2" is positioned higher and to the right within the nest region. The orientation arrows (E single-head, N double-head) indicate the left nest orientation.

**Subfigure 1 to 2:** The first rhombus to migrate (processing in standard order -- the "2"-labelled pink rhombus, which is maximal) encounters a hexagon to its right. The pink rhombus visibly shifts rightward. Examining the orientation change: the rhombus rotates and moves toward the migration direction. The rhombus edges change orientation, consistent with a **120-degree rotation**. The pink rhombus moves from one orientation to another, with its long axis rotating. Comparing before and after: the rhombus appears to undergo a rotation that brings an edge toward the horizontal-in-migration-frame. This is a **3-fold case, and the direction chosen makes an edge horizontal (parallel to migration)**.

**Subfigure 2 to 3:** The pink "2" rhombus continues moving right. Another rotation step. The rhombus shifts further. Again a visible orientation change in the pink rhombus. The surrounding tiles rearrange. This step also appears to be a **120-degree rotation** (3-fold hexagon).

**Subfigure 3 to 4:** The pink "2" rhombus reaches the right side of the mosaic, near the right nest. One more rotation step. The pink rhombus has now moved significantly to the right. In subfigure 4, the pink rhombus appears to be at or near the right nest boundary. The two cyan "1" rhombi remain on the left side, largely unmoved.

### Row 2: Subfigures 5-8 (Second rhombus migration begins)

**Subfigure 5:** The pink "2" is now settled in the right nest region. One of the cyan "1" rhombi begins to migrate. The cyan rhombus that was next-to-maximal in standard order starts its journey rightward.

**Subfigure 5 to 6:** The active cyan "1" rhombus encounters a hexagon to its right and rotates. The other cyan "1" remains on the left. The pink "2" is stationary on the right. The rotation step shows the cyan rhombus changing orientation -- this appears to be a **3-fold (120-degree) rotation**. The cyan rhombus edges shift to become more aligned with the migration direction.

**Subfigure 6 to 7:** Another rotation step for the cyan rhombus. It continues moving rightward through the interior. The orientation changes again. Another **120-degree rotation step**.

**Subfigure 7 to 8:** The cyan rhombus reaches the right side. In subfigure 8, it appears to be near or at the right nest, positioned adjacent to the pink "2". The cyan rhombus has visibly changed orientation from its starting configuration. Notably, the remaining cyan "1" rhombus is still on the left.

### Row 3: Subfigures 9-10 (Third rhombus migration)

**Subfigure 9:** The last cyan "1" rhombus begins migrating. Both the pink "2" and the first cyan "1" are now settled on the right side. The remaining cyan "1" starts from the left nest.

**Subfigure 9 to 10:** The last cyan rhombus rotates and moves rightward. This appears to be a **120-degree step**, with the rhombus edges reorienting. The rhombus moves into the interior of the mosaic.

### Row 4: Subfigures 11-14 (Final steps and arrival)

**Subfigures 11-14:** The final row shows the last rhombus completing its migration. Through subfigures 11, 12, 13, 14, the last cyan "1" rhombus progresses rightward through several more rotation steps until all three rhombi are packed in the right nest. In the final subfigure (14), all three rhombi are on the right side, forming a flock in the right nest with the shape of the LR-tableau given on page 3:

```
| 1 |
| 2 |
| 1 |
```

## Key Empirical Observations on Direction Choices

### Observation 1: Predominance of 3-fold rotations

Nearly every visible step in Figure 2 involves a 3-fold (120-degree) rotation rather than a 2-fold (180-degree) rotation. The 2-fold case (where the hexagon has only 2-fold symmetry) appears to occur rarely or not at all in this example. This makes sense geometrically: the 2-fold case arises when the rhombus sits inside a hexagon made of two triangles and a square in a specific 2-fold symmetric configuration, which is less common in a typical mosaic interior.

### Observation 2: The direction rule -- edges becoming horizontal

At each 120-degree step, the rhombus orientation changes. The rule says to choose "whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical." Tracing the orientations visible in Figure 2:

- The rhombi in the LEFT nest have edges aligned with the left nest orientation (E_A, N_A).
- After migration, the rhombi in the RIGHT nest have edges aligned with (E_B, N_B).
- During transit, the orientation cycles through intermediate values.

The critical insight: in Purbhoo's rotated frame (where migration points right), "horizontal" means parallel to the migration direction and "vertical" means perpendicular. A 30/150-degree rhombus has three distinct edge orientations (each pair of parallel edges). After a 120-degree rotation, one pair of edges rotates to a new orientation. The rule selects the rotation that brings SOME edge pair to horizontal or vertical alignment in the migration frame.

### Observation 3: Consistent rightward drift

Every rotation step moves the rhombus to the right (in the migration direction). This is guaranteed by the rule, as stated on page 13-14: "the rhombus moves roughly to the right (although sometimes it moves strictly vertically)." The direction choice ensures net rightward progress.

### Observation 4: Orientation cycle pattern

Tracking the pink "2" rhombus across Row 1 (subfigures 1-4), its orientation appears to cycle through three states. In the equilateral-triangle grid used by hex_mosaic.py, these correspond to the three edge orientations H, NE, NW. The observed cycle, consistent with rightward migration, appears to follow:

**NW -> H -> NE -> NW -> ...**

This is exactly the cycle hypothesized in HANDOFF.md. Each 120-degree rotation advances the orientation one step in this cycle. For the hex_mosaic.py spoke system:
- spoke 2 (NW): should go to spoke 0 (H), which is CCW (shift -2)
- spoke 0 (H): should go to spoke 4 (NE), which is also CCW (shift -2, i.e., 0-2 = -2 mod 6 = 4)
- spoke 4 (NE): should go to spoke 2 (NW), which is also CCW (shift -2, i.e., 4-2 = 2)

Wait -- this gives a uniform CCW rule for spokes 0, 2, 4. But from CONTEXT.md:
- spoke 0 (H): CCW gives spoke 4 (NE) -- MATCHES the cycle NW->H->NE
- spoke 2 (NW): CCW gives spoke 0 (H) -- MATCHES the cycle NW->H->NE  
- spoke 4 (NE): CCW gives spoke 2 (NW) -- MATCHES the cycle ...->NE->NW->H->...

But CONTEXT.md also lists spokes 1, 3, 5 (the odd spokes). These have:
- spoke 1 (NE): CW gives spoke 3 (H), CCW gives spoke 5 (NW)
- spoke 3 (H): CW gives spoke 5 (NW), CCW gives spoke 1 (NE)
- spoke 5 (NW): CW gives spoke 1 (NE), CCW gives spoke 3 (H)

For the same NW->H->NE cycle on odd spokes:
- spoke 5 (NW): go to H = spoke 3, which is CW (shift +2)
- spoke 3 (H): go to NE = spoke 1, which is CCW (shift -2)
- spoke 1 (NE): go to NW = spoke 5, which is CW (shift +2)

This is NOT uniform! The direction alternates between CW and CCW depending on even/odd spoke parity.

### Observation 5: Resolution of the NW spoke ambiguity

The CONTEXT.md identified spokes 2 and 5 (both NW orientation) as ambiguous. From the figure trace and the cycle analysis:
- **Spoke 2 (NW, even):** CCW rotation -> spoke 0 (H). This is the correct choice for rightward migration.
- **Spoke 5 (NW, odd):** CW rotation -> spoke 1 (NE). This is the correct choice for rightward migration.

Both choices advance the orientation in the NW->H->NE direction, but they require DIFFERENT rotation directions because of the even/odd spoke indexing. This resolves the ambiguity: the direction rule is NOT simply "always CW" or "always CCW" -- it depends on the spoke parity.

## Summary Table: Direction Rule for A->B Migration

| spoke | orient | correct rotation | new spoke | new orient | rule |
|-------|--------|-----------------|-----------|------------|------|
| 0     | H      | CCW (-2)        | 4         | NE         | advance cycle |
| 1     | NE     | CW (+2)         | 3         | H          | advance cycle |
| 2     | NW     | CCW (-2)        | 0         | H          | advance cycle |
| 3     | H      | CCW (-2)        | 1         | NE         | advance cycle |
| 4     | NE     | CW (+2)         | 0         | H          | advance cycle |
| 5     | NW     | CW (+2)         | 1         | NE         | advance cycle |

Wait -- spokes 3 and 4 need re-examination. Let me reconsider.

The cycle NW->H->NE means: from NW, we want H; from H, we want NE; from NE, we want NW (to restart the cycle and continue advancing rightward).

- spoke 0 (H) -> want NE: CCW gives spoke 4 (NE). Correct: CCW.
- spoke 1 (NE) -> want NW: CW gives spoke 3 (H)... no, want NW. CCW gives spoke 5 (NW). Correct: CCW.
- spoke 2 (NW) -> want H: CW gives spoke 4 (NE), CCW gives spoke 0 (H). Correct: CCW.
- spoke 3 (H) -> want NE: CW gives spoke 5 (NW), CCW gives spoke 1 (NE). Correct: CCW.
- spoke 4 (NE) -> want NW: CW gives spoke 0 (H), CCW gives spoke 2 (NW). Correct: CCW.
- spoke 5 (NW) -> want H: CW gives spoke 1 (NE), CCW gives spoke 3 (H). Correct: CCW.

**Revised conclusion: the direction is uniformly CCW (shift -2) for ALL spokes when migrating rightward (A->B direction).** My earlier analysis had the cycle direction wrong for some spokes. The cycle is NW->H->NE->NW and EVERY spoke achieves this by shifting -2 (CCW).

This would mean `_MIGRATION_SHIFT` should be uniformly {0: -2, 1: -2, 2: -2, 3: -2, 4: -2, 5: -2} for A->B migration, or equivalently, always use CCW 120-degree rotation.

## Remaining Uncertainty

I am approximately 70% confident in the uniform-CCW conclusion. The figure is detailed but small, and precisely tracking individual edge orientations through 14 subfigures at PDF resolution is difficult. Two sources of doubt remain:

1. **The 2-fold vs 3-fold distinction**: I may be missing some 180-degree steps in the figure that would change the analysis.
2. **The frame rotation**: Purbhoo's warning about relative directions means "horizontal" and "vertical" depend on which way migration points. If I have the frame rotation wrong by 60 degrees, the preferred orientations would shift.

The strongest test would be to implement the uniform-CCW rule and check whether it resolves the collision bug for the test case in HANDOFF.md.

## Note on Vanilla Center

The vanilla_findings.md file does not yet exist. Once the center agent has posted their synthesis, the uniform-CCW hypothesis should be cross-checked against any algebraic analysis they provide. In particular, the question of whether the direction rule depends on the migration direction (north vs east vs south vs west) or whether it is always CCW in some canonical sense would benefit from the vanilla center's perspective.
