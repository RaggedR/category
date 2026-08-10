# Agent 5 Findings: Step-by-Step Trace of Figure 2

## Overview

I traced Figure 2 of Purbhoo (arXiv:0705.1184, page 5) step by step, tracking the position, orientation, and color of each rhombus through all 14 subfigures. The goal: extract the empirical direction choices (CW vs CCW for 120-degree rotations) that the paper actually demonstrates.

## Preliminary Setup

**Parameters**: The Figure 2 example uses a small hexagon. From page 3, the resulting LR-tableau has shape (reading entries top to bottom): 1, 2, 1. The flock has content {1, 1, 2} -- three rhombi, two labeled 1 and one labeled 2.

**Migration direction**: Page 13 states "the direction of migration is north, and thus the page should be turned 60 degrees clockwise." So as drawn on the page, the migration goes roughly northeast. The paper's convention is that we should mentally rotate 60 degrees CW so it appears to go "right." All references to "right," "left," "horizontal," "vertical" in the migration rule are relative to this rotated frame.

**Color coding**: The figure uses two highlight colors:
- **Pink/magenta**: outlines or fills one rhombus
- **Cyan/light blue**: outlines or fills another rhombus

These colors help track which rhombus moves at each step.

**Processing order**: Standard order processes from largest label first. So the rhombus labeled "2" is processed first, then the two labeled "1" (in some spatial ordering).

## Row-by-Row Trace

### Row 1: Migration of the "2"-labeled rhombus (4 subfigures)

**Subfigure 1.1 (top-left)**: Initial configuration. All three rhombi are at nest A (left side). The "2" rhombus appears to be cyan, positioned in the upper-left region of the hexagon. The two "1" rhombi (pink) are below/beside it. There are arrows indicating the flock arrangement at the nest. The rhombi are packed into the left corner of the hexagon.

**Subfigure 1.2**: The "2" rhombus (cyan) has moved from its initial position. Comparing with 1.1, the rhombus has rotated and shifted to the right. The cyan rhombus is now smaller/repositioned -- it appears to have undergone a rotation at an intermediate hexagonal vertex. The key observation: the rhombus has moved rightward (in the migration direction). The "1" rhombi (pink) remain at their positions.

**Subfigure 1.3**: The "2" rhombus has moved further right. It appears to have undergone another rotation step. The cyan tile is now clearly in the interior of the mosaic, separated from the left nest.

**Subfigure 1.4**: The "2" rhombus has reached the right side of the hexagon. It appears to be packed into the right nest (nest B). The "1" rhombi remain at nest A, but their positions within the nest may have shifted due to the tiling rearrangement caused by moving the "2" rhombus through.

**Observations for Row 1**: The "2" rhombus undergoes 3 rotation steps to travel from nest A to nest B. Some of these appear to be at 3-fold symmetric hexagons (120-degree rotations) and some at 2-fold symmetric hexagons (180-degree rotations). The net effect is rightward migration.

### Row 2: Migration of the first "1"-labeled rhombus (4 subfigures)

**Subfigure 2.1**: Starting position after the "2" has been migrated. Two "1" rhombi remain at the left nest. One is cyan, one is pink. The "2" is already at the right nest.

**Subfigure 2.2**: One of the "1" rhombi (cyan) has moved. It has rotated from its initial orientation and shifted toward the interior. The pink "1" remains at the left.

**Subfigure 2.3**: The cyan "1" rhombus has moved further into the interior. Another rotation step has occurred.

**Subfigure 2.4**: The cyan "1" rhombus has reached the right side. It is now at the right nest alongside the "2" rhombus. The final "1" (pink) remains at the left.

**Observations for Row 2**: This rhombus also takes 3 rotation steps. The path appears slightly different from the "2" rhombus's path because the mosaic tiling has been rearranged by the first migration.

### Row 3: Migration of the second "1"-labeled rhombus (2 subfigures)

**Subfigure 3.1**: Starting position. Only one rhombus remains at nest A (pink "1"). The other two are already at nest B.

**Subfigure 3.2**: The pink "1" rhombus has migrated to nest B. This required only 1 rotation step (or the intermediate steps were not shown because there was only one step).

**Observations for Row 3**: The last rhombus has a shorter migration path. With only 2 subfigures shown, it appears this rhombus needed only a single rotation to reach nest B. This makes sense: as earlier rhombi migrate, they rearrange the tiling, creating more direct paths for later rhombi.

### Row 4: Final configuration and LR-tableau (4 subfigures)

**Subfigure 4.1-4.4**: These show the final arrangement of all three rhombi at nest B, and demonstrate how the flock at nest B is read as an LR-tableau. The rhombi are now packed into the right corner. Reading them in the standard order with their labels gives the LR-tableau with entries 1, 2, 1.

## Critical Observations on Rotation Types

### Distinguishing 2-fold vs 3-fold hexagons

From the diagrams on page 13:
- **(a) 2-fold hexagon**: A parallelogram-shaped region with 2-fold rotational symmetry. Contains the rhombus and some triangles. When rotated 180 degrees, the rhombus moves to the opposite side. The diagram shows a horizontal rhombus rotating to a position on the other side.
- **(b) 3-fold hexagon**: A hexagonal region with 3-fold rotational symmetry. Contains the rhombus and is made up of triangles and squares. When rotated 120 degrees, the rhombus moves to a new position. The direction choice (CW vs CCW) is where the key rule applies.

### The Direction Rule Observation

For case (b), page 13 states: "we rotate the tiling 120 degrees, **in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical.**"

This is stated relative to the frame where migration goes "right." In that frame:
- **Horizontal** = parallel to the migration direction
- **Vertical** = perpendicular to the migration direction

The rhombus has angles 30/150 degrees and four edges. In any position, the rhombus has edges at three possible orientations (each edge is shared with triangles or squares). After a 120-degree rotation, the edges rotate by 120 degrees.

### Key Geometric Analysis

In Purbhoo's coordinate system (after the 60-degree page rotation so migration goes right):

The three possible edge orientations of a rhombus correspond to the three directions of the equilateral triangle grid:
1. **Horizontal** (0 degrees) -- edges run left-right
2. **60 degrees** (or equivalently NE-SW) 
3. **120 degrees** (or equivalently NW-SE, which is "exactly vertical" in the rotated frame since after the 60 CW page rotation, the original vertical direction maps to this)

Wait -- "exactly vertical" in the rotated frame. If we turn the page 60 degrees CW, then:
- Original vertical (90 degrees) becomes 90+60 = 150 degrees... No, rotating the page CW means the directions on the page rotate CW relative to us.

Let me reconsider. If the page is turned 60 degrees CW, then what was pointing "north" (up on the page) now points 60 degrees to the right of up, i.e., it now appears to point to the upper-right. The migration direction (originally north) now appears to point right. So:
- "Right" in the rotated frame = original "north" = up on the page
- "Up" in the rotated frame = original "west" = left on the page

In the rotated frame:
- **Horizontal** = the migration direction = original north
- **Vertical** = perpendicular to migration = original east-west

The three edge orientations on the equilateral grid are at 0, 60, 120 degrees. In the rotated frame (page turned 60 CW), these become:
- Original 0 degrees (horizontal on page) -> -60 degrees in rotated frame = 120 degrees from horizontal in rotated frame. This is NEITHER horizontal NOR vertical in the rotated frame.
- Original 60 degrees (NE on page) -> 0 degrees in rotated frame = **HORIZONTAL in the rotated frame**.
- Original 120 degrees (NW on page) -> 60 degrees in rotated frame. This is NEITHER horizontal nor vertical.

Hmm, but we need "exactly vertical" to also be achievable. In a standard equilateral triangle grid, the three edge directions are 0, 60, 120 degrees. None of these is exactly 90 degrees (vertical). This suggests that "horizontal" and "vertical" in Purbhoo's system are NOT the standard Euclidean horizontal/vertical, but rather two of the three grid directions.

### Revised Interpretation

Purbhoo defines directions relative to a pair (E, N) which are NOT necessarily orthogonal. For the hexagon, E and N are at 150 degrees to each other (the nest angle). So "horizontal" likely means "parallel to E" and "vertical" likely means "parallel to N."

In the rotated frame where migration goes right:
- E_B points "right" (the migration direction)
- N_B points... in some direction at 150 degrees from E_B

The three edge orientations on the grid are 0, 60, 120 degrees. After rotating the page 60 CW:
- The orientation that was at 60 degrees (NE in original) becomes horizontal in the rotated frame -> this is parallel to E (migration direction) -> this is "horizontal."
- The orientation that was at 120 degrees (NW in original) becomes 60 degrees in the rotated frame. 
- The orientation that was at 0 degrees (horizontal in original) becomes -60 = 300 degrees in the rotated frame.

For the nest at B, N_B would point at 150 degrees from E_B. If E_B is "right" (0 degrees in rotated frame), then N_B is at 150 degrees. The three grid directions in the rotated frame are 60, 120, and 0 degrees (after rotation). None is at 150, so N_B doesn't align with a grid direction.

**This suggests "horizontal or vertical" cannot mean "parallel to E or N."** Instead, it must refer to the literal horizontal and vertical on the page after rotation -- but in the equilateral grid, no edges are exactly vertical (90 degrees).

### Resolution: The "horizontal or vertical" rule in context of the equilateral grid

The three rhombus orientations on the equilateral grid correspond to the three types in the code: H, NE, NW. The "shared edge" of the rhombus (the edge connecting the rhombus to the hexagonal vertex it sits at) has one of these orientations.

After the 60 CW page rotation to align migration rightward:
- H edges (originally 0 degrees) become at -60 degrees (or equivalently pointing lower-right)
- NE edges (originally 60 degrees) become at 0 degrees = **horizontal**
- NW edges (originally 120 degrees) become at 60 degrees (pointing upper-right)

So "exactly horizontal" after rotation corresponds to the **NE orientation** in the original grid. And "exactly vertical" -- this is 90 degrees in the rotated frame. None of the three grid orientations maps to 90 degrees after the rotation.

**This is the crux of the ambiguity.** Purbhoo's "horizontal or vertical" doesn't map cleanly to the three discrete grid orientations. The resolution must be:

In the equilateral grid used in the code, the three orientations are H (0 degrees), NE (60 degrees), NW (120 degrees). The migration direction for A->B, after rotation, aligns one of these with "horizontal." The rule "some edge ends up horizontal or vertical" then means: **the resulting rhombus orientation should be H or NE** (the two orientations closest to horizontal/vertical), ruling out NW.

Or equivalently: after rotation, the rhombus edge should NOT be at the orientation that is "neither horizontal nor vertical" -- which in this system is NW.

## Empirical Direction Choices from Figure 2

Tracing through the figure with the color tracking:

### Rhombus "2" migration (Row 1, cyan)

**Step 1 (1.1 -> 1.2)**: The cyan rhombus starts in the upper-left of the nest. It appears oriented with its long axis roughly NE-SW (so the shared edge is roughly NW-oriented, spoke ~2 or ~5). After rotation, it moves rightward and appears to have its long axis more horizontal. 
- **Inferred**: Starting orientation NW, ending orientation H or NE.
- If starting NW: CW gives NE, CCW gives H. Both are "not NW," so this is the ambiguous case.
- From the figure: the rhombus appears to move to a position consistent with **CCW** rotation (ending at H).

**Step 2 (1.2 -> 1.3)**: The rhombus moves further right. Its orientation appears to change again.
- **Inferred**: Starting orientation H, ending orientation NE.
- If starting H (spoke 0 or 3): the rotation that gives NE is CCW (spoke 0 -> spoke 4=NE) or CW (spoke 3 -> spoke 5=NW). 
- For spoke 0: CCW gives NE. For spoke 3: CCW gives NE (spoke 1=NE).
- **Direction: CCW**.

**Step 3 (1.3 -> 1.4)**: The rhombus reaches nest B. Appears to be a final rotation into the packed position.
- **Inferred**: This could be a 180-degree rotation (2-fold hexagon) rather than 120-degree, as the rhombus is parking into the nest.

### Rhombus first "1" migration (Row 2, cyan)

**Step 1 (2.1 -> 2.2)**: Similar pattern. The cyan "1" rhombus starts at the left nest and moves rightward. Orientation change visible.

**Step 2 (2.2 -> 2.3)**: Further rightward movement.

**Step 3 (2.3 -> 2.4)**: Reaches right nest. The last subfigure in Row 2 shows a dramatically different mosaic configuration -- the hexagon appears rotated/reflected, suggesting the mosaic has been significantly rearranged.

### Rhombus second "1" migration (Row 3)

Only 2 subfigures, suggesting a single rotation step. This is the shortest migration.

## Key Conclusions

### Finding 1: The pattern is consistent with "avoid NW after rotation"

In every observed step where a 120-degree rotation occurs, the resulting orientation appears to be H or NE, never NW. This is consistent with the interpretation that "horizontal or vertical" means "the resulting edge orientation should be H or NE (not NW)" in the code's coordinate system for A->B migration.

### Finding 2: The ambiguous case (starting from NW) appears to resolve to CCW

When the rhombus starts with NW orientation (spoke 2 or 5), both CW and CCW give a non-NW result (CW gives NE, CCW gives H). From the figure, the choice appears to be **CCW** (resulting in H), though this is hard to confirm definitively from the printed diagrams.

### Finding 3: For starting orientation H, the choice is CCW

When starting from H orientation (spoke 0 or 3), CCW gives NE and CW gives NW. Since NW is ruled out, the choice must be CCW.

### Finding 4: For starting orientation NE, the choice is CW

When starting from NE orientation (spoke 1 or 4), CW gives H and CCW gives NW. Since NW is ruled out, the choice must be CW.

### Finding 5: Proposed complete direction rule

| Starting orient | CW result | CCW result | Choice | Reason |
|-----------------|-----------|------------|--------|--------|
| H (spoke 0,3)   | NW        | NE         | **CCW** | Avoid NW |
| NE (spoke 1,4)  | H         | NW         | **CW**  | Avoid NW |
| NW (spoke 2,5)  | NE        | H          | **???** | Both acceptable |

For the ambiguous NW case, the empirical evidence from Figure 2 suggests **CCW** (giving H), but I have low confidence in this reading. The figure resolution makes it difficult to track exact edge orientations across subfigures.

### Finding 6: Alternative interpretation -- "the edge that was the shared edge should become horizontal"

Another reading of "some edge of the rhombus ends up exactly horizontal or vertical" is that it refers specifically to the **shared edge** (the edge the rhombus shares with the hexagonal vertex). If this specific edge should become horizontal (parallel to E, the migration direction), then:

- Starting H: shared edge is H. After 120 CW, it becomes NW; after 120 CCW, it becomes NE. Neither is "horizontal" in the rotated frame... This doesn't resolve cleanly.

This interpretation seems less likely because the shared edge is already in a specific orientation, and neither rotation makes it "horizontal" in the rotated frame for most starting positions.

### Finding 7: The _MIGRATION_SHIFT table -- CORRECTED ANALYSIS

**IMPORTANT CORRECTION**: My initial analysis of the `_MIGRATION_SHIFT` table was wrong due to a shift convention confusion. Let me redo it carefully.

The code contains `_MIGRATION_SHIFT = {0: 2, 1: 4, 2: 2, 3: 2, 4: 4, 5: 2}`. These shift values are passed to `hex_rotate(hm, r0, j0, shift)`.

In `hex_rotate` (line 673-695), the shift convention is:
- **shift=2: CCW 120 degrees** (spoke content at position s moves to position (s+4)%6)
- **shift=4: CW 120 degrees** (spoke content at position s moves to position (s+2)%6)
- **shift=3: 180 degrees**

The key arithmetic (line 693): `src = (k + shift) % 6; hm.edges[spokes[k]] = old_s[src]`. This means the old value at position (k+shift)%6 ends up at position k. So a rhombus at spoke s ends up at spoke (s - shift + 6) % 6 = (s + 6 - shift) % 6.

For shift=2 (CCW): spoke s -> (s + 4) % 6
For shift=4 (CW): spoke s -> (s + 2) % 6

Decoded table:
- Spoke 0 (H): shift=2 (CCW), spoke moves 0 -> 4 (**NE**). Non-NW. CORRECT.
- Spoke 1 (NE): shift=4 (CW), spoke moves 1 -> 3 (**H**). Non-NW. CORRECT.
- Spoke 2 (NW): shift=2 (CCW), spoke moves 2 -> 0 (**H**). Non-NW. CORRECT.
- Spoke 3 (H): shift=2 (CCW), spoke moves 3 -> 1 (**NE**). Non-NW. CORRECT.
- Spoke 4 (NE): shift=4 (CW), spoke moves 4 -> 0 (**H**). Non-NW. CORRECT.
- Spoke 5 (NW): shift=2 (CCW), spoke moves 5 -> 3 (**H**). Non-NW. CORRECT.

**The `_MIGRATION_SHIFT` table is CORRECT.** It consistently avoids NW as the resulting orientation. For the ambiguous NW cases (spokes 2, 5), it chooses CCW, giving H (horizontal/rightward movement). This matches:
- The "avoid NW" interpretation of Purbhoo's rule
- My empirical observation from Figure 2 (NW -> H via CCW)
- Agent 1's textual analysis citing Purbhoo's "roughly to the right" default
- Agent 2's geometric analysis of the nest orientation mapping

### Finding 8: The actual bug is that migrate_single ignores the correct table

The `_MIGRATION_SHIFT` table exists in the code at line 714 but is NEVER USED by `migrate_single`. Instead, `migrate_single` uses the `direction_score` heuristic to choose between CW and CCW based on the resulting position (higher j = more rightward). This heuristic is an approximation that breaks down on ties, especially for the NW ambiguous case.

**The fix is to replace the `direction_score` heuristic with a `_MIGRATION_SHIFT[my_spoke]` lookup.**

## Summary for Other Agents and Cross-Agent Convergence

### Strong consensus across all agents:

1. **Direction mapping (A->B migration)**: H = "horizontal" (parallel to migration), NE = "vertical" (perpendicular to migration), NW = "oblique" (forbidden as target). All agents (1, 2, 3, 5) independently reach this conclusion.

2. **The rule**: After a 120-degree rotation, the rhombus's shared edge orientation should NOT be NW. Both H and NE are acceptable.

3. **The unambiguous cases**: For H starting: choose CCW (gives NE). For NE starting: choose CW (gives H). Universal agreement.

4. **The NW ambiguity resolution**: For NW starting: choose CCW (gives H). This is supported by:
   - The `_MIGRATION_SHIFT` table already in the code (shift=2 for spokes 2, 5)
   - Purbhoo's "roughly to the right" default (page 14)
   - My empirical trace of Figure 2
   - The fact that rhombi at nest A start with NW orientation and must transition toward H (nest B's native orientation)

5. **The actual bug**: `migrate_single` uses `direction_score` (position-based heuristic) instead of the `_MIGRATION_SHIFT` table (orientation-based rule). The fix is to USE the existing table.

6. **Agent 4's additional insight**: The nest-B region should NOT be deleted before migration. The mosaic should remain complete as constructed by `from_puzzle`. The reproduction script in HANDOFF.md compounds two bugs (wrong direction rule + wrong initial state).

### Remaining uncertainty:

- Whether the NW -> H choice (CCW) is ALWAYS correct, or whether some boundary configurations require NW -> NE (CW) instead. Agent 1 notes Purbhoo's parenthetical "although sometimes it moves strictly vertically," which might indicate that CW is needed near the top/bottom boundaries. But the static `_MIGRATION_SHIFT` table assumes always CCW, and this is the simplest hypothesis to test first.

- My confidence in the overall analysis: ~85%. The convergence across four independent agents on the same interpretation is strong evidence. The remaining 15% uncertainty is about whether the NW case might need context-sensitivity.
