# Agent 6 (Synthesis): Unified Fix for the Purbhoo Migration Direction Rule

## Status of Inputs

At the time of writing, no findings from Agents 1-5 have appeared. This synthesis is therefore derived entirely from primary sources: Purbhoo's paper (arXiv:0705.1184, Sections 2-5), the existing code in `hex_mosaic.py`, the CONTEXT.md analysis, the HANDOFF.md, and the FINDINGS.md hypotheses S1 and S5. If agent findings appear later, a revised synthesis should incorporate any conflicting evidence.

---

## A. The Direction Rule Fix

### The Core Problem

Purbhoo's direction rule (page 13) states: when a hexagon has 3-fold rotational symmetry (case (b)), "rotate the tiling 120 degrees, *in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical.*"

The current code ignores this geometric rule entirely and instead uses a position-based heuristic (`direction_score`), which picks whichever rotation moves the rhombus to a higher `j` coordinate (more rightward), breaking ties on `r`. This heuristic sometimes gives the wrong answer, particularly for NW-oriented spokes.

### What "horizontal" and "vertical" mean

Purbhoo warns explicitly (page 4): "all directions in this paper are used in a relative, rather than an absolute sense." And on page 13: "we will assume that our diagrams are turned so that the direction of migration is pointing directly to our right."

So "horizontal" and "vertical" are defined relative to the migration direction. For A-to-B migration (which the code treats as rightward/eastward in the grid), "horizontal" means parallel to the migration direction (east-west) and "vertical" means perpendicular to it (north-south).

In the equilateral triangle grid used by `hex_mosaic.py`, there are three edge orientations:
- **H** (horizontal): edges like `((r,j), (r,j+1))` -- these run east-west
- **NE** (northeast): edges like `((r,j), (r+1,j))` -- these run at 60 degrees
- **NW** (northwest): edges like `((r,j), (r+1,j-1))` -- these run at 120 degrees

Now comes the critical geometric question. In Purbhoo's mosaic coordinate system (not the equilateral grid), the hexagon A'AB'BC'C has a specific geometry. The mosaic tiles include 30/150-degree rhombi. When the mosaic is "straightened" into the equilateral triangle grid, the three edge orientations H, NE, NW of the grid correspond to three particular orientations in the mosaic.

The comment at line 710-713 of the code states the intended interpretation:
```
# "Rotate in whichever way causes an edge to end up horizontal or vertical."
# In mosaic geometry: H = horizontal, NE = vertical, NW = neither.
# After rotation, shared edge must be H or NE, not NW.
```

This is the key interpretive claim: that "horizontal" in Purbhoo maps to H-edges in the grid, "vertical" maps to NE-edges, and NW-edges are the "neither" case that should be avoided. This interpretation is geometrically sound: when we orient the page so migration points right, H-edges are literally horizontal, NE-edges are the closest to vertical (at 60 degrees from horizontal), and NW-edges point backward relative to the migration direction.

### The Spoke-by-Spoke Direction Table

Using the analysis from CONTEXT.md (which is correct), for each starting spoke, we can determine which rotation yields an H or NE result (good) vs NW (bad):

| Starting Spoke | Orient | CW -> Spoke (Orient) | CCW -> Spoke (Orient) | Correct Choice | Reasoning |
|---|---|---|---|---|---|
| 0 | H | 2 (NW) | 4 (NE) | **CCW** | NE is good, NW is bad |
| 1 | NE | 3 (H) | 5 (NW) | **CW** | H is good, NW is bad |
| 2 | NW | 4 (NE) | 0 (H) | **AMBIGUOUS** | Both NE and H are good |
| 3 | H | 5 (NW) | 1 (NE) | **CCW** | NE is good, NW is bad |
| 4 | NE | 0 (H) | 2 (NW) | **CW** | H is good, NW is bad |
| 5 | NW | 1 (NE) | 3 (H) | **AMBIGUOUS** | Both NE and H are good |

The unambiguous cases (spokes 0, 1, 3, 4) are clear: always pick the rotation that avoids NW. The ambiguous cases (spokes 2 and 5, both starting from NW orientation) are where the bug lives.

### Resolving the NW-Spoke Ambiguity

Several lines of evidence converge on the resolution:

**Evidence 1: The `_MIGRATION_SHIFT` table.** The code contains a pre-computed table (line 714) that was intended to encode the correct direction:
```python
_MIGRATION_SHIFT = {0: 2, 1: 4, 2: 2, 3: 2, 4: 4, 5: 2}
```
A shift of 2 in `hex_rotate` means CCW 120 degrees. A shift of 4 means CW 120 degrees. Converting: spokes 0,2,3,5 use CCW; spokes 1,4 use CW. For the ambiguous cases: spoke 2 uses CCW (giving H), spoke 5 uses CCW (giving H). This table says: **when starting from NW, always rotate to H**.

**Evidence 2: Migration direction analysis.** For rightward migration, the rhombus should advance to the right at each step. An NW-oriented rhombus has its long axis pointing upper-left to lower-right. To advance it rightward:
- Rotating to H (horizontal) moves the long axis to point east-west -- directly aligned with the migration direction. This is the natural "next step" toward the target nest.
- Rotating to NE moves the long axis to point from lower-left to upper-right -- partially perpendicular to migration direction.

The H rotation is the more natural "advance" step. This supports spoke 2 -> CCW (giving H) and spoke 5 -> CCW (giving H).

**Evidence 3: Orientation cycling.** The HANDOFF.md proposed an orientation cycle NW -> H -> NE -> NW for rightward migration. Under this cycle:
- Spoke 2 (NW): next should be H -> choose CCW (giving H). Correct.
- Spoke 5 (NW): next should be H -> choose CCW (giving H). Correct.
- Spoke 0 (H): next should be NE -> choose CCW (giving NE). Correct.
- Spoke 1 (NE): next should be NW... but we said NW is bad?

This reveals a subtlety: the cycle NW -> H -> NE -> NW doesn't mean we always advance one step in the cycle. It means that when facing a 3-fold choice, the rule picks the non-NW option, and the cycle emerges as a consequence of the geometry. The "next" orientation after NE would naturally be NW (returning to start), but at that point the rhombus encounters a new hexagon and the process continues.

**Evidence 4: Proposition 3.1.** Purbhoo proves that all rhombi in the migrated flock end up in the same nest. The proof (page 14) notes: "the rhombus can only ever be in two of three orientations during the migration process, thereby ruling out one of the nests." For A-to-B migration (rightward), the two allowed orientations during transit should be H and NE (the "horizontal" and "vertical" ones). NW is the orientation the rhombus starts with in nest A but should not return to during transit -- it would correspond to the rhombus trying to go back toward nest A or toward nest C.

### The Definitive Direction Lookup

Combining all evidence, the correct direction for each spoke is:

| Spoke | Starting Orient | Correct Rotation | Resulting Orient | Shift Value |
|---|---|---|---|---|
| 0 | H | CCW (-1) | NE | 4 (as shift) or -2 (as step) |
| 1 | NE | CW (+1) | H | 2 (as shift) or +2 (as step) |
| 2 | NW | CCW (-1) | H | 4 (as shift) or -2 (as step) |
| 3 | H | CCW (-1) | NE | 4 (as shift) or -2 (as step) |
| 4 | NE | CW (+1) | H | 2 (as shift) or +2 (as step) |
| 5 | NW | CCW (-1) | H | 4 (as shift) or -2 (as step) |

Wait -- I need to re-examine the shift values. The `_MIGRATION_SHIFT` table says `{0: 2, 1: 4, 2: 2, 3: 2, 4: 4, 5: 2}`. But `hex_rotate` with shift=2 is CCW, and `rotate_120` with direction=+1 is CW (spoke k moves to spoke (k+2)%6). These are DIFFERENT conventions.

Let me reconcile. In `rotate_120`:
- `direction=+1` (CW): value at spoke k moves to spoke (k+2)%6. So if rhombus is at spoke s, after CW it's at spoke (s+2)%6.
- `direction=-1` (CCW): value at spoke k moves to spoke (k-2)%6. So rhombus at spoke s moves to (s-2)%6 = (s+4)%6.

In `hex_rotate`:
- `shift=2` (labeled CCW): spoke k gets value from spoke (k+2)%6. This means value from old spoke (k+2) goes to new spoke k. So rhombus at old spoke s ends up at new spoke (s-2)%6 = (s+4)%6. This is the SAME as `rotate_120(direction=-1)`.

So `_MIGRATION_SHIFT = {0: 2, 1: 4, ...}` with `hex_rotate`:
- shift=2 (CCW in hex_rotate) = direction=-1 in rotate_120. Rhombus at spoke s -> spoke (s+4)%6.
- shift=4 (CW in hex_rotate) = direction=+1 in rotate_120. Rhombus at spoke s -> spoke (s+2)%6.

Re-examining `_MIGRATION_SHIFT`:
- spoke 0, shift=2: new spoke = (0+4)%6 = 4. Orient of spoke 4 = NE. Good.
- spoke 1, shift=4: new spoke = (1+2)%6 = 3. Orient of spoke 3 = H. Good.
- spoke 2, shift=2: new spoke = (2+4)%6 = 0. Orient of spoke 0 = H. Good (resolves ambiguity toward H).
- spoke 3, shift=2: new spoke = (3+4)%6 = 1. Orient of spoke 1 = NE. Good.
- spoke 4, shift=4: new spoke = (4+2)%6 = 0. Orient of spoke 0 = H. Good.
- spoke 5, shift=2: new spoke = (5+4)%6 = 3. Orient of spoke 3 = H. Good (resolves ambiguity toward H).

**The `_MIGRATION_SHIFT` table is correct.** The problem is that `migrate_single` does not use it.

The equivalent lookup in terms of `rotate_120` direction:

```python
_DIRECTION_FOR_SPOKE = {0: -1, 1: +1, 2: -1, 3: -1, 4: +1, 5: -1}
```

Pattern: NE spokes (1, 4) use CW (+1). All others (H and NW: spokes 0, 2, 3, 5) use CCW (-1).

---

## B. Boundary and Termination Fix

### Nest B Preparation

Purbhoo's algorithm (Corollary 3.2, page 14) gives the flock on nest A the orientation (-E_A, -N_A) and migrates from A to B. The flock at nest B is **created by the migration itself** -- it does not pre-exist.

This means:
1. **Remove all pre-existing rhombi from nest B before migration.** The nest-B region should contain only triangular tiles (all edges = 1 in the triangular grid representation) before migration starts.
2. **Keep nest-A rhombi** (they are the source flock) and **keep nest-C rhombi** (they are part of the mosaic boundary data gamma).
3. The hexagonal domain must extend far enough to include the nest-B positions where rhombi will land. All edges in this region should be initialized to 1 (empty/triangular tiles).

In the code, this means: before calling `migrate_flock`, ensure the nest-B area has all edges defined as 1, with no 10-edges (rhombus edges). The HANDOFF.md already identifies this: "Don't fill nest B with pre-placed rhombi -- that causes 180 degree swap cascades."

### Termination Condition

The current `_is_stuck` function checks whether the pivot vertex has no 10 on the opposite spoke and both 120-degree spokes are undefined. This is too simplistic.

Purbhoo (page 14): "Repeat this process, until the rhombus is packed into a new nest." A rhombus is "packed" when it is in the nest-B region and is part of a valid transformed diagram. Concretely:

1. A rhombus has entered nest B when its position satisfies the nest-B geometric criterion (the existing `_is_in_nest_b` function).
2. A rhombus is "settled" (packed) when it can no longer advance. This happens when the pivot vertex's spokes in the 120-degree directions are at the boundary of the domain (undefined edges).

The current `_is_stuck` implementation is almost correct but needs one fix: it should also check that the 120-degree target positions would be **outside the hexagonal domain** (not just that the spokes are undefined). A spoke might be undefined because the domain wasn't properly set up, not because we've reached the boundary.

A more robust termination: a rhombus is settled when it is in nest B AND the next rotation (using the direction rule) would place it outside the defined domain, OR when no further rotation is possible because all target spokes are outside the domain.

---

## C. Remaining Uncertainties

### C1: The `rotate_120` vs `flip_3` question

The code has two rotation primitives: `rotate_120` (permutes all 6 spokes) and `flip_3` (permutes only the 3 alternating spokes containing the rhombus). Purbhoo's description says to "rotate the tiling 120 degrees" around the hexagon, which suggests all 6 spokes should be permuted (i.e., `rotate_120` is correct). However, there is a question of whether the 3 non-participating spokes (those not containing the rhombus or its neighbors in the 3-fold symmetry) should be left alone. If the hexagon truly has 3-fold symmetry, all 6 spokes are related by the symmetry and all should be permuted.

**Recommendation**: Use `rotate_120` (all 6 spokes), not `flip_3`. The 3-fold symmetry of the hexagon means all edges participate. Using `flip_3` would break the tile validity of the faces not directly involved.

### C2: Orientation after multiple rotations

When a rhombus undergoes several 120-degree and 180-degree rotations in sequence, the orientation flips between H, NE, and NW. The direction rule must be applied fresh at each step based on the CURRENT spoke, not based on any memory of previous orientations. The code does this correctly -- it recomputes the direction at each iteration of the while loop.

### C3: Flock ordering

Purbhoo (page 13): "For north or east migration, we begin by locating the rhombus in zeta which is maximal in the standard order." For A-to-B migration, the direction depends on whether it's north, east, south, or west. The compass direction must be consistent with the flock orientation. The code uses `flock_standard_order` which should return rhombi from maximal to minimal for north/east migration, or minimal to maximal for south/west.

This should be verified empirically but is a separate issue from the direction rule bug.

### C4: Empirical validation against Figure 2

Figure 2 shows 14 subfigures of a complete migration (3 rhombi, labels 1, 1, 2). Tracing the orientation changes visible in each subfigure would provide empirical ground truth for the direction rule. The pink (magenta) and cyan colors distinguish the two rhombus orientations. I cannot read exact edge orientations from the PDF rendering at this resolution, but the overall trajectory (left to right, roughly horizontal path) is consistent with the H/NE orientation cycling predicted by the analysis above.

---

## D. Concrete Implementation Sketch

### D1: Replace the direction heuristic in `migrate_single`

The fix is surgical: replace lines 967-987 (the `direction_score` heuristic) with a spoke-based lookup.

```python
# REPLACE this block (lines 967-987):
#   # 3-fold: rotate 120 degrees in direction that maximizes j (rightward)
#   cw_spoke = (my_spoke + 2) % 6
#   ccw_spoke = (my_spoke - 2) % 6
#   cw_pos = _spoke_to_position(...)
#   ccw_pos = _spoke_to_position(...)
#   def direction_score(p): ...
#   if direction_score(cw_pos) >= direction_score(ccw_pos): ...

# WITH this:
            # 3-fold: Purbhoo's direction rule.
            # "Rotate in whichever way causes an edge to end up horizontal or vertical."
            # For A->B (rightward) migration:
            #   H-edges are "horizontal", NE-edges are "vertical", NW-edges are "neither".
            #   Always rotate AWAY from NW.
            #   NE spokes (1, 4): use CW (+1). All others (H: 0,3 and NW: 2,5): use CCW (-1).
            _DIRECTION_FOR_SPOKE = {0: -1, 1: +1, 2: -1, 3: -1, 4: +1, 5: -1}
            direction = _DIRECTION_FOR_SPOKE[my_spoke]

            if direction == 1:
                rotate_120(hm, cr, cj, direction=1)
                new_spoke = (my_spoke + 2) % 6
            else:
                rotate_120(hm, cr, cj, direction=-1)
                new_spoke = (my_spoke - 2) % 6
```

This is equivalent to using the existing `_MIGRATION_SHIFT` table (which gives the correct answers) but translates it into the `rotate_120` convention used by `migrate_single`.

### D2: Nest-B domain preparation

Before calling `migrate_flock`, add a preparation step:

```python
def prepare_nest_b(hm, N, d):
    """Clear nest B of pre-existing rhombi and ensure all edges are defined as 1.

    Nest B is the right region where j + r > N (approximately).
    All edges in this region must be 1 (triangular tiles only) before migration.
    """
    # Remove any 10-edges in the nest-B region
    to_fix = []
    for ek, label in hm.edges.items():
        if label == 10:
            # Check if this edge is in the nest-B zone
            (r1, j1), (r2, j2) = ek
            if _is_in_nest_b_region(r1, j1, r2, j2, N):
                to_fix.append(ek)
    for ek in to_fix:
        hm.edges[ek] = 1

    # Ensure all edges in the nest-B region exist and are set to 1
    for r in range(0, N):
        for j in range(N - r, N + d + 1):
            for ek in _up_edges(r, j) + _down_edges(r, j):
                hm.edges.setdefault(ek, 1)
```

The exact boundary of the nest-B region depends on the parameters n, d and the partition shapes, but the principle is: enough room for all migrating rhombi to pack into their target positions.

### D3: Summary of changes

1. **`migrate_single` lines 967-987**: Replace `direction_score` heuristic with `_DIRECTION_FOR_SPOKE` lookup table. (~10 lines changed.)
2. **Add `prepare_nest_b` function**: Clear the nest-B region of pre-existing rhombi and fill with 1-edges. (~15 lines.)
3. **Call `prepare_nest_b` before `migrate_flock`**: One line addition.
4. **Optional**: Move `_DIRECTION_FOR_SPOKE` to module level (next to `_MIGRATION_SHIFT`) and add a comment explaining the equivalence.

### D4: Verification plan

1. Run the test case: lambda = mu = (2,1), nu = (3,2,1), N = 6. Expect 6 rhombi at 6 distinct positions.
2. Check that the 6 positions match `right_nest_positions((2,1), 3, 3)` = [(1,4),(1,5),(1,6),(3,2),(3,3),(5,0)].
3. Verify tile validity after migration (all faces satisfy DELTA/NABLA constraints).
4. Test a second case: the Purbhoo Figure 2 example (3 rhombi, labels 1,1,2) to confirm the migration path matches the paper's illustration.
5. Test commutativity: migrate A->B, then migrate B->A on the result; verify we recover the original mosaic (Theorem 2).

---

## E. Why the Bug Went Undetected for Rhombus 1

The HANDOFF.md notes that rhombus 1 migrates correctly. This is because rhombus 1 (the maximal one in standard order) is the first to travel through the interior. Before any migration has occurred, the interior mosaic is "pristine" -- it has the original tile structure from `from_puzzle()`. In this pristine interior, the position-based heuristic `direction_score` happens to agree with the orientation-based rule for the sequence of hexagons that rhombus 1 traverses.

After rhombus 1 migrates, the "wake" it leaves behind (Purbhoo's terminology, page 17) changes the tile structure along its path. Subsequent rhombi (2-6) traverse partially-modified terrain. The position-based heuristic, which depends only on coordinates and not on the actual tile structure, gives the same answers as before -- but the correct direction may have changed because the tiles have changed. The orientation-based lookup, by contrast, depends only on the spoke index (which is deterministic from the rhombus's current position relative to the pivot vertex) and is invariant to the tile modifications.

This explains why the bug manifests only for rhombi 3-6: after enough rhombi have migrated and modified the interior, the position heuristic starts disagreeing with the correct orientation rule in the ambiguous NW-spoke cases.
