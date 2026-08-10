# Handoff: Migration Algorithm — 120° Direction Bug

## What was accomplished this session

1. **Star topology experiment** — ran successfully, results in `star_debate.json`. Confirms GECCO paper's ring > star diversity ordering for LLM debate.
2. **Writeups** — `RESULTS.md` (full experiment comparison), `FINDINGS.md` (all tested/untested hypotheses), `linkedin-post.md`.
3. **Migration algorithm** — significant progress on `~/git/puzzles/hex_mosaic.py`. The algorithm now correctly migrates rhombus 1 through the interior and parks it at nest B position (1,4). Remaining: 2 collisions out of 6 rhombi.

## The remaining bug

**File**: `~/git/puzzles/hex_mosaic.py`, function `migrate_single` (around line 750)

**Symptom**: Rhombi 3-6 sometimes land at the same position as an earlier rhombus (collision). Expected: all 6 rhombi land at 6 distinct nest-B positions.

**Root cause**: The 120° rotation direction choice uses a heuristic:
```python
cw_score = direction_score(cw_pos)
ccw_score = direction_score(ccw_pos)
```
where `direction_score(p) = (1 if r>=0 else 0, j, r)`. When both options have the same j-value (rightward extent), the tiebreaker (r) makes an arbitrary choice that's sometimes wrong.

**What Purbhoo's paper actually says** (Section 3.1, arXiv:0705.1184): "Rotate 120° in the direction that makes the rhombus's shared edge change orientation toward the migration direction." For A→B (south) migration, this means the orientation should cycle: NW → H → NE → NW, advancing rightward at each step.

**The fix**: Instead of scoring by position, score by **resulting orientation**. For A→B migration:
- If current orient is NW: prefer the rotation that gives H (not NE)
- If current orient is H: prefer the rotation that gives NE (not NW)  
- If current orient is NE: prefer the rotation that gives NW (not H)

This deterministic cycle (NW→H→NE→NW) is Purbhoo's rule. The "rightward" heuristic approximates it but fails on ties.

Actually — re-reading the code: `_spoke_to_orient` maps spoke indices to orientations. After a CW rotation: spoke moves +2, after CCW: spoke moves -2. The orientation cycle depends on which spoke we start from. So:
- my_spoke=2 (NW): CW→spoke 4 (NE), CCW→spoke 0 (H). For rightward migration, prefer CCW (→H).
- my_spoke=3 (H): CW→spoke 5 (NW), CCW→spoke 1 (NE). Prefer CCW (→NE).  
- my_spoke=4 (NE): CW→spoke 0 (H), CCW→spoke 2 (NW). Prefer CW (→H)?? 

Wait — this doesn't give a clean cycle. Need to check against the paper more carefully. The spoke numbering and orientation mapping needs to be verified.

**Alternative fix**: The simpler approach — just check if the resulting orientation advances the cycle NW→H→NE rightward. Implement `_preferred_direction(my_spoke)` that returns +1 (CW) or -1 (CCW) based on the migration direction.

## How to reproduce

```bash
cd ~/git/puzzles
python -c "
from hex_mosaic import *
from hex_mosaic import _up_edges, _down_edges
from lr_coefficients import lr_puzzles

lam, mu, nu = (2, 1), (2, 1), (3, 2, 1)
N = 6
puzzles = lr_puzzles(lam, mu, nu)

# Prepare mosaic: remove nest-B, fill left/bottom nests
hm = HexMosaic.from_puzzle(puzzles[0])
for ek in [ek for ek in hm.edges if any(rv >= 0 and jv > N-rv for rv,jv in ek)]:
    del hm.edges[ek]
for r in range(-3, N+1):
    for j in range(-4, 0):
        for ek in _up_edges(r,j) + _down_edges(r,j):
            hm.edges.setdefault(ek, 1)
for r in range(-3, 0):
    for j in range(-1, N+2):
        for ek in _up_edges(r,j) + _down_edges(r,j):
            hm.edges.setdefault(ek, 1)

flock = build_flock_left(hm)
for pos, label in flock_standard_order(flock):
    new_pos, orient, _ = migrate_single(hm, pos, 'nw', label, verbose=True)
    print(f'label={label}: {pos} -> {new_pos} ({orient})')
"
```

Expected: 6 distinct positions (matching `right_nest_positions((2,1), 3, 3)` = [(1,4),(1,5),(1,6),(3,2),(3,3),(5,0)])

## What NOT to do

- Don't fill nest B with 1-edges — that causes overshooting
- Don't fill nest B with pre-placed rhombi — that causes 180° swap cascades
- Don't use position-based heuristics for 120° direction — use orientation cycling

## Context

- The migration is implementing Purbhoo arXiv:0705.1184, Section 3.1
- This is the only known bijection: KTW puzzle → LR tableau
- Successfully implementing it would test hypothesis S1 from `FINDINGS.md`
- The `~/git/puzzles/mosaic.py` file has the basic Mosaic class; `hex_mosaic.py` has the extended hexagonal version with nests and migration
