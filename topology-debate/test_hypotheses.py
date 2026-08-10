#!/usr/bin/env python3
"""Test the four hypotheses from the Opus ring-topology debate.

Hypotheses concern the open question: does the KTW puzzle transfer matrix
define a growth diagram whose face labels biject puzzles with LR tableaux?

Test cases:
  H1-H4 on c^(3,2,1)_{(2,1),(2,1)} = 2  (12 free ∇, 3 tableau entries)
  H1-H4 on c^(4,4,3,1)_{(3,2,1),(3,2,1)} = 3  (22 free ∇, 6 tableau entries)
"""

import sys
from pathlib import Path
from itertools import combinations

# Import puzzle code
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "puzzles"))
from lr_coefficients import (
    lr_puzzles, lr_count, _DELTA_FUNC, _NABLA_BRANCH, _NABLA_FUNC,
    partition_to_boundary, boundary_to_partition, _bounding_box,
)


# -- Utilities ----------------------------------------------------------------

def extract_nabla_info(puzzle):
    """Extract all ∇ tile positions, their inputs, and whether they were free.

    Returns a list of dicts, one per ∇ tile, ordered row-by-row left-to-right:
        {
            'row': r, 'col': j,
            'left_diag': left╲ input,
            'top': actual top─ chosen,
            'right_diag': actual right╱ chosen,
            'free': bool (True if 2 choices were available),
            'choice': 0 or 1 (index into _NABLA_BRANCH[left_diag]),
        }
    """
    N = puzzle.N
    levels = puzzle.levels
    nablas = []

    for r in range(N):
        w = N - r  # number of △ tiles in this row
        state = levels[r]   # bottom edges (length w)
        nxt = levels[r + 1] if r + 1 <= N else ()  # top edges (length w-1)

        # Walk left-to-right through the row, tracking the diagonal
        # The left boundary diagonal is the ╱ edge entering from the left
        # We need to reconstruct it from boundary conditions
        lam_bits = partition_to_boundary(puzzle.lam, puzzle.n, puzzle.m)
        mu_bits = partition_to_boundary(puzzle.mu, puzzle.n, puzzle.m)

        diag = lam_bits[r]  # left boundary ╱ for row r

        for j in range(w):
            # △_j: deterministic
            key = (state[j], diag)
            if key not in _DELTA_FUNC:
                break
            r_diag = _DELTA_FUNC[key]  # right╲ output of △_j

            if j < w - 1:
                # ∇_j: r_diag is left╲ input
                left_diag = r_diag
                branches = _NABLA_BRANCH.get(left_diag, ())
                is_free = len(branches) > 1

                # Actual choice: look at the top edge that was placed
                actual_top = nxt[j]
                actual_right = None
                choice_idx = None

                for idx, (t_val, n_diag) in enumerate(branches):
                    if t_val == actual_top:
                        actual_right = n_diag
                        choice_idx = idx
                        break

                nablas.append({
                    'row': r,
                    'col': j,
                    'left_diag': left_diag,
                    'top': actual_top,
                    'right_diag': actual_right,
                    'free': is_free,
                    'choice': choice_idx,
                })

                # Next diagonal input is right╱ of this ∇
                diag = actual_right
            else:
                # Last △ in row — no ∇ follows
                pass

    return nablas


def get_choice_vector(puzzle):
    """Get the binary choice vector for free ∇ tiles only."""
    nablas = extract_nabla_info(puzzle)
    return tuple(n['choice'] for n in nablas if n['free'])


def get_all_nabla_positions(puzzle):
    """Get positions of all free ∇ tiles."""
    nablas = extract_nabla_info(puzzle)
    return [(n['row'], n['col']) for n in nablas if n['free']]


# -- Hypothesis Tests ---------------------------------------------------------

def test_h1_active_passive(lam, mu, nu):
    """H1: Active ∇ tiles (those affecting downstream) == |ν/λ| in count.

    Method: For each valid puzzle, for each free ∇, flip the choice and
    check if subsequent forced △ tiles produce different outputs.
    """
    print(f"\n{'='*60}")
    print(f"H1: Active/Passive Decomposition")
    print(f"    c^{nu}_{{{lam},{mu}}}")
    print(f"{'='*60}")

    puzzles = lr_puzzles(lam, mu, nu)
    c = len(puzzles)

    # Compute skew shape size
    nu_padded = nu + (0,) * (len(lam) if len(lam) > len(nu) else 0)
    lam_padded = lam + (0,) * (max(len(nu) - len(lam), 0))
    skew_size = sum(nu_padded[i] - lam_padded[i] for i in range(min(len(nu_padded), len(lam_padded))))
    # More simply: |ν| - |λ|
    skew_size = sum(nu) - sum(lam)

    print(f"  c = {c}")
    print(f"  |ν/λ| = |{nu}| - |{lam}| = {sum(nu)} - {sum(lam)} = {skew_size}")

    for pi, puzzle in enumerate(puzzles):
        nablas = extract_nabla_info(puzzle)
        free_nablas = [n for n in nablas if n['free']]
        total_nablas = len(nablas)

        print(f"\n  Puzzle {pi+1}:")
        print(f"    Total ∇ tiles: {total_nablas}")
        print(f"    Free ∇ tiles:  {len(free_nablas)}")

        # For each free ∇, check if it's "active" by examining whether
        # the two branches lead to different downstream states
        active_positions = []
        for fn in free_nablas:
            branches = _NABLA_BRANCH[fn['left_diag']]
            # Both branches give different (top, right╱) pairs
            # A ∇ is "active" if the two branches produce different top values
            # (which then affect the next row's state)
            tops = [b[0] for b in branches]
            if len(set(tops)) > 1:
                active_positions.append((fn['row'], fn['col']))

        print(f"    Active ∇ tiles (different top values): {len(active_positions)}")
        print(f"    Positions: {active_positions}")

        # Also check: do both branches always give different top values?
        for fn in free_nablas:
            branches = _NABLA_BRANCH[fn['left_diag']]
            tops = [b[0] for b in branches]
            right_diags = [b[1] for b in branches]
            is_active = len(set(tops)) > 1
            print(f"      ∇({fn['row']},{fn['col']}): left╲={fn['left_diag']} "
                  f"→ branches={branches} → active={is_active}")

    # H1 prediction: number of active ∇ = |ν/λ|
    print(f"\n  H1 PREDICTION: active ∇ count == {skew_size}")
    print(f"  H1 RESULT: See above — compare active count to {skew_size}")


def test_h3_exhaustive_search(lam, mu, nu):
    """H3: Enumerate all 2^F choice vectors, count valid puzzles.

    This doesn't literally try all 2^F vectors (the tree is path-dependent).
    Instead, we modify the enumeration to count ALL branching paths (valid or not)
    and compare to the number of valid puzzles.
    """
    print(f"\n{'='*60}")
    print(f"H3: Exhaustive Search / Bulk Invariance")
    print(f"    c^{nu}_{{{lam},{mu}}}")
    print(f"{'='*60}")

    puzzles = lr_puzzles(lam, mu, nu)
    c = len(puzzles)

    # Count free ∇ per puzzle
    for pi, puzzle in enumerate(puzzles):
        cv = get_choice_vector(puzzle)
        positions = get_all_nabla_positions(puzzle)
        print(f"  Puzzle {pi+1}: choice vector = {cv} (len={len(cv)})")
        print(f"    Free ∇ positions: {positions}")

    # Check if all puzzles have the same number of free ∇
    free_counts = [len(get_choice_vector(p)) for p in puzzles]
    all_same = len(set(free_counts)) == 1
    F = free_counts[0] if all_same else None
    print(f"\n  Free ∇ count per puzzle: {free_counts}")
    print(f"  All same? {all_same}")

    if F is not None:
        print(f"  F = {F}")
        print(f"  2^F = {2**F}")
        print(f"  Valid puzzles (c) = {c}")
        print(f"  Valid/Total ratio = {c}/{2**F} = {c/2**F:.6f}")

    # Compute Hamming distances between choice vectors
    if c > 1:
        vectors = [get_choice_vector(p) for p in puzzles]
        print(f"\n  Hamming distances between valid choice vectors:")
        for i, j in combinations(range(c), 2):
            dist = sum(a != b for a, b in zip(vectors[i], vectors[j]))
            print(f"    Puzzle {i+1} vs {j+1}: distance = {dist}")
            # Show which positions differ
            diffs = [k for k in range(len(vectors[i])) if vectors[i][k] != vectors[j][k]]
            print(f"    Differing positions: {diffs}")

    print(f"\n  H3 PREDICTION: Exactly c={c} of 2^{F} vectors are valid")
    print(f"  H3 NOTE: The tree is path-dependent, so not all 2^{F} vectors")
    print(f"           are reachable. Need deeper analysis to check if")
    print(f"           every reachable vector is valid.")


def test_h2_geometric(lam, mu, nu):
    """H2: Active ∇ positions border the skew shape ν/λ."""
    print(f"\n{'='*60}")
    print(f"H2: Geometric Identification (Top Faces)")
    print(f"    c^{nu}_{{{lam},{mu}}}")
    print(f"{'='*60}")

    puzzles = lr_puzzles(lam, mu, nu)

    # Compute skew shape boxes
    skew_boxes = []
    for i in range(len(nu)):
        lam_i = lam[i] if i < len(lam) else 0
        for j in range(lam_i, nu[i]):
            skew_boxes.append((i, j))

    print(f"  Skew shape ν/λ = {nu}/{lam}")
    print(f"  Skew boxes: {skew_boxes}")
    print(f"  |ν/λ| = {len(skew_boxes)}")

    # For each puzzle, show all free ∇ positions and check for geometric pattern
    for pi, puzzle in enumerate(puzzles):
        nablas = extract_nabla_info(puzzle)
        free_positions = [(n['row'], n['col']) for n in nablas if n['free']]
        print(f"\n  Puzzle {pi+1}:")
        print(f"    Free ∇ positions (row, col in triangle): {free_positions}")
        print(f"    Count: {len(free_positions)}")

    # Check if the free positions are the same across all puzzles
    all_free = [set(get_all_nabla_positions(p)) for p in puzzles]
    same_positions = all(fp == all_free[0] for fp in all_free)
    print(f"\n  Same free positions across all puzzles? {same_positions}")
    if same_positions:
        print(f"  Common free ∇ positions: {sorted(all_free[0])}")

    print(f"\n  H2 PREDICTION: Active ∇ positions correspond geometrically")
    print(f"                 to the skew boxes {skew_boxes}")
    print(f"  H2 NOTE: The coordinate mapping between triangle grid and")
    print(f"           partition coordinates needs to be worked out to")
    print(f"           verify this. For now, compare counts.")


def test_h4_code_structure(lam, mu, nu):
    """H4: Error-correcting code structure."""
    print(f"\n{'='*60}")
    print(f"H4: Error-Correcting Code Structure")
    print(f"    c^{nu}_{{{lam},{mu}}}")
    print(f"{'='*60}")

    puzzles = lr_puzzles(lam, mu, nu)
    c = len(puzzles)
    skew_size = sum(nu) - sum(lam)

    vectors = [get_choice_vector(p) for p in puzzles]
    F = len(vectors[0]) if vectors else 0

    print(f"  Claimed code: ({F}, {skew_size}) binary code")
    print(f"  Codewords (valid puzzles): {c}")
    print(f"  Information bits: {skew_size}")
    print(f"  Parity bits: {F - skew_size}")

    # Check if differing positions between codewords are consistent
    if c > 1:
        # Find positions where ANY two codewords differ
        all_diff_positions = set()
        for i, j in combinations(range(c), 2):
            diffs = {k for k in range(F) if vectors[i][k] != vectors[j][k]}
            all_diff_positions |= diffs
            print(f"  Puzzles {i+1} vs {j+1}: differ at positions {sorted(diffs)}")

        print(f"\n  Total positions that ever differ: {sorted(all_diff_positions)}")
        print(f"  Count of distinguishing positions: {len(all_diff_positions)}")
        print(f"  Expected (|ν/λ|): {skew_size}")

        # H4 prediction: distinguishing positions == information bits
        if len(all_diff_positions) == skew_size:
            print(f"\n  *** H4 CONFIRMED: exactly {skew_size} positions distinguish "
                  f"the {c} codewords ***")

            # Extract just the information bits
            info_positions = sorted(all_diff_positions)
            print(f"  Information bit positions: {info_positions}")
            for pi, v in enumerate(vectors):
                info_bits = tuple(v[k] for k in info_positions)
                print(f"    Puzzle {pi+1}: info bits = {info_bits}")
        else:
            print(f"\n  H4 PARTIAL: {len(all_diff_positions)} distinguishing "
                  f"positions != {skew_size} expected")

    # Check if non-distinguishing positions are identical across all codewords
    if c > 1 and F > 0:
        parity_positions = [k for k in range(F) if k not in all_diff_positions]
        parity_values = [tuple(v[k] for k in parity_positions) for v in vectors]
        all_parity_same = len(set(parity_values)) == 1
        print(f"\n  Parity positions ({len(parity_positions)}): {parity_positions}")
        print(f"  Parity values identical across all puzzles? {all_parity_same}")
        if all_parity_same:
            print(f"  Parity value: {parity_values[0]}")


# -- Main --------------------------------------------------------------------

def run_all_tests(lam, mu, nu):
    """Run all four hypothesis tests on one triple."""
    print(f"\n{'#'*60}")
    print(f"# Testing c^{nu}_{{{lam},{mu}}}")
    print(f"# c = {lr_count(lam, mu, nu)}")
    print(f"{'#'*60}")

    test_h1_active_passive(lam, mu, nu)
    test_h2_geometric(lam, mu, nu)
    test_h3_exhaustive_search(lam, mu, nu)
    test_h4_code_structure(lam, mu, nu)


if __name__ == "__main__":
    # Minimal test case: c = 2
    run_all_tests((2, 1), (2, 1), (3, 2, 1))

    # Second test case: c = 3
    run_all_tests((3, 2, 1), (3, 2, 1), (4, 4, 3, 1))
