# Experimental Design

## What we are comparing

Two communication structures for multi-agent debate:

1. **Directed cycle** (ring): 4 specialist agents, each sees only its predecessor's most recent output. Communication is one-way. 4 arcs.
2. **Bidirectional star**: 5 agents (4 specialists + 1 generalist hub). Hub sees all specialists' outputs, broadcasts synthesis back. 8 arcs.

## The spiral, not the cycle

The "ring" is not a closed cycle in any single round. Agents execute sequentially:

- **Round 1**: Agent A starts with no input (chain, not cycle). A->B->C->D.
- **Round 2**: Agent A sees Agent D's round-1 output. The cycle closes for the first time. But the information is staggered: A is responding to D's *previous-round* output.
- **Round 3+**: Same staggering. Each revolution is offset by one round.

The ring is better understood as a **spiral** — a chain that wraps forward in time. The first and last rounds are incomplete: round 1 is a pure chain (no feedback), and the final round's last agent never has its output fed back.

The star does not have this problem. Each round: hub synthesises all previous-round specialist outputs, then all specialists receive the synthesis simultaneously. The star is temporally symmetric within each round.

## Known confounds

| Confound | Ring | Star | Notes |
|----------|------|------|-------|
| Arc count (bandwidth) | 4 | 8 | Star has 2x bandwidth |
| Node count | 4 | 5 | Star has an extra agent |
| Parallelism | Sequential | Parallel peripherals | Star is faster |
| Temporal structure | Spiral (chain in round 1) | Symmetric per round | Different dynamical systems |
| Hub bottleneck | None | Yes | Hub homogenises |

These confounds mean we cannot attribute diversity differences to any single factor. The comparison tests ring-vs-star *as practitioners would use them*, not any isolated topological property.

## Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Rounds | 3 | Minimum for cycle to close (round 2) + one more |
| Coupling | Negative ("challenge and extend") | Maximises diversity signal |
| Model | Claude Opus | Rich mathematical reasoning |
| Max tokens per response | 800 | Enough for substantive argument |
| RAG chunks per agent | 3 | ~4500 chars, keeps prompt under 6K tokens |
| Seeds | 1 (current), **10 (target)** | Need multiple seeds for statistical power |

## What we observe

Diversity measured as mean pairwise cosine distance (all-MiniLM-L6-v2) per round.

Current results (single run):

| Round | Ring | Star |
|-------|------|------|
| 1 | 0.431 | **0.489** |
| 2 | **0.526** | 0.448 |
| 3 | **0.437** | 0.362 |

- Star wins round 1 (higher bandwidth? more diverse starting input? temporal symmetry?)
- Ring wins rounds 2-3 (cycle has closed, local disagreements persist)
- Star converges monotonically; ring is non-monotone (peaks round 2)

This temporal pattern is consistent with GA Experiment 3, but with n=1 we cannot draw conclusions.

## What we need

1. **10+ seeds** per topology to compute confidence intervals
2. **Per-agent diversity** to understand whether the pattern is driven by specific agents
3. **Round-by-round qualitative coding** to assess whether diversity reflects genuine intellectual variety or just surface-level rephrasing
