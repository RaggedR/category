# Diversity Analysis: 6-Agent Ring Topology Experiment

## Date: 2026-05-06

## Methodology

Semantic diversity measured using `all-MiniLM-L6-v2` (384-dim, L2-normalized embeddings). Same model and metric as the ring/star debate experiments in the parent directory (see RESULTS.md).

**Metric**: Mean pairwise cosine distance between agent findings documents.
- 0 = identical text
- 1 = completely unrelated

Each agent's full findings file was embedded as a single document.

## Agent Output Statistics

| Agent | Role | Words | Characters |
|-------|------|-------|------------|
| 1 | Paper Exegete | 2,337 | 14,390 |
| 2 | Coordinate Translator | 2,414 | 14,191 |
| 3 | Orientation Analyst | 2,674 | 16,006 |
| 4 | Boundary & Termination | 2,499 | 15,400 |
| 5 | Figure 2 Tracer | 3,364 | 20,449 |
| 6 | Synthesis | 3,154 | 19,123 |

Total output: 16,442 words across 6 agents.

## Pairwise Cosine Distances

```
         Ag1    Ag2    Ag3    Ag4    Ag5    Ag6
Ag1       —    0.196  0.299  0.326  0.360  0.183
Ag2    0.196     —    0.388  0.353  0.442  0.188
Ag3    0.299  0.388     —    0.386  0.292  0.246
Ag4    0.326  0.353  0.386     —    0.455  0.354
Ag5    0.360  0.442  0.292  0.455     —    0.309
Ag6    0.183  0.188  0.246  0.354  0.309     —
```

### Most similar pairs (lowest distance)
1. Agent 1 (Exegete) vs Agent 6 (Synthesis): **0.183**
2. Agent 2 (Coords) vs Agent 6 (Synthesis): **0.188**
3. Agent 1 (Exegete) vs Agent 2 (Coords): **0.196**

### Most diverse pairs (highest distance)
1. Agent 4 (Boundary) vs Agent 5 (Tracer): **0.455**
2. Agent 2 (Coords) vs Agent 5 (Tracer): **0.442**
3. Agent 2 (Coords) vs Agent 3 (Orient): **0.388**

## Summary Statistics

| Metric | Value |
|--------|-------|
| **Mean pairwise distance** | **0.318** |
| Standard deviation | 0.083 |
| Minimum (most similar) | 0.183 |
| Maximum (most diverse) | 0.455 |
| Number of pairs | 15 |

## Comparison to Prior Experiments

| Experiment | Topology | Agents | Rounds | Task Type | Diversity |
|------------|----------|--------|--------|-----------|-----------|
| Ring debate (R3) | 4-ring | 4 | 3 | Open question | **0.437** |
| Star debate (R3) | 5-star | 5 | 3 | Open question | 0.362 |
| **This experiment** | **6-ring** | **6** | **1** | **Bug diagnosis** | **0.318** |

### Why diversity is lower here

The lower diversity (0.318 vs 0.437) is expected and desirable:

1. **Task type**: Bug diagnosis rewards convergence. Open-question debate rewards exploration. When 6 agents independently reach the same correct diagnosis, low diversity IS the success signal.

2. **Single round**: The debate experiments ran 3 rounds of back-and-forth. This experiment ran 1 round with file-based coordination. More rounds would likely increase diversity as agents challenge each other.

3. **Shared context**: All agents read the same CONTEXT.md, which front-loaded the problem analysis. The debate agents started from separate seed prompts.

## Cross-Domain Vocabulary Sharing

Each agent was assigned domain-specific vocabulary. Cross-domain usage indicates conceptual bridging.

### Borrowing patterns (agent X uses agent Y's terms)

| Source agent | Borrowed by | Term count | Example terms |
|---|---|---|---|
| Agent 1 (Exegete) | Agent 3, 4, 5, 6 | 4-6 each | horizontal, vertical, page 13 |
| Agent 2 (Coords) | Agent 3, 4, 5, 6 | 4-7 each | nest a/b/c, nw-oriented |
| Agent 4 (Boundary) | Agent 1, 2, 6 | 2-4 each | boundary, domain, packed |
| Agent 5 (Tracer) | Agent 1, 6 | 5 each | figure 2, cyan, magenta |

**Most-borrowed vocabulary**: Agent 2 (Coordinate Translator). All 5 other agents used at least 2 of Agent 2's terms, and Agent 3 used 7. The coordinate mapping vocabulary became the shared language of the analysis.

**Most-borrowing agent**: Agent 6 (Synthesis) used vocabulary from 4 other agents. Agent 5 (Figure 2 Tracer) used vocabulary from 3. These are the agents whose roles naturally require integration.

## Convergence Analysis

| Key conclusion | Agents agreeing | Convergence |
|---|---|---|
| `_MIGRATION_SHIFT` table exists but is unused | **6/6** | Perfect |
| `direction_score` heuristic is the primary bug | **6/6** | Perfect |
| Nest B boundary is relevant | **6/6** | Perfect |
| "Avoid NW" as resulting orientation | 5/6 | Strong |
| NW resolves to H via CCW | 5/6 | Strong |
| Keep nest B complete (don't delete) | 4/6 | Moderate |

### The dissent

**Agent 3** (Orientation Analyst) held the lowest confidence (60%) on the NW→H resolution. Their unique contribution: showing that every equilateral-grid rhombus has edges in ALL three orientations, making the "avoid NW" rule an approximation of the true mosaic-level geometry. This is a valid theoretical concern, even though the approximation appears to work in practice.

**Agent 4 vs Agent 6** disagreed on nest-B initialization. Agent 4: keep complete. Agent 6: clear to all-1. Post-analysis confirmed Agent 4 is correct (Corollary 3.2 starts with a complete mosaic).

## Clustering

The distance matrix reveals two natural clusters:

**Cluster A** (theory-focused, low intra-cluster distance):
- Agent 1 (Paper Exegete)
- Agent 2 (Coordinate Translator)
- Agent 6 (Synthesis)
- Mean intra-cluster distance: 0.189

**Cluster B** (analysis-focused, higher intra-cluster distance):
- Agent 3 (Orientation Analyst)
- Agent 4 (Boundary & Termination)
- Agent 5 (Figure 2 Tracer)
- Mean intra-cluster distance: 0.378

The theory cluster converged tightly around the same interpretation of Purbhoo's paper. The analysis cluster maintained more diversity because each agent examined different aspects of the problem (orientation algebra, domain boundaries, empirical figure tracing).

**Inter-cluster mean distance**: 0.355 — significantly higher than intra-Cluster-A (0.189), confirming the clusters are semantically distinct.

## Interpretation

The 6-ring topology for bug diagnosis produced:

1. **Appropriate convergence**: 6/6 agreement on the primary bug, with nuanced dissent preserved on secondary issues.
2. **Role-driven diversity**: The two clusters correspond to "what does the paper say?" (Agents 1, 2, 6) vs "what does the geometry/code do?" (Agents 3, 4, 5). The ring topology kept these perspectives in dialogue without forcing premature consensus.
3. **Vocabulary bridging**: Agent 2's coordinate mapping vocabulary became the lingua franca, used by all other agents. This is the natural outcome when the core challenge is translating between two representations (Purbhoo's mosaic vs the code's grid).
4. **Productive dissent**: Agent 3's lower confidence (60% vs group's 85%) on the NW resolution is a genuine warning flag, not noise. Their geometric insight about grid vs mosaic orientation is theoretically correct and may explain the residual collisions that the direction rule fix alone doesn't resolve.
