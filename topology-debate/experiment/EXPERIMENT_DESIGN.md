# Experiment: 6-Agent Ring Topology for Mathematical Problem Solving

## Date: 2026-05-06

## Motivation

The GECCO paper (Langer & Turing, 2026) establishes that ring topology preserves diversity better than star topology in genetic algorithm island models (Kendall's W = 1.0). The topology-debate experiments in this directory confirmed the same ordering for LLM multi-agent debate: ring maintained 21% more diversity than star by round 3.

This experiment extends the ring topology approach to **mathematical problem solving** — specifically, fixing a bug in the Purbhoo mosaic migration bijection (arXiv:0705.1184). Rather than debating an open question, the 6 agents collaborate on diagnosing and fixing a concrete implementation bug.

## The Problem

File: `~/git/puzzles/hex_mosaic.py`, function `migrate_single` (line ~915)

The migration algorithm moves labeled rhombi from nest A to nest B through a hexagonal mosaic. At each step, when encountering a hexagon with 3-fold rotational symmetry, the algorithm must choose between clockwise (CW) and counterclockwise (CCW) 120-degree rotation. Purbhoo's paper states the rule: "rotate in whichever way causes some edge of the rhombus to end up exactly horizontal, or exactly vertical."

The current implementation uses a position-based heuristic (`direction_score`) that maximizes the j-coordinate (rightward position). This fails on ties, causing rhombi 3-6 to collide.

## Experimental Design

### Topology: 6-Node Ring

```
Agent 1 (Paper Exegete)
    ↓ writes agent1_findings.md
Agent 2 (Coordinate Translator) reads agent1
    ↓ writes agent2_findings.md
Agent 3 (Orientation Analyst) reads agent2
    ↓ writes agent3_findings.md
Agent 4 (Boundary & Termination) reads agent3
    ↓ writes agent4_findings.md
Agent 5 (Figure 2 Tracer) reads agent4
    ↓ writes agent5_findings.md
Agent 6 (Synthesis) reads agent5
    ↓ writes agent6_findings.md → Agent 1 reads agent6
```

Each agent writes partial results to its file, which the next agent in the ring can read. The circular dependency (Agent 6 → Agent 1) creates a feedback loop.

### Agent Roles

| Agent | Role | Analytical Angle |
|-------|------|-----------------|
| 1 | **Paper Exegete** | Parse Section 3.1's exact direction rule; what do "horizontal" and "vertical" mean in Purbhoo's relative coordinate system? |
| 2 | **Coordinate Translator** | Map between Purbhoo's mosaic geometry (hexagonal tiles, 150-degree nests) and hex_mosaic.py's equilateral grid (spoke labels, face coordinates) |
| 3 | **Orientation Analyst** | Analyze the orientation cycle NW→H→NE; determine CW vs CCW for each spoke index 0-5 |
| 4 | **Boundary & Termination** | How should nest B be prepared? What is the correct stopping condition? |
| 5 | **Figure 2 Tracer** | Empirically trace the worked example to extract direction choices from the paper's diagrams |
| 6 | **Synthesis** | Read all findings, propose a unified fix |

### Model and Configuration

| Parameter | Value |
|-----------|-------|
| Model | Claude Opus 4.6 (1M context) |
| Execution | All 6 agents launched in parallel (background) |
| Communication | File-based, ring topology |
| Task type | Mathematical analysis only (no coding) |
| Shared context | CONTEXT.md with extracted problem statement, code snippets, and spoke-orientation tables |

### Shared Input Materials

All agents read:
1. `experiment/CONTEXT.md` — Problem statement, current code analysis, direction rule table
2. `HANDOFF.md` — Bug description and reproduction steps
3. `FINDINGS.md` — Prior hypotheses S1, S5 about migration
4. Purbhoo paper (PDF): `~/git/puzzles/seed-papers/purbhoo-2007-puzzles-tableaux-mosaics.pdf`
5. Source code: `~/git/puzzles/hex_mosaic.py` (relevant sections)

### Why This Design

1. **Ring over star**: Preserves diversity of analytical perspectives. No central hub to homogenize.
2. **6 agents (not 4)**: The problem has 6 distinct analytical dimensions (paper parsing, coordinate mapping, orientation analysis, boundary conditions, empirical tracing, synthesis).
3. **Specialist roles**: Each agent focuses on one analytical angle, preventing the "jack of all trades" diffusion seen in star topologies.
4. **File-based communication**: Asynchronous, persistent. Later agents can incorporate earlier agents' insights if timing permits.
5. **Parallel launch**: All agents start simultaneously. The ring dependency is "soft" — agents write their own analysis first, then check for predecessor findings.

### Differences from Prior Experiments

| Aspect | Ring Debate (Experiment 1) | This Experiment |
|--------|--------------------------|----------------|
| Goal | Generate diverse hypotheses about an open question | Diagnose and fix a specific bug |
| Agents | 4 domain specialists (Growth Diagram, Yang-Baxter, Bijection, Computational) | 6 analytical roles (Exegete, Translator, Analyst, Boundary, Tracer, Synthesis) |
| Rounds | 3 rounds of sequential responses | 1 round, parallel with file-based coordination |
| Coupling | Negative ("challenge and extend") | Collaborative (build on predecessor's findings) |
| Output | Debate transcript (JSON) | Individual findings files (Markdown) |
| Tools | None (text-only) | File read/write only |
| Evaluation | Semantic diversity metric | Convergence/divergence analysis |

## Execution Timeline

All 6 agents launched simultaneously at ~same time. Completion order:
1. Agent 1 (Paper Exegete) — ~3 min
2. Agent 3 (Orientation Analyst) — ~3 min
3. Agent 2 (Coordinate Translator) — ~4 min
4. Agent 4 (Boundary & Termination) — ~4 min
5. Agent 6 (Synthesis) — ~4 min
6. Agent 5 (Figure 2 Tracer) — ~6 min

Agent 5 took longest due to the detailed step-by-step tracing of Figure 2's 14 subfigures.

## Results

See RESULTS.md in this directory for the full analysis of findings, convergences, and disagreements.
