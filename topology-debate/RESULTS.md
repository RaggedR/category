# Topology-Debate Experiment Results

## The Problem

We have an open question in algebraic combinatorics: does the KTW puzzle transfer matrix define a growth diagram whose face labels biject puzzles with Littlewood-Richardson tableaux? The question sits at the intersection of four specialisms — Fomin growth diagrams, Yang-Baxter integrability, Purbhoo's migration bijection, and computational enumeration — none of which alone has cracked it.

We use this problem as a testbed for **multi-agent debate**, where specialist LLM agents discuss an open mathematical problem under different communication topologies. The question: does the topology of information flow between agents affect the diversity and quality of ideas produced?

## Motivation: The GECCO Paper

This experiment is a direct extension of:

> Langer, R. and Turing, C. (2026). *Composition Determines Diversity: Categorical Fingerprints of Genetic Algorithms.* GECCO 2026.

That paper establishes a **universal topology ordering** for genetic algorithm island models:

```
none > ring > star > random > fully connected
```

This ordering holds with perfect rank concordance (Kendall's W = 1.0, p = 0.00008) across six unrelated domains. Topology explains 23.9× more variance than domain. The key finding: **migration topology determines diversity dynamics independently of what is being optimised.**

We test whether this ordering transfers from genetic algorithms to LLM multi-agent debate. If it does, it confirms the categorical structure is truly domain-independent — it governs information flow regardless of whether the "individuals" are bitstrings or language model outputs.

## Experimental Setup

### Model and Configuration

| Parameter | Value |
|-----------|-------|
| Model | Claude Opus (`claude-opus-4-20250514`) |
| Max tokens | 800 per response |
| Tools | None (text-only generation) |
| Coupling | Negative ("challenge and extend") |
| Rounds | 3 |
| RAG | 3 chunks per agent from puzzle-rag (ChromaDB, 721 chunks from ~15 papers) |

### Why Opus, No Tools

- **Opus** produces richer mathematical reasoning than Sonnet. For a diversity comparison, we need responses with enough substance to differentiate semantically.
- **No tools** eliminates a confounding variable. Tool use causes convergence (agents seeing the same computational results agree more). We want to isolate the effect of topology on diversity of *reasoning*, not diversity of *data access*.

### Cost and Duration

| Topology | Duration | Responses | Approx. Cost |
|----------|----------|-----------|-------------|
| Ring | 266s | 12 | ~$0.45 |
| Star | 134s | 15 | ~$0.55 |

The star is faster despite more responses because its 4 peripheral agents run **in parallel** (fan-out from center), while the ring is fully sequential.

## Topologies

### Ring (Cycle) — 4 nodes, 4 edges

```
Growth Diagram → Yang-Baxter → Bijection → Computational → (loop)
```

Each agent sees **only its predecessor's** most recent output. Information propagates sequentially around the cycle. In each round, all 4 agents speak once in order.

The ring ordering maximises productive tension: each agent's predecessor represents a complementary (not similar) perspective.

### Star — 5 nodes, 4 edges

```
         Growth Diagram
              ↑↓
Yang-Baxter ←→ VANILLA ←→ Bijection
              ↑↓
         Computational
```

A **vanilla agent** (generalist mathematician, no domain specialism) acts as the central hub. It sees all 4 peripheral specialists' outputs; each specialist sees only the hub's synthesis.

The star has the same edge count as the ring (4 edges) but fundamentally different information flow: all cross-pollination must pass through the center. Peripherals never see each other directly.

### Why a 5th Agent for the Star

A star on 4 nodes has only 3 edges (one center connected to 3 peripherals). To match the ring's 4 edges, we need either:
- 5 nodes with a 4-edge star (our choice), or
- 4 nodes with 3 star edges + 1 peripheral-to-peripheral edge (no longer a pure star)

The vanilla hub also avoids the confound of giving one specialist a structural advantage. In a 4-node star with a specialist as center, that specialist's framing would dominate.

## Diversity Metric

### Definition

We measure **semantic diversity** as mean pairwise cosine distance between agent responses within each round:

1. Embed each response using `all-MiniLM-L6-v2` (384-dim, L2-normalised)
2. Compute pairwise cosine similarity: `sim(a,b) = embed(a) · embed(b)`
3. Convert to distance: `dist(a,b) = 1 - sim(a,b)`
4. Report mean distance across all pairs in the round

**Interpretation**: Higher distance = more diverse perspectives. A score of 0 means all agents said the same thing; a score near 1 means responses are semantically unrelated.

We also track **cross-domain connections**: when agent X uses vocabulary from agent Y's domain (a heuristic for conceptual bridging). This is keyword-based, not embedding-based.

### Why This Metric

This is the same metric used in `~/git/orchestration/experiment.py` and validated against the GECCO paper's diversity measurements. Using `all-MiniLM-L6-v2` gives fast, cheap embeddings adequate for relative comparison. We don't need a math-specific embedding model — we're measuring diversity of *expression*, not mathematical similarity.

## Results

### Diversity Over Time

| Round | Ring (cycle) | Star | Δ (ring − star) |
|-------|-------------|------|-----------------|
| 1 | 0.431 | **0.489** | −0.058 |
| 2 | **0.526** | 0.448 | +0.078 |
| 3 | **0.437** | 0.362 | +0.075 |

### Cross-Domain Connections

| Round | Ring | Star |
|-------|------|------|
| 1 | 9 | 6 |
| 2 | 5 | 6 |
| 3 | 7 | 7 |

### Key Observations

1. **The star converges monotonically** (0.489 → 0.448 → 0.362). Diversity decreases every round. The central synthesizer homogenises perspectives by broadcasting a unified framing that all peripherals respond to.

2. **The ring has non-monotone diversity** (0.431 → 0.526 → 0.437). Diversity *increases* in round 2 before declining. The sequential information flow allows local disagreements to amplify before eventually damping.

3. **Ring maintains higher diversity from round 2 onward.** By round 3, the ring has 21% more diversity than the star (0.437 vs 0.362).

4. **Star has higher initial diversity** (round 1: 0.489 vs 0.431). This is because in round 1 the star's peripherals all respond to the same generic prompt from the vanilla center, while ring agents each respond to a different predecessor. The ring's initial diversity is constrained by the sequential chain within round 1 itself.

5. **Cross-domain connections are similar.** The star's vanilla agent accounts for many connections (it naturally borrows from all domains since it sees all outputs). Specialist-to-specialist bridging is lower in the star because specialists never see each other.

## Confirmation of the GECCO Result

The GECCO paper predicts: **ring > star** for diversity maintenance.

Our results confirm this for LLM multi-agent debate:

| | GECCO (GA islands) | This experiment (LLM debate) |
|---|---|---|
| ring > star? | Yes (W = 1.0) | Yes (rounds 2-3) |
| Star converges faster? | Yes | Yes (monotonic decline) |
| Ring sustains diversity? | Yes | Yes (non-monotone, higher floor) |

The mechanism is the same in both settings: **the star's hub acts as an information bottleneck that smooths out productive disagreements.** In GAs, the star topology allows the best individual to spread everywhere via the hub, collapsing diversity. In LLM debate, the vanilla synthesizer broadcasts a unified framing that pulls specialists toward consensus.

The ring's sequential propagation preserves local variation — each agent only hears one voice, so disagreements between non-adjacent agents persist longer.

## LangGraph Implementation Details

### Framework

Built on [LangGraph](https://github.com/langchain-ai/langgraph) (not raw API calls) for:

1. **Compilable graph objects** — the topology is visible in code as explicit edges, not implicit in context-passing logic
2. **Typed state** via `TypedDict` — clear schema for what flows through the graph
3. **Built-in cycle support** — the round loop is a conditional edge back to the first node
4. **Concurrent fan-out** — the star topology's parallel peripheral execution comes free from LangGraph's graph semantics

### State Design

```python
class DebateState(TypedDict):
    problem_statement: str
    max_rounds: int
    coupling: str              # "negative" or "positive"
    use_tools: bool
    current_round: int
    current_agent_idx: int
    responses: Annotated[list[dict], operator.add]  # fan-in compatible
    metrics: Annotated[list[dict], operator.add]
```

The `Annotated[list, operator.add]` reducer allows multiple nodes writing in parallel (fan-in from star peripherals) to merge their outputs automatically.

### Ring Graph Structure

```
START → agent_growth_diagram → agent_yang_baxter → agent_bijection → agent_computational
                                                                           ↓
                                                                    [should_continue]
                                                                     ↙            ↘
                                                              round_transition     END
                                                                     ↓
                                                              agent_growth_diagram (loop)
```

### Star Graph Structure

```
START → agent_vanilla → [agent_growth_diagram, agent_yang_baxter,
                          agent_bijection, agent_computational]  (parallel)
                                    ↓ (all)
                                  join
                                    ↓
                             [should_continue]
                              ↙            ↘
                       round_transition     END
                              ↓
                        agent_vanilla (loop)
```

### RAG Integration

Each agent queries a ChromaDB vector store (721 chunks from ~15 papers on puzzles, integrability, and growth diagrams). Queries adapt each round: the seed query is augmented with a snippet of the predecessor's output, so agents retrieve different literature as the debate evolves.

### API Integration

Direct Anthropic SDK (`anthropic` Python package), not the LangChain-Anthropic adapter. The SDK gives full control over the message format and avoids an unnecessary abstraction layer. LangGraph is used only for graph orchestration, not for LLM calls.

## Files

| File | Purpose |
|------|---------|
| `opus_debate.json` | Ring topology results (4 agents × 3 rounds = 12 responses) |
| `star_debate.json` | Star topology results (5 agents × 3 rounds = 15 responses) |
| `state.py` | Agent definitions, vanilla agent, state schema |
| `agents.py` | Ring agent nodes (sequential) |
| `agents_star.py` | Star agent nodes (center + peripheral) |
| `graph.py` | Ring topology graph construction |
| `graph_star.py` | Star topology graph construction |
| `evaluate.py` | DiversityMeter + cross-domain connection detection |
| `main.py` | Ring experiment CLI |
| `run_star.py` | Star experiment CLI |

## Future Work

1. **Additional topologies**: fully connected (all-to-all), none (isolated agents), random graph — to test the full GECCO ordering
2. **Positive coupling comparison**: does "build on" (positive) vs "challenge" (negative) interact with topology?
3. **Tool-augmented star**: does giving agents computational tools change the topology effect? (Hypothesis: tools cause convergence in both topologies, narrowing the gap)
4. **Statistical significance**: run each topology 10+ times with different random seeds to compute confidence intervals
5. **LLM-as-judge**: have Claude assess whether each round produced genuine mathematical insight, not just diverse text
