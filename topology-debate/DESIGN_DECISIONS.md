# Design Decisions

## Architecture

### Why LangGraph (not raw API calls like orchestration/)

The `~/git/orchestration/` project uses raw OpenAI API calls in a single 540-line file. That worked for measuring diversity, but this project has a second purpose: **interview demonstration**. LangGraph gives us:

1. A compilable graph object we can visualize and explain
2. Typed state with explicit edges — the topology is visible in code, not implicit in a `get_context()` function
3. Built-in support for cycles (the round loop) without manual counter logic
4. LangSmith tracing for free if we add an env var later

Trade-off: heavier dependency, slightly more boilerplate. Worth it for demonstrability.

### Why sequential within a round (not parallel)

A ring topology is inherently sequential: agent B sees agent A's output, agent C sees B's output. Parallel execution would mean agents only see the *previous round's* outputs, turning the ring into something more like a synchronous cellular automaton. We chose sequential because:

1. It matches the ring semantics exactly — information propagates around the ring within a single round
2. The orchestration project used sequential too (`solutions[:agent_idx]`)
3. It's easier to reason about and debug

Alternative we rejected: parallel agents with round-boundary sharing. This would be the "synchronous island model" from the GA paper. Could add as a separate topology mode later.

### Why 4 agents (not 3 or 6)

The open problem in OPEN_QUESTION.md identifies exactly 4 perspectives:
1. Growth diagrams (Fomin/Langer)
2. Yang-Baxter / integrability (Zinn-Justin)
3. Bijection theory (Purbhoo)
4. Computational verification

A 4-node ring also has nice properties: each agent is exactly 2 hops from its opposite, creating maximal "distance" between complementary perspectives (Growth Diagram ↔ Computational, Yang-Baxter ↔ Bijection).

### Why Claude (not OpenAI)

1. The mathematical reasoning quality is noticeably better on Claude for algebraic combinatorics
2. Robin's interview is for an AI engineering role — showing comfort with the Anthropic SDK is relevant
3. The `anthropic` SDK was already installed

Trade-off: OAuth tokens from Max plan don't work for scripts (Anthropic disabled this). Must use prepaid API credits (~$0.02/round Sonnet, $0.15/round Opus).

## RAG Integration

### Why query puzzle-rag each round (not pre-load context once)

Each agent's RAG query adapts based on its predecessor's output. Round 1 uses a seed query; later rounds append keywords from the predecessor. This means:

1. Agents retrieve *different* literature as the debate evolves
2. The system is genuinely adaptive, not just prompting with a fixed context window
3. It demonstrates dynamic RAG in an interview setting

Trade-off: 4 extra OpenAI embedding API calls per round (for the queries). Cost is negligible (~$0.0001).

### Why top_k=3 chunks per agent

Each agent gets 3 retrieved passages (~4500 chars). This keeps the prompt under 6K tokens total while providing meaningful literature grounding. With top_k=5, we'd risk exceeding useful context (and most of the extra chunks are low-relevance noise from a 721-chunk corpus).

### Why import puzzle-rag directly (not copy the ChromaDB)

`sys.path.insert(0, .../puzzle-rag)` and `import rag` reuses the existing infrastructure without duplication. The ChromaDB lives at a fixed path. If we moved this to a different machine, we'd need to parameterize the path — but for now, simplicity wins.

## State Design

### Why flat list of responses (not nested dict)

```python
responses: list[dict]  # [{round, agent_id, text, rag_sources}, ...]
```

We considered `responses[round][agent_id] = text` (nested), but a flat list with metadata is:
1. Easier to append to in LangGraph (no reducer complexity)
2. Easier to filter (`[r for r in responses if r["round"] == 2]`)
3. Preserves ordering information naturally
4. JSON-serializable without custom logic

### Why replace full state (not use LangGraph reducers)

LangGraph supports `Annotated[list, add]` reducers that automatically append. We initially planned to use these but found that the state replacement pattern (`return {"responses": state["responses"] + [new_entry]}`) is:
1. More explicit about what's happening
2. Avoids subtle bugs where reducers accumulate duplicates on retries
3. Easier to reason about in an interview explanation

Trade-off: slightly more memory (copying the list each time). Irrelevant at our scale (max 12 entries for a 3-round debate).

## Topology

### Why the specific ring ordering

```
growth_diagram -> yang_baxter -> bijection -> computational -> (loop)
```

Each agent's predecessor represents a *complementary* perspective:
- **Growth Diagram** sees Computational (grounds abstract theory in enumeration)
- **Yang-Baxter** sees Growth Diagram (challenges local rules with spectral structure)
- **Bijection** sees Yang-Baxter (connects integrability to concrete bijection steps)
- **Computational** sees Bijection (tests bijection claims numerically)

This maximizes productive tension. An alternative ordering (GD→Bij→YB→Comp) would pair similar perspectives, reducing the "challenge" effect.

### Why negative coupling as default

From the orchestration project's sign-flip result: negative coupling ("challenge and extend") + ring topology produces maximum diversity. Since this is a research problem where we want diverse approaches, negative coupling is the right default.

Positive coupling ("build on") is available via `--coupling positive` for convergence/synthesis phases.

## Evaluation

### Why DiversityMeter with all-MiniLM-L6-v2

Directly reuses the proven pattern from `~/git/orchestration/experiment.py`. The model is:
1. Fast (384-dim, tiny model)
2. Already validated against the orchestration results
3. Good enough for relative comparisons between topologies

We don't need a math-specific embedding model here — we're measuring *diversity of expression*, not *mathematical similarity*.

### Why cross-domain keyword detection (heuristic)

The `detect_cross_domain_connections()` function is a crude but informative proxy. When agent X uses agent Y's domain vocabulary, it suggests conceptual bridging. This is:
1. Cheap (no API call)
2. Interpretable (you can see exactly which terms triggered)
3. A reasonable interview talking point ("here's how we detect emergent cross-domain insight")

A better approach would use LLM-as-judge to assess whether genuine novel connections were made. That's a natural extension but costs extra API calls per evaluation.

## File Structure

### Why 5 files (not 1 or 10)

- `state.py` — data definitions, no logic. Read this first to understand the system.
- `agents.py` — the "hot path" (RAG + LLM calls). Where the money is spent.
- `graph.py` — topology is visible at a glance. ~50 lines.
- `evaluate.py` — completely separable from execution. Could swap in different metrics.
- `main.py` — CLI + orchestration + output formatting.

This structure maps to interview explanation: "Here's the state, here's how agents work, here's the graph, here's how I evaluate, here's how you run it."

## Experiment 2: Tool Use (added 2026-05-04)

### Why give agents tools

Experiment 1 (text-only debate) produced 4 testable hypotheses. We tested them ourselves in `test_hypotheses.py` and found that agents had fabricated computational claims (e.g., specific choice vectors that were never computed). Giving agents `run_python`, `read_file`, and `write_file` tools lets them verify their own claims in real time.

### Tool design

Three tools via the Anthropic SDK's native tool use:

| Tool | Scope | Purpose |
|------|-------|---------|
| `run_python(code)` | Executes in `~/git/puzzles/` | Test conjectures, enumerate puzzles, compute examples |
| `read_file(path)` | `~/git/puzzles/` + `~/git/category/topology-debate/` | Inspect source code, read prior results |
| `write_file(path, content)` | `~/git/category/topology-debate/` only | Save findings, analysis scripts |

### Why native Anthropic tool use (not LangChain tools)

LangChain/LangGraph have their own tool abstraction (`@tool` decorator, `ToolNode`). We use the Anthropic SDK's native `tools` parameter instead because:

1. No extra dependency (langchain-anthropic not needed)
2. The tool-use loop is explicit — easier to debug and explain in interviews
3. Full control over the loop: max rounds, logging, truncation

### Why max 5 tool rounds per agent turn

Prevents runaway loops where an agent keeps calling tools without producing a final response. 5 rounds is enough for: read code → run experiment → check result → write findings → summarize. In practice, agents typically use 1-3 rounds.

### Why all agents get tools (not just Computational)

Every agent benefits from verification. The Growth Diagram Agent can check if a local rule actually satisfies the axioms. The Bijection Agent can trace Purbhoo's algorithm on concrete examples. Restricting tools to one agent would artificially limit the others.

### Why `cwd=puzzles/` for Python execution

All the mathematical code (`lr_coefficients.py`, `solve_puzzle.py`, `transfer_operators.py`) lives in `~/git/puzzles/`. Running from that directory means agents can `from lr_coefficients import lr_puzzles` with no path manipulation.

### CLI

```bash
# Text-only debate (Experiment 1)
python main.py --rounds 3 --coupling negative --output debate.json

# Tool-augmented debate (Experiment 2)
python main.py --rounds 3 --coupling negative --tools --output debate_tools.json
```

## Future Directions (not implemented)

1. **Alternative topologies**: Star (hub agent coordinates), fully-connected, none (isolation baseline). Change only `graph.py`.
2. **Haskell specification layer**: Type-safe topology definitions that compile to JSON config. Prevents invalid graphs.
3. **Multi-round coupling switching**: Negative coupling for rounds 1-2 (explore), positive for round 3 (converge).
4. **LLM-as-judge evaluation**: Have Claude assess whether each round produced genuine mathematical insight.
5. **Computational agent with tool use**: Actually run `lr_coefficients.py` and report results back to the debate.
