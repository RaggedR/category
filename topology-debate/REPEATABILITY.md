# Repeatability Guide

Instructions for reproducing the topology-debate case study.

## Prerequisites

- Python 3.10+
- Anthropic API key with credits (set `ANTHROPIC_API_KEY`)
- OpenAI API key for RAG embeddings (set `OPENAI_API_KEY`)

## Dependencies

```bash
cd ~/git/category/topology-debate
pip install -r requirements.txt
```

Key packages: `anthropic`, `langgraph`, `sentence-transformers`, `chromadb`.

## RAG Vector Database

The experiment uses a ChromaDB vector store at `~/git/puzzles/puzzle-rag/chroma_db/`.

### Papers in the corpus (9 papers, 721 chunks)

**Core puzzles:**
- Knutson-Tao (2003) — Puzzles and equivariant cohomology
- Zinn-Justin (2009) — LR coefficients via integrable tilings

**Extensions:**
- Knutson-Zinn-Justin (2017) — Schubert puzzles
- Wheeler-Zinn-Justin (2016) — Hall-inverse Kostka puzzles
- Wheeler-Zinn-Justin (2019) — Grothendieck LR

**Vertex models:**
- Motegi-Sakai (2015) — Partition functions and symmetric functions
- Wheeler-Zinn-Justin (2015) — Refined Cauchy and Hall-Littlewood
- Zinn-Justin (2009) — Six-vertex model survey
- Zinn-Justin (2024) — Integrability and combinatorics

### Rebuilding the vector store

If the ChromaDB needs rebuilding:

```bash
cd ~/git/puzzles/puzzle-rag
python rag.py --rebuild
```

This re-embeds all PDFs using `text-embedding-3-small` (OpenAI). Cost: ~$0.01.

## Running the experiments

### Ring topology (directed cycle, 4 agents, sequential)

```bash
python main.py --rounds 3 --coupling negative --output ring_run_1.json
```

### Star topology (5 agents: 1 hub + 4 specialists, hub bidirectional)

```bash
python run_star.py --rounds 3 --coupling negative --output star_run_1.json
```

### Multiple seeds

Run each command multiple times with different `--output` filenames. Each run uses a fresh API call sequence; non-determinism comes from LLM sampling.

## Evaluating results

```bash
python evaluate.py ring_run_1.json star_run_1.json
```

This computes per-round diversity (mean pairwise cosine distance, all-MiniLM-L6-v2) and cross-domain connection counts.

## Scripts

| File | Purpose |
|------|---------|
| `main.py` | Ring experiment CLI |
| `run_star.py` | Star experiment CLI |
| `evaluate.py` | Diversity + cross-domain metrics |
| `state.py` | Agent definitions, state schema |
| `agents.py` | Ring agent nodes (sequential, RAG-augmented) |
| `agents_star.py` | Star agent nodes (center + parallel peripherals) |
| `graph.py` | Ring LangGraph construction |
| `graph_star.py` | Star LangGraph construction |
| `tools.py` | Optional tool definitions (run_python, read_file, write_file) |

## Existing results

| File | Topology | Rounds | Model |
|------|----------|--------|-------|
| `opus_debate.json` | Ring | 3 | Claude Opus |
| `star_debate.json` | Star | 3 | Claude Opus |

## Cost estimate

- Ring (3 rounds, 4 agents): ~$0.45 per run
- Star (3 rounds, 5 agents): ~$0.55 per run
- 10 seeds of each: ~$10 total
