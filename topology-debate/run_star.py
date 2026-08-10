#!/usr/bin/env python3
"""Run the star-topology debate experiment.

Usage:
    python run_star.py --rounds 3 --coupling negative --output star_debate.json
"""

import argparse
import json
import time
from pathlib import Path

from state import DebateState, RING_ORDER
from graph_star import build_star_graph, CENTER_ID, PERIPHERAL_IDS
from evaluate import DiversityMeter, detect_cross_domain_connections


def run_star_debate(max_rounds: int = 3, coupling: str = "negative",
                    output_path: str | None = None):
    """Run the star-topology debate and evaluate."""
    print(f"Star-Topology Debate: {max_rounds} rounds, {coupling} coupling, text only")
    print(f"Center: {CENTER_ID} | Peripherals: {', '.join(PERIPHERAL_IDS)}")
    print(f"Edges: {CENTER_ID} <-> each peripheral (4 edges)")
    print(f"{'='*60}\n")

    graph = build_star_graph(max_rounds=max_rounds)

    # Load problem statement
    problem_path = Path(__file__).resolve().parent.parent.parent / "puzzles" / "OPEN_QUESTION.md"
    problem_statement = problem_path.read_text()

    # Initial state
    initial_state: DebateState = {
        "problem_statement": problem_statement,
        "max_rounds": max_rounds,
        "coupling": coupling,
        "use_tools": False,
        "current_round": 0,
        "current_agent_idx": 0,
        "responses": [],
        "metrics": [],
    }

    # Execute
    t0 = time.time()
    print("Running debate...\n")
    final_state = graph.invoke(initial_state)
    elapsed = time.time() - t0

    print(f"\nDebate complete in {elapsed:.1f}s")
    print(f"Total responses: {len(final_state['responses'])}")

    # -- Evaluation --
    print(f"\n{'='*60}")
    print("EVALUATION")
    print(f"{'='*60}")

    meter = DiversityMeter()
    metrics = []

    for r in range(max_rounds):
        round_responses = [e for e in final_state["responses"] if e["round"] == r]
        if not round_responses:
            continue

        texts = [e["text"] for e in round_responses]
        labels = [e["agent_id"] for e in round_responses]
        diversity = meter.measure(texts, labels)
        connections = detect_cross_domain_connections(round_responses)

        metrics.append({
            "round": r + 1,
            "diversity": diversity["mean_distance"],
            "connections": connections,
        })

        print(f"\nRound {r + 1}:")
        print(f"  Diversity (mean cosine distance): {diversity['mean_distance']:.4f}")
        print(f"  Cross-domain connections: {len(connections)}")
        for c in connections:
            print(f"    - {c}")

    # -- Save results --
    result = {
        "config": {
            "topology": "star",
            "max_rounds": max_rounds,
            "coupling": coupling,
            "center": CENTER_ID,
            "peripherals": PERIPHERAL_IDS,
        },
        "responses": final_state["responses"],
        "metrics": metrics,
        "elapsed_seconds": elapsed,
    }

    out = output_path or "star_debate.json"
    Path(out).write_text(json.dumps(result, indent=2))
    print(f"\nResults saved to {out}")

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Star-topology debate experiment")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--coupling", choices=["negative", "positive"], default="negative")
    parser.add_argument("--output", type=str, default="star_debate.json")
    args = parser.parse_args()

    run_star_debate(
        max_rounds=args.rounds,
        coupling=args.coupling,
        output_path=args.output,
    )
