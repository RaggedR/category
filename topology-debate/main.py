#!/usr/bin/env python3
"""CLI entry point for the ring-topology multi-agent debate.

Usage:
    python main.py --rounds 3 --coupling negative --output debate.json
    python main.py --rounds 1 --coupling positive   # quick test
"""

import argparse
import json
import time
from pathlib import Path

from state import DebateState, RING_ORDER
from graph import build_graph
from evaluate import DiversityMeter, detect_cross_domain_connections


def run_debate(max_rounds: int = 3, coupling: str = "negative",
               output_path: str | None = None, use_tools: bool = False):
    """Run the full ring-topology debate and evaluate."""
    mode = "with tools" if use_tools else "text only"
    print(f"Ring-Topology Debate: {max_rounds} rounds, {coupling} coupling, {mode}")
    print(f"Ring: {' -> '.join(RING_ORDER)} -> (loop)")
    print(f"{'='*60}\n")

    graph = build_graph(max_rounds=max_rounds)

    # Load the problem statement
    problem_path = Path(__file__).resolve().parent.parent.parent / "puzzles" / "OPEN_QUESTION.md"
    problem_statement = problem_path.read_text()

    # Initial state
    initial_state: DebateState = {
        "problem_statement": problem_statement,
        "max_rounds": max_rounds,
        "coupling": coupling,
        "use_tools": use_tools,
        "current_round": 0,
        "current_agent_idx": 0,
        "responses": [],
        "metrics": [],
    }

    # Execute the graph
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
    all_metrics = []

    for round_num in range(max_rounds):
        round_responses = [r for r in final_state["responses"] if r["round"] == round_num]
        if len(round_responses) < 2:
            continue

        texts = [r["text"] for r in round_responses]
        labels = [r["agent_id"] for r in round_responses]
        diversity = meter.measure(texts, labels)
        connections = detect_cross_domain_connections(round_responses)

        round_metrics = {
            "round": round_num + 1,
            "diversity": diversity["mean_distance"],
            "connections": connections,
        }
        all_metrics.append(round_metrics)

        print(f"\n--- Round {round_num + 1} ---")
        print(f"  Diversity (mean pairwise distance): {diversity['mean_distance']:.4f}")
        print(f"  Range: [{diversity['min_distance']:.4f}, {diversity['max_distance']:.4f}]")

        if connections:
            print(f"  Cross-domain connections:")
            for c in connections:
                print(f"    - {c}")

        # Print agent summaries
        for r in round_responses:
            first_line = r["text"].split("\n")[0][:80]
            print(f"  [{r['agent_id']}] {first_line}...")
            print(f"    RAG sources: {r.get('rag_sources', [])}")
            if r.get("tool_calls"):
                print(f"    Tool calls: {len(r['tool_calls'])}")
                for tc in r["tool_calls"]:
                    print(f"      - {tc['tool']}()")

    # Diversity trend
    if len(all_metrics) > 1:
        print(f"\n--- Diversity Trend ---")
        for m in all_metrics:
            bar = "#" * int(m["diversity"] * 50)
            print(f"  Round {m['round']}: {m['diversity']:.4f} {bar}")

    # Save full output
    if output_path:
        output = {
            "config": {
                "topology": "ring",
                "max_rounds": max_rounds,
                "coupling": coupling,
                "ring_order": RING_ORDER,
            },
            "responses": final_state["responses"],
            "metrics": all_metrics,
            "elapsed_seconds": elapsed,
        }
        Path(output_path).write_text(json.dumps(output, indent=2))
        print(f"\nFull output saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Ring-topology multi-agent debate on KTW puzzle conjecture",
    )
    parser.add_argument(
        "--rounds", type=int, default=3,
        help="Number of debate rounds (default: 3)",
    )
    parser.add_argument(
        "--coupling", choices=["negative", "positive"], default="negative",
        help="Coupling sign: negative=challenge, positive=build (default: negative)",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Save full JSON output to file",
    )
    parser.add_argument(
        "--tools", action="store_true",
        help="Give agents Python execution, file read/write tools",
    )
    args = parser.parse_args()
    run_debate(max_rounds=args.rounds, coupling=args.coupling,
               output_path=args.output, use_tools=args.tools)


if __name__ == "__main__":
    main()
