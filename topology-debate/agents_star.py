"""Agent node functions for star topology.

Star topology information flow:
- Center (vanilla) sees ALL peripheral outputs from previous round
- Peripherals see ONLY the center's most recent output
"""

import os
import sys
from pathlib import Path

from anthropic import Anthropic
from dotenv import load_dotenv

from state import AGENT_BY_ID, COUPLING_PROMPTS, DebateState

# Load API key (same logic as agents.py)
if not os.environ.get("ANTHROPIC_API_KEY"):
    creds_path = Path.home() / ".claude" / ".credentials.json"
    if creds_path.exists():
        import json
        creds = json.loads(creds_path.read_text())
        if "claudeAiOauth" in creds:
            os.environ["ANTHROPIC_API_KEY"] = creds["claudeAiOauth"]["accessToken"]
    if not os.environ.get("ANTHROPIC_API_KEY"):
        load_dotenv(Path.home() / "git" / "mary" / ".env")

# Import puzzle-rag
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "puzzles" / "puzzle-rag"))
import rag as puzzle_rag

DEFAULT_MODEL = "claude-opus-4-20250514"


def get_peripheral_outputs(state: DebateState, peripheral_ids: list[str]) -> list[dict]:
    """Get most recent output from each peripheral agent."""
    latest = {}
    for entry in state["responses"]:
        if entry["agent_id"] in peripheral_ids:
            latest[entry["agent_id"]] = entry
    return [latest[pid] for pid in peripheral_ids if pid in latest]


def get_center_output(state: DebateState, center_id: str) -> str | None:
    """Get the most recent output from the center agent."""
    for entry in reversed(state["responses"]):
        if entry["agent_id"] == center_id:
            return entry["text"]
    return None


def make_star_center_node(center_id: str, peripheral_ids: list[str]):
    """Create the center (vanilla) agent node.

    The center sees ALL peripheral outputs from the previous round.
    """
    center_config = AGENT_BY_ID[center_id]

    def center_node(state: DebateState) -> dict:
        client = Anthropic()
        round_num = state["current_round"]
        model = os.environ.get("DEBATE_MODEL", DEFAULT_MODEL)

        # RAG retrieval
        rag_query = center_config.get("rag_seed_query", "")
        rag_results = puzzle_rag.query(rag_query, top_k=3)
        docs = rag_results["documents"][0] if rag_results["documents"][0] else []
        rag_context = "\n---\n".join(docs)

        # Build prompt
        user_parts = [f"## Research Problem\n{state['problem_statement']}"]

        if rag_context:
            user_parts.append(f"\n## Relevant Literature (retrieved for you)\n{rag_context}")

        user_parts.append(f"\n## Round {round_num + 1} of {state['max_rounds']}")

        # Include all peripheral outputs from previous round
        peripheral_outputs = get_peripheral_outputs(state, peripheral_ids)
        if peripheral_outputs:
            user_parts.append("\n## Specialist Analyses (all from previous round)")
            for entry in peripheral_outputs:
                user_parts.append(f"\n### {entry['agent_name']}\n{entry['text']}")
            user_parts.append(
                f"\n## Your Directive\n"
                "You are the central synthesizer. You see ALL specialist perspectives. "
                "Identify the strongest thread, the deepest disagreement, and the most "
                "promising direction. Challenge weak claims. Propose a unified framing "
                "that the specialists should respond to."
            )
        else:
            user_parts.append(
                "\n## Your Directive\n"
                "Provide your initial framing of this problem. What are the key "
                "questions? What would different specialists need to address?"
            )

        user_parts.append(
            "\nRespond in 200-400 words. Be precise and technical. "
            "Cite specific theorems, equations, or computational evidence."
        )

        print(f"  [{center_config['name']}] (center, round {round_num + 1})")
        response = client.messages.create(
            model=model,
            max_tokens=800,
            system=center_config["system_prompt"],
            messages=[{"role": "user", "content": "\n".join(user_parts)}],
        )
        text = response.content[0].text

        return {"responses": [{
            "round": round_num,
            "agent_id": center_id,
            "agent_name": center_config["name"],
            "text": text,
            "rag_query": rag_query,
            "rag_sources": [m.get("source", "") for m in (rag_results["metadatas"][0] or [])],
            "tool_calls": [],
        }]}

    center_node.__name__ = f"agent_{center_id}"
    return center_node


def make_star_peripheral_node(peripheral_id: str, center_id: str):
    """Create a peripheral specialist node.

    Each peripheral sees ONLY the center's most recent output.
    """
    agent_config = AGENT_BY_ID[peripheral_id]

    def peripheral_node(state: DebateState) -> dict:
        client = Anthropic()
        round_num = state["current_round"]
        model = os.environ.get("DEBATE_MODEL", DEFAULT_MODEL)

        # RAG retrieval (adaptive based on center's output)
        center_text = get_center_output(state, center_id)
        base_query = agent_config["rag_seed_query"]
        if center_text and round_num > 0:
            snippet = center_text[-300:].replace("\n", " ")[:100]
            rag_query = f"{base_query} {snippet}"
        else:
            rag_query = base_query

        rag_results = puzzle_rag.query(rag_query, top_k=3)
        docs = rag_results["documents"][0] if rag_results["documents"][0] else []
        rag_context = "\n---\n".join(docs)

        # Build prompt
        user_parts = [f"## Research Problem\n{state['problem_statement']}"]

        if rag_context:
            user_parts.append(f"\n## Relevant Literature (retrieved for you)\n{rag_context}")

        user_parts.append(f"\n## Round {round_num + 1} of {state['max_rounds']}")

        if center_text:
            center_name = AGENT_BY_ID[center_id]["name"]
            user_parts.append(f"\n## {center_name}'s Synthesis\n{center_text}")
            user_parts.append(f"\n## Your Directive\n{COUPLING_PROMPTS[state['coupling']]}")
        else:
            user_parts.append(
                "\n## Your Directive\n"
                "Provide your initial analysis from your specialist perspective. "
                "What do you see in this problem that others might miss?"
            )

        user_parts.append(
            "\nRespond in 200-400 words. Be precise and technical. "
            "Cite specific theorems, equations, or computational evidence."
        )

        print(f"  [{agent_config['name']}] (peripheral, round {round_num + 1})")
        response = client.messages.create(
            model=model,
            max_tokens=800,
            system=agent_config["system_prompt"],
            messages=[{"role": "user", "content": "\n".join(user_parts)}],
        )
        text = response.content[0].text

        return {"responses": [{
            "round": round_num,
            "agent_id": peripheral_id,
            "agent_name": agent_config["name"],
            "text": text,
            "rag_query": rag_query,
            "rag_sources": [m.get("source", "") for m in (rag_results["metadatas"][0] or [])],
            "tool_calls": [],
        }]}

    peripheral_node.__name__ = f"agent_{peripheral_id}"
    return peripheral_node
