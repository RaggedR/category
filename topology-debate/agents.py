"""Agent node functions with RAG-augmented generation and tool use."""

import os
import sys
from pathlib import Path

from anthropic import Anthropic
from dotenv import load_dotenv

from state import AGENTS, AGENT_BY_ID, COUPLING_PROMPTS, DebateState
from tools import TOOL_DEFINITIONS, execute_tool

# Load API key: check env, then ~/.claude/.credentials.json OAuth, then ~/git/mary/.env
if not os.environ.get("ANTHROPIC_API_KEY"):
    creds_path = Path.home() / ".claude" / ".credentials.json"
    if creds_path.exists():
        import json
        creds = json.loads(creds_path.read_text())
        if "claudeAiOauth" in creds:
            os.environ["ANTHROPIC_API_KEY"] = creds["claudeAiOauth"]["accessToken"]
    if not os.environ.get("ANTHROPIC_API_KEY"):
        load_dotenv(Path.home() / "git" / "mary" / ".env")

# Import the puzzle-rag query function
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "puzzles" / "puzzle-rag"))
import rag as puzzle_rag

# Model selection — override via DEBATE_MODEL env var
DEFAULT_MODEL = "claude-opus-4-20250514"

# Max tool-use rounds per agent turn (prevents runaway loops)
MAX_TOOL_ROUNDS = 5


def get_predecessor_output(state: DebateState, agent_config: dict) -> str | None:
    """Get the most recent output from this agent's predecessor in the ring."""
    pred_id = agent_config["predecessor"]
    for entry in reversed(state["responses"]):
        if entry["agent_id"] == pred_id:
            return entry["text"]
    return None


def build_rag_query(agent_config: dict, predecessor_text: str | None, round_num: int) -> str:
    """Construct a RAG query that adapts based on predecessor's output."""
    base_query = agent_config["rag_seed_query"]
    if predecessor_text and round_num > 0:
        snippet = predecessor_text[-300:].replace("\n", " ")[:100]
        return f"{base_query} {snippet}"
    return base_query


def run_tool_use_loop(client, model: str, system: str, messages: list,
                      agent_name: str) -> tuple[str, list[dict]]:
    """Run a tool-use loop until the agent produces a final text response.

    Returns (final_text, tool_calls_made).
    """
    tool_calls = []

    for round_num in range(MAX_TOOL_ROUNDS):
        response = client.messages.create(
            model=model,
            max_tokens=2000,
            system=system,
            messages=messages,
            tools=TOOL_DEFINITIONS,
        )

        # Collect text and tool_use blocks from the response
        text_parts = []
        tool_uses = []
        for block in response.content:
            if block.type == "text":
                text_parts.append(block.text)
            elif block.type == "tool_use":
                tool_uses.append(block)

        # If no tool calls, we're done — return the text
        if not tool_uses:
            return "\n".join(text_parts), tool_calls

        # Process each tool call
        # First, add the assistant's response (with tool_use blocks) to messages
        messages.append({"role": "assistant", "content": response.content})

        # Then add tool results
        tool_results = []
        for tool_use in tool_uses:
            print(f"    [{agent_name}] calling {tool_use.name}()")
            result = execute_tool(tool_use.name, tool_use.input)
            tool_calls.append({
                "tool": tool_use.name,
                "input": tool_use.input,
                "output": result[:500],
            })
            tool_results.append({
                "type": "tool_result",
                "tool_use_id": tool_use.id,
                "content": result,
            })

        messages.append({"role": "user", "content": tool_results})

    # Hit max rounds — force a final response WITHOUT tools
    print(f"    [{agent_name}] max tool rounds reached, forcing summary")
    messages.append({
        "role": "user",
        "content": (
            "You've used all your tool rounds. Now write your final analysis "
            "based on everything you computed and read. Be concise (200-400 words). "
            "Include specific numbers and results from your computations."
        ),
    })
    response = client.messages.create(
        model=model,
        max_tokens=1000,
        system=system,
        messages=messages,
        # No tools — forces text output
    )
    final_text = "\n".join(b.text for b in response.content if b.type == "text")
    return final_text, tool_calls


def make_agent_node(agent_idx: int):
    """Factory: create a LangGraph node function for the agent at position agent_idx."""
    agent_config = AGENTS[agent_idx]
    agent_id = agent_config["id"]

    def agent_node(state: DebateState) -> dict:
        client = Anthropic()
        round_num = state["current_round"]
        use_tools = state.get("use_tools", False)

        # 1. Get predecessor's output (ring topology constraint)
        pred_output = get_predecessor_output(state, agent_config)

        # 2. RAG retrieval (adaptive query)
        rag_query = build_rag_query(agent_config, pred_output, round_num)
        rag_results = puzzle_rag.query(rag_query, top_k=3)

        docs = rag_results["documents"][0] if rag_results["documents"][0] else []
        metas = rag_results["metadatas"][0] if rag_results["metadatas"][0] else []
        rag_context = "\n---\n".join(docs)

        # 3. Build the prompt
        user_parts = [
            f"## Research Problem\n{state['problem_statement']}",
        ]

        if rag_context:
            user_parts.append(f"\n## Relevant Literature (retrieved for you)\n{rag_context}")

        user_parts.append(f"\n## Round {round_num + 1} of {state['max_rounds']}")

        if pred_output:
            pred_name = AGENT_BY_ID[agent_config["predecessor"]]["name"]
            user_parts.append(f"\n## {pred_name}'s Analysis\n{pred_output}")
            user_parts.append(f"\n## Your Directive\n{COUPLING_PROMPTS[state['coupling']]}")
        else:
            user_parts.append(
                "\n## Your Directive\n"
                "Provide your initial analysis from your specialist perspective. "
                "What do you see in this problem that others might miss?"
            )

        if use_tools:
            user_parts.append(
                "\nYou have access to tools: run_python (execute code), "
                "read_file (inspect source), write_file (save results). "
                "Use run_python to TEST your claims with real computation — "
                "don't just theorize, verify! The working directory is ~/git/puzzles/ "
                "so you can `from lr_coefficients import lr_puzzles` directly. "
                "After computing, write a concise summary of your findings."
            )
        else:
            user_parts.append(
                "\nRespond in 200-400 words. Be precise and technical. "
                "Cite specific theorems, equations, or computational evidence."
            )

        # 4. Call Claude (with or without tools)
        model = os.environ.get("DEBATE_MODEL", DEFAULT_MODEL)
        messages = [{"role": "user", "content": "\n".join(user_parts)}]

        if use_tools:
            print(f"  [{agent_config['name']}] (with tools)")
            text, tool_calls = run_tool_use_loop(
                client, model,
                agent_config["system_prompt"],
                messages,
                agent_config["name"],
            )
        else:
            response = client.messages.create(
                model=model,
                max_tokens=800,
                system=agent_config["system_prompt"],
                messages=messages,
            )
            text = response.content[0].text
            tool_calls = []

        # 5. Return state update (responses is Annotated[list, add])
        return {
            "responses": [{
                "round": round_num,
                "agent_id": agent_id,
                "agent_name": agent_config["name"],
                "text": text,
                "rag_query": rag_query,
                "rag_sources": [m.get("source", "") for m in metas],
                "tool_calls": tool_calls,
            }],
            "current_agent_idx": (agent_idx + 1) % len(AGENTS),
        }

    agent_node.__name__ = f"agent_{agent_id}"
    return agent_node
