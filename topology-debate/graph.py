"""LangGraph StateGraph construction with ring topology."""

from langgraph.graph import StateGraph, START, END

from state import DebateState, RING_ORDER
from agents import make_agent_node


def should_continue(state: DebateState) -> str:
    """Route after the last agent: loop for another round or finish."""
    # current_agent_idx wraps to 0 after the last agent completes
    if state["current_agent_idx"] == 0:
        next_round = state["current_round"] + 1
        if next_round >= state["max_rounds"]:
            return "done"
        return "next_round"
    return "done"


def round_transition(state: DebateState) -> dict:
    """Increment the round counter between rounds."""
    return {"current_round": state["current_round"] + 1}


def build_graph(max_rounds: int = 3):
    """Construct the ring-topology debate graph.

    Ring order: growth_diagram -> yang_baxter -> bijection -> computational
    Each agent sees only its predecessor's output.

    After the last agent, either loop (increment round) or end.
    """
    builder = StateGraph(DebateState)

    # -- Nodes --
    node_names = []
    for i, agent_id in enumerate(RING_ORDER):
        name = f"agent_{agent_id}"
        builder.add_node(name, make_agent_node(i))
        node_names.append(name)

    builder.add_node("round_transition", round_transition)

    # -- Edges: linear chain within a round --
    builder.add_edge(START, node_names[0])
    for i in range(len(node_names) - 1):
        builder.add_edge(node_names[i], node_names[i + 1])

    # -- Conditional: after last agent, loop or finish --
    builder.add_conditional_edges(
        node_names[-1],
        should_continue,
        {
            "next_round": "round_transition",
            "done": END,
        },
    )

    # Round transition feeds back to first agent
    builder.add_edge("round_transition", node_names[0])

    return builder.compile()
