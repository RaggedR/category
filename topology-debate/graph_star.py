"""LangGraph StateGraph construction with star topology.

Star topology (5 nodes, 4 edges):
- Center: vanilla agent (synthesizer, no domain bias)
- Peripherals: 4 specialist agents (same as ring experiment)
- Edges: center <-> each peripheral (4 edges, matching ring's 4 edges)

Information flow per round:
  1. Center speaks (sees all 4 peripherals from previous round)
  2. All 4 peripherals respond in parallel (each sees only center)
"""

from langgraph.graph import StateGraph, START, END

from state import DebateState, RING_ORDER
from agents_star import make_star_center_node, make_star_peripheral_node

# Star config: vanilla is the hub, 4 specialists are peripherals
CENTER_ID = "vanilla"
PERIPHERAL_IDS = RING_ORDER  # all 4 specialists


def should_continue(state: DebateState) -> str:
    """After peripherals complete, loop or finish."""
    next_round = state["current_round"] + 1
    if next_round >= state["max_rounds"]:
        return "done"
    return "next_round"


def round_transition(state: DebateState) -> dict:
    """Increment round counter."""
    return {"current_round": state["current_round"] + 1}


def build_star_graph(max_rounds: int = 3):
    """Construct the star-topology debate graph.

    Structure per round:
      center -> [peripheral_1, peripheral_2, peripheral_3, peripheral_4] (fan-out)
      all peripherals -> join -> loop/end (fan-in)
    """
    builder = StateGraph(DebateState)

    # -- Nodes --
    center_name = f"agent_{CENTER_ID}"
    builder.add_node(center_name, make_star_center_node(CENTER_ID, PERIPHERAL_IDS))

    peripheral_names = []
    for pid in PERIPHERAL_IDS:
        name = f"agent_{pid}"
        builder.add_node(name, make_star_peripheral_node(pid, CENTER_ID))
        peripheral_names.append(name)

    def join_node(state: DebateState) -> dict:
        """No-op join after all peripherals complete."""
        return {}

    builder.add_node("join", join_node)
    builder.add_node("round_transition", round_transition)

    # -- Edges --
    builder.add_edge(START, center_name)

    # Center -> all peripherals (fan-out)
    for pname in peripheral_names:
        builder.add_edge(center_name, pname)

    # All peripherals -> join (fan-in)
    for pname in peripheral_names:
        builder.add_edge(pname, "join")

    # Join -> conditional: loop or end
    builder.add_conditional_edges(
        "join",
        should_continue,
        {
            "next_round": "round_transition",
            "done": END,
        },
    )

    # Round transition -> center (next round)
    builder.add_edge("round_transition", center_name)

    return builder.compile()
