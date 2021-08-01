import itertools
from typing import List, Dict, Callable, Tuple, Set

from matplotlib import pyplot as plt

from src.model.fragment import Fragment
from src.model.precursor import Precursor
from src.utilities.constants import get_golden_bonds
from src.utilities.dataloading import load_protein, load_fragment_matches

import networkx as nx


def intersects(t1, t2):
    return t1[0] in t2 or t1[1] in t2


def normalize(x, xs):
    bottom = abs(max(xs) - min(xs))
    return abs(x - min(xs)) / (bottom if bottom else 1)


def draw(graph, node_scores, edge_scores, ax, cmap=plt.get_cmap("PiYG")):
    node_colors = [cmap(float(s)) for s in node_scores]
    edge_colors = [cmap(float(s)) for s in edge_scores]
    return nx.draw_circular(
        graph,
        ax=ax,
        with_labels=True,
        node_size=1200,
        font_size=12,
        node_color="#E7DBB7",
        linewidths=2,
        arrowsize=10,
        edgecolors=node_colors,
        edge_color=edge_colors,
        width=[0.3 if sc < 0.5 else 2 for sc in edge_scores],
        connectionstyle="arc3,rad=0.2",
    )


def build_graph(cysteines: List[int]):
    nodes = list(reversed(cysteines[3:] + cysteines[:3]))
    return nx.DiGraph(nx.complete_graph(nodes))


def nodes_edges_from_golden(graph, golden_bonds: List):
    node_scores = [not any(n in b for b in golden_bonds) for n in graph.nodes()]
    edge_scores = [tuple(sorted(e)) in golden_bonds for e in graph.edges()]
    return node_scores, edge_scores


def nodes_edges_from_data(graph, positive_evidence: Dict, alkylation_evidence: Dict):
    node_scores = [(n, alkylation_evidence[n]) for n in graph.nodes()]
    edge_scores = [(e, positive_evidence[tuple(sorted(e))]) for e in graph.edges()]

    edge_scores_norm = [
        normalize(
            esc,
            [n for f, n in edge_scores if f[0] == e[0]]
            + [nsc for node, nsc in node_scores if node == e[0]],
        )
        for e, esc in edge_scores
    ]
    node_scores_norm = [
        normalize(s, [n for f, n in edge_scores if f[0] == node])
        for node, s in node_scores
    ]
    return node_scores_norm, edge_scores_norm


def calculate_scores(
    scored_fragment_matches: List[Dict],
    cysteines: List[int],
    calc_weight: Callable,
    positive_ev_w: float = 1,
    indirect_positive_ev_w: float = 1,
    alkylation_w: float = 1,
):
    possible_bonds = list(itertools.combinations(cysteines, 2))
    positive_evidence = {b: 0 for b in possible_bonds}
    alkylation_evidence = {c: 0 for c in cysteines}

    for match in scored_fragment_matches:
        fragment: Fragment = match["fragment"]

        if fragment is None:
            continue

        match_weight = calc_weight(match)
        cys_in_bonds = []

        for b in fragment.connected_bonds:
            positive_evidence[b] += match_weight * positive_ev_w
            cys_in_bonds += list(b)

        # TODO: Check if only interesting matter
        for c in fragment.disconnected_cys:
            cys_in_bonds.append(c)
            for b in positive_evidence:
                if c in b:
                    positive_evidence[b] += match_weight * indirect_positive_ev_w

        for c in cysteines:
            if c not in cys_in_bonds and any(
                b <= c < e for b, e in fragment.residue_ranges
            ):
                alkylation_evidence[c] += match_weight * alkylation_w

    return positive_evidence, alkylation_evidence


def numerator(match: Dict):
    lowest = {
        (63, 114),
        (5, 75),
        (5, 79),
        (75, 114),
        (75, 126),
        (79, 114),
        (79, 126),
        (93, 114),
        (93, 126),
    }
    middle = {(5, 126), (5, 63), (5, 93)}

    frag: Fragment = match["fragment"]
    bonds = set(frag.connected_bonds)

    if not bonds.isdisjoint(lowest):
        return 0.1
    elif not bonds.isdisjoint(middle):
        return 0.4
    else:
        return 1


def score_match(match: Dict, numerator=1):
    prec: Precursor = match["precursor"]
    return numerator / (prec.error_ppm + prec.to_dict()["prec_max_mc_count"] + 1)


if __name__ == "__main__":
    import argparse

    args = argparse.ArgumentParser(
        description="Visualize disulhpide bonds for the given protein"
    )
    args.add_argument(
        "--protein",
        type=str,
        required=True,
        help="protein code (usually three letters)",
    )
    args.add_argument(
        "--kind",
        type=str,
        choices=["AT", "RAT"],
        required=True,
        help="measurement type (AT/RAT)",
    )
    args.add_argument(
        "--prec_error",
        type=int,
        required=True,
        help="allowed measurement error in ppm",
    )
    args.add_argument(
        "--frag_error",
        type=int,
        required=True,
        help="allowed measurement error in ppm",
    )
    args.add_argument(
        "--prec_segments",
        type=int,
        required=True,
        help="upper bound of segment count in matched precursors",
    )
    args.add_argument(
        "--frag_breaks",
        type=int,
        required=True,
        help="upper bound of break count in matched fragments",
    )
    args.add_argument(
        "--code",
        type=str,
        default=None,
        required=False,
        help="code to append to the output file name",
    )
    args = args.parse_args()

    prot_sequence = load_protein(args.protein)
    cysteines = [i for i, res in enumerate(prot_sequence) if res == "C"]
    golden_bonds = get_golden_bonds(args.protein)

    fragment_matches = load_fragment_matches(
        protein=args.protein,
        kind=args.kind,
        segments=args.prec_segments,
        breaks=args.frag_breaks,
        error=args.frag_error,
        code=args.code,
    )

    if args.kind == "AT":
        scoring_fun = lambda m: score_match(m, numerator(m))
    else:
        scoring_fun = score_match

    positive_evidence, alkylation_evidence = calculate_scores(
        scored_fragment_matches=fragment_matches,
        cysteines=cysteines,
        calc_weight=scoring_fun,
        indirect_positive_ev_w=0,
        alkylation_w=0.5,
    )

    graph = build_graph(cysteines=cysteines)
    gn, ge = nodes_edges_from_golden(graph, golden_bonds)
    n, e = nodes_edges_from_data(graph, positive_evidence, alkylation_evidence)

    fig, axs = plt.subplots(1, ncols=2, figsize=(10, 5), dpi=300)

    for ax in axs:
        ax.set_aspect("equal")

    draw(graph, gn, ge, axs[0])
    draw(graph, n, e, axs[1])

    image_path = "../out/plots/{}_{}_segments={}_breaks={}_perr={}_ferr={}{}".format(
        args.protein,
        args.kind,
        args.prec_segments,
        args.frag_breaks,
        args.prec_error,
        args.frag_error,
        "" if args.code is None else f"_{args.code}",
    )

    fig.suptitle(f"{args.protein} {args.kind}", fontsize=24)

    fig.savefig(image_path)
