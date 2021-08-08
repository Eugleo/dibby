import itertools
from typing import List, Dict, Callable
import matplotlib.ticker as mticker

import numpy as np
import numpy.lib.recfunctions as rf
import pandas as pd
from matplotlib import pyplot as plt, cm

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


def draw(graph, node_scores, edge_scores, ax, as_matrix: bool):
    if as_matrix:
        cmap = plt.get_cmap("GnBu")
    else:
        cmap = plt.get_cmap("PiYG")

    if as_matrix:
        ns = len(graph.nodes())
        matrix = np.zeros([ns, ns, 4])
        encode = {n: i for i, n in enumerate(graph.nodes())}
        encoded_edge_scores = {
            tuple(encode[x] for x in edge): score
            for score, edge in zip(edge_scores, graph.edges())
        }
        for index, _ in np.ndenumerate(matrix):
            x, y, i = index
            if x == y:
                matrix[index] = cmap(float(node_scores[x]))[i]
            else:
                matrix[index] = cmap(float(encoded_edge_scores[(x, y)]))[i]
        im = ax.imshow(matrix)
        ax.set_xticks(np.arange(ns))
        ax.set_xticks(np.arange(ns) + 0.5, minor=True)
        ax.set_yticks(np.arange(ns))
        ax.set_yticks(np.arange(ns) + 0.5, minor=True)
        ax.set_xticklabels(graph.nodes())
        ax.set_yticklabels(graph.nodes())
        ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
        ax.spines[:].set_visible(False)
        # ax.minorticks_on()
        ax.grid(which="minor", color="white", linestyle="-", linewidth=3)
        ax.tick_params(which="minor", bottom=False, left=False)
        return im
    else:
        node_colors = [cmap(float(s)) for s in node_scores]
        edge_colors = [cmap(float(s)) for s in edge_scores]
        return nx.draw_circular(
            graph,
            ax=ax,
            with_labels=True,
            node_size=1200,
            font_size=12,
            node_color="#E7DBB7",
            linewidths=4,
            arrowsize=10,
            edgecolors=node_colors,
            edge_color=edge_colors,
            width=[0.5 if sc < 0.5 else 4 for sc in edge_scores],
            connectionstyle="arc3,rad=0.2",
        )


def build_graph(cysteines: List[int], as_matrix: bool):
    if as_matrix:
        nodes = sorted(cysteines)
    else:
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
        if "frag_sequence" not in match:
            continue

        match_weight = calc_weight(match)
        cys_in_bonds = []

        for b in match["frag_connected_bonds"]:
            positive_evidence[b] += match_weight * positive_ev_w
            cys_in_bonds += list(b)

        # TODO: Check if only interesting matter
        for c in match["frag_interesting_disconnected_cys"]:
            cys_in_bonds.append(c)
            for b in positive_evidence:
                if c in b:
                    positive_evidence[b] += match_weight * indirect_positive_ev_w

        for c in cysteines:
            if c not in cys_in_bonds and any(
                b <= c < e for b, e in match["frag_residue_ranges"]
            ):
                alkylation_evidence[c] += match_weight * alkylation_w

    return positive_evidence, alkylation_evidence


def numerator(match: Dict):
    lowest = {(75, 114), (79, 114), (79, 126), (93, 114)} | {
        (29, 119),
        (29, 381),
        (72, 381),
    }
    middle = {(5, 126), (93, 126), (75, 126)}

    bonds = set(match["frag_connected_bonds"])

    if not bonds.isdisjoint(lowest):
        return 0.1
    elif not bonds.isdisjoint(middle):
        return 0.8
    else:
        return 1


def score_match(match: Dict, numerator=1):
    return (
        numerator / match["prec_antiscore"]
        + 0.5 * numerator / match["prec_median_frag_antiscore"]
    )


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
    args.add_argument(
        "--matrix",
        type=bool,
        default=False,
        required=False,
        help="render as a matrix",
        action=argparse.BooleanOptionalAction,
    )
    args = args.parse_args()

    prot_sequence = load_protein(args.protein)
    cysteines = [i for i, res in enumerate(prot_sequence) if res == "C"]
    golden_bonds = get_golden_bonds(args.protein)

    fig, axs = plt.subplots(2, 2, figsize=(10, 10), dpi=300, constrained_layout=True)
    for ax in axs.flat:
        ax.set_aspect("equal")

    axs[0, 0].set_title("RAT, gold", size=16)
    axs[0, 1].set_title("RAT, computed", size=16)
    axs[1, 0].set_title("AT, gold", size=16)
    axs[1, 1].set_title("AT, computed", size=16)

    for kind, axes in [("RAT", axs[0, :]), ("AT", axs[1, :])]:
        fragment_matches = load_fragment_matches(
            protein=args.protein,
            kind=kind,
            segments=args.prec_segments,
            breaks=args.frag_breaks,
            error=args.frag_error,
            code=args.code,
        )

        if kind == "AT":
            scoring_fun = lambda m: score_match(m, numerator(m))
        else:
            scoring_fun = score_match

        fragment_records = [
            fm["scan"].to_dict()
            | fm["precursor"].to_dict()
            | fm["variant"].to_dict()
            | (fm["fragment"].to_dict())
            | {"prec_variant_count": fm["variant_count"]}
            for fm in fragment_matches
            if fm["fragment"] is not None
        ]
        df = pd.DataFrame(fragment_records)

        df["frag_mod_count"] = [len(r) for r in df["frag_mods"]]
        cols = [
            "frag_charge",
            "frag_error_ppm",
            "frag_mod_count",
            "frag_break_count",
            "prec_mass",
            "prec_max_mc_count",
            "prec_error",
            "prec_variant_count",
        ]
        for col in cols:
            c = df[col]
            df[col + "_norm"] = (c - np.min(c)) / (np.max(c) - np.min(c))

        df["frag_antiscore"] = (
            1
            + 16 * df["frag_charge_norm"]
            + 4 * df["frag_error_ppm_norm"]
            + 4 * df["frag_mod_count_norm"]
        )

        df["prec_antiscore"] = (
            1
            + 32 * df["prec_variant_count_norm"]
            + 4 * df["prec_max_mc_count_norm"]
            + 4 * df["prec_mass_norm"]
            + 4 * df["prec_error_norm"]
        )

        df["var_bonds"] = [str(r) for r in df["var_bonds"]]
        df["prec_median_frag_antiscore"] = df.groupby(
            ["scan_id", "prec_sequence", "var_bonds"], as_index=False
        )["frag_antiscore"].transform(np.median)

        positive_evidence, alkylation_evidence = calculate_scores(
            scored_fragment_matches=df.to_dict(orient="records"),
            cysteines=cysteines,
            calc_weight=scoring_fun,
            indirect_positive_ev_w=0,
            alkylation_w=1,
        )

        graph = build_graph(cysteines=cysteines, as_matrix=args.matrix)
        gn, ge = nodes_edges_from_golden(graph, golden_bonds if kind == "AT" else [])
        n, e = nodes_edges_from_data(graph, positive_evidence, alkylation_evidence)

        im = draw(graph, gn, ge, axes[0], as_matrix=args.matrix)
        draw(graph, n, e, axes[1], as_matrix=args.matrix)
        if args.matrix:
            cmap = plt.get_cmap("GnBu")
        else:
            cmap = plt.get_cmap("PiYG")
    image_path = "../out/plots/{}_segments={}_breaks={}_perr={}_ferr={}{}.pdf".format(
        args.protein,
        args.prec_segments,
        args.frag_breaks,
        args.prec_error,
        args.frag_error,
        "" if args.code is None else f"_{args.code}",
    )

    fig.colorbar(
        cm.ScalarMappable(cmap=cmap), ax=axs, shrink=0.3, orientation="horizontal"
    )
    fig.suptitle(f"Positions of disulfide bridges in {args.protein}", fontsize=24)
    fig.savefig(
        image_path,
    )
