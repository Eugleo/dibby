import math
import pickle
from typing import Dict, Tuple, List, Any

import tqdm

from src.utilities.constants import (
    PROTON,
    H2,
    SULPHUR,
    H2O,
    OH,
    NH3,
)
from src.utilities.dataloading import load_precursor_matches, cleave_protein
from src.utilities.error import err_margin, compute_error
from src.model.modification import (
    combine_modifications,
    Modification,
    SEVERED_CYS_BOND_BOTH,
    SEVERED_CYS_BOND_1,
    SEVERED_CYS_BOND_2,
    SEVERED_CYS_BOND_3,
    SEVERED_CYS_BOND_4,
    SEVERED_CYS_BOND_MIN_MASS,
    SEVERED_CYS_BOND_MAX_MASS,
)
from src.model.fragment import Fragment
from src.model.peptide import (
    Peptide,
)
from src.model.precursor import Precursor
from typing import NamedTuple

from src.model.scan import Scan
import sys

from src.model.variant import Variant

sys.setrecursionlimit(10000)


def sort_into(x, xs: Tuple):
    return tuple(sorted(xs + (x,)))


Y_ION_MOD = PROTON
B_ION_MOD = -PROTON

H2O_NEUTRAL_LOSS = Modification("-H2O (neutral loss)", -H2O)
NH3_NEUTRAL_LOSS = Modification("-NH3 (neutral loss)", -NH3)

SegmentCuts = Dict[int, Tuple[Tuple[int, int], ...]]
Target = NamedTuple(
    "Target",
    [
        ("id", int),
        ("mass", float),
        ("mz", float),
        ("intensity", float),
        ("charge", int),
        ("i_ratio", float),
    ],
)


def first(xs, pred):
    return next(filter(pred, xs), None)


# TODO: Neumím matchovat věci z variant, kde jdou ty segmenty hned po sobě, a to jen u segmentu vpravo
def _fragments_matching_targets(
    targets: List[Target],
    variant: Variant,
    allowed_breaks: int,
    error_ppm: float,
) -> List[Fragment]:
    target_masses = [t.mass for t in targets]
    result = []

    def continue_run(
        i: int,
        target_i: int,
        min_end: int,
        max_end: int,
        current_mass: float,
        breaks_left: int,
        pivots: Tuple[int, ...],
        selection: Tuple[int, ...],
        disconnected_cys: Tuple[int, ...],
        connected_cys: Tuple[int, ...],
        neutral_losses_count: int,
        segment_cuts: SegmentCuts,
        fragment_start: int,
        modded_residues: Dict[str, int],
    ):
        if i > max_end:
            raise AssertionError("This should never happen, i > max_end")

        if target_i == len(targets):
            # We've ran out of targets to fulfil, we're done
            return

        # We can end the run
        if min_end <= i <= max_end:
            segment = variant.segment(i)
            # We are ending the run even though we don't need to
            premature_end = i < max_end
            open_end = i < variant.segment_end(segment)

            if not premature_end or breaks_left > 0:
                new_current_mass = current_mass + (B_ION_MOD if open_end else OH)
                new_neutral_loss_count = neutral_losses_count + open_end
                # TODO: Add support for negative-mass modifications
                # TODO: Tighten the lower bound
                # - add must-have modifications
                # - properly count Cys modifications
                min_mass = (
                    new_current_mass
                    + (len(disconnected_cys)) * SEVERED_CYS_BOND_MIN_MASS
                    + new_neutral_loss_count * H2O_NEUTRAL_LOSS.mass
                )
                lower_bound = min_mass - err_margin(min_mass, error_ppm)

                next_target = first(
                    range(target_i, len(target_masses)),
                    lambda ti: target_masses[ti] > lower_bound,
                )

                # Our mass is too high, all the targets are below us
                if next_target is None:
                    return

                new_segment_cuts = segment_cuts.copy()
                new_cuts = []
                for cb, ce in segment_cuts[segment]:
                    if ce == max_end:  # The 'current' cut we're now bisecting
                        current_beg = selection[-1]
                        if cb < current_beg:
                            new_cuts.append((cb, current_beg))
                        if i < ce:
                            new_cuts.append((i, ce))
                    else:
                        new_cuts.append((cb, ce))
                new_segment_cuts[segment] = tuple(new_cuts)

                # End the run
                start_new_run(
                    segment_cuts=new_segment_cuts,
                    pivots=pivots,
                    breaks_left=breaks_left - premature_end,
                    fragment_start=fragment_start,
                    target_i=next_target,
                    current_mass=new_current_mass,
                    selection=selection + (i,),
                    modded_residues=modded_residues,
                    disconnected_cys=disconnected_cys,
                    connected_cys=connected_cys,
                    neutral_losses_count=new_neutral_loss_count,
                )

            if not premature_end:
                # No point in continuing this run, there's nowhere to continue to
                return

        residue = variant[i]

        if variant.can_be_modified(i):
            new_modded_residues = modded_residues.copy()
            new_modded_residues.setdefault(residue.name, 0)
            new_modded_residues[residue.name] += 1
        else:
            new_modded_residues = modded_residues

        in_bond = (j := variant.bond_partner(i)) is not None and residue.name == "C"
        connected = in_bond and j in connected_cys
        disconnected = in_bond and j in disconnected_cys
        bond_is_new = in_bond and not (j in connected_cys or j in disconnected_cys)

        if bond_is_new:
            # We create two branches: one where the bond is kept...
            if j > fragment_start:
                # ...else we've already seen this fragment from the symmetric side
                continue_run(
                    i + 1,
                    target_i=target_i,
                    min_end=min_end,
                    max_end=max_end,
                    current_mass=current_mass + residue.mass - H2,
                    breaks_left=breaks_left,
                    pivots=sort_into(j, pivots),
                    selection=selection,
                    disconnected_cys=disconnected_cys,
                    connected_cys=connected_cys + (i,),
                    neutral_losses_count=neutral_losses_count,
                    segment_cuts=segment_cuts,
                    fragment_start=fragment_start,
                    modded_residues=modded_residues,
                )

            # ...and one where it's not
            if breaks_left > 0:
                continue_run(
                    i + 1,
                    target_i=target_i,
                    min_end=min_end,
                    max_end=max_end,
                    current_mass=current_mass + residue.mass,
                    breaks_left=breaks_left - 1,
                    pivots=pivots,
                    selection=selection,
                    disconnected_cys=disconnected_cys + (i,),
                    connected_cys=connected_cys,
                    neutral_losses_count=neutral_losses_count,
                    segment_cuts=segment_cuts,
                    fragment_start=fragment_start,
                    modded_residues=modded_residues,
                )
        else:
            # Add current residue, continue the run
            # If there was a bond into us, add us to broken or unbroken cysteines
            continue_run(
                i + 1,
                target_i=target_i,
                min_end=min_end,
                max_end=max_end,
                current_mass=current_mass + residue.mass,
                breaks_left=breaks_left,
                pivots=pivots,
                selection=selection,
                disconnected_cys=(disconnected_cys + (i,))
                if disconnected
                else disconnected_cys,
                connected_cys=(connected_cys + (i,)) if connected else connected_cys,
                neutral_losses_count=neutral_losses_count,
                segment_cuts=segment_cuts,
                fragment_start=fragment_start,
                modded_residues=new_modded_residues,
            )

    def start_new_run(
        segment_cuts: SegmentCuts,
        pivots: Tuple[int, ...],
        breaks_left: int,
        fragment_start: int,
        target_i: int = 0,
        current_mass: float = 0,
        selection: Tuple[int, ...] = (),
        modded_residues=None,
        disconnected_cys: Tuple[int, ...] = (),
        connected_cys=(),
        neutral_losses_count: int = 0,
    ):

        if modded_residues is None:
            modded_residues = {}

        if len(pivots) == 0:
            potential_mods = []

            for res, seen in modded_residues.items():
                peptide_mod, must_have = variant.modification_on(res)

                minimum_mods = max(must_have - (variant.count(res) - seen), 0)

                for _ in range(minimum_mods):
                    potential_mods.append([peptide_mod])

                # How many can I have
                maximum_mods = min(must_have, seen)
                # Optional mods
                for _ in range(maximum_mods - minimum_mods):
                    potential_mods.append([None, peptide_mod])

            for _ in range(neutral_losses_count):
                # MAYBE: Make this more granular? Or ditch this altogether
                potential_mods.append([H2O_NEUTRAL_LOSS, NH3_NEUTRAL_LOSS, None])

            for c in disconnected_cys:
                symmetric = (
                    (j := variant.bond_partner(c)) is not None
                ) and j in disconnected_cys

                # Symmetry breaking
                if symmetric and c > j:
                    continue

                if symmetric:
                    potential_mods.append([SEVERED_CYS_BOND_BOTH(c, j)])
                else:
                    potential_mods.append(
                        [
                            SEVERED_CYS_BOND_1(c),
                            SEVERED_CYS_BOND_2(c),
                            SEVERED_CYS_BOND_3(c),
                            SEVERED_CYS_BOND_4(c),
                        ]
                    )

            residue_ranges = list(zip(selection[::2], selection[1::2]))
            sequence = []
            for rb, re in residue_ranges:
                s = ""
                for k in range(rb, re):
                    s += variant[k].name
                sequence.append(s)
            sequence = "+".join(sequence)

            max_mass = current_mass + len(disconnected_cys) * SEVERED_CYS_BOND_MAX_MASS
            upper_bound = max_mass + err_margin(max_mass, error_ppm)

            valid_targets = [
                i
                for i in range(target_i, len(target_masses))
                if target_masses[i] <= upper_bound
            ]

            for target_i in valid_targets:
                combinations = combine_modifications(
                    potential_mods,
                    starting_mass=current_mass,
                    target_mass=target_masses[target_i],
                    error_ppm=error_ppm,
                )
                target = targets[target_i]

                mod_combination_set = set(
                    tuple(sorted(c, key=lambda m: m.mass)) for c in combinations
                )

                for mod_combination in mod_combination_set:
                    total_mass = current_mass + sum(m.mass for m in mod_combination)
                    result.append(
                        Fragment(
                            id=target.id,
                            sequence=sequence,
                            residue_ranges=residue_ranges,
                            intensity=target.intensity,
                            intensity_ratio=target.i_ratio,
                            target_mass=target.mass,
                            mass=total_mass + target.charge * PROTON,
                            target_mz=target.mz,
                            mz=(total_mass / target.charge) + PROTON,
                            charge=target.charge,
                            break_count=allowed_breaks - breaks_left,
                            error_ppm=compute_error(
                                total_mass, target_masses[target_i]
                            ),
                            modifications=mod_combination,
                            connected_bonds=[
                                (i, variant.bond_partner(i))
                                for i in connected_cys
                                if i < variant.bond_partner(i)
                            ],
                            disconnected_cys=list(disconnected_cys),
                        )
                    )

            return

        pivot = pivots[0]
        segment = variant.segment(pivot)

        applicable_cuts = [
            (cb, ce) for cb, ce in segment_cuts[segment] if cb <= pivot < ce
        ]
        if not applicable_cuts:
            return
        cut_beginning, cut_end = applicable_cuts[0]

        min_beginning = max(cut_beginning, fragment_start)
        max_beginning = pivot
        min_end = pivot + 1
        max_end = variant.segment_end(segment)

        segment_beginning = variant.segment_beginning(segment)
        first_in_segment = cut_beginning == segment_beginning
        shift_optim = min_beginning != pivot and not first_in_segment

        for next_beginning in range(min_beginning + shift_optim, max_beginning + 1):
            is_break = next_beginning > cut_beginning
            is_open = next_beginning > segment_beginning

            if not is_break or breaks_left > 0:
                continue_run(
                    next_beginning,
                    target_i=target_i,
                    min_end=min_end,
                    max_end=max_end,
                    current_mass=current_mass + (Y_ION_MOD if is_open else PROTON),
                    breaks_left=breaks_left - is_break,
                    pivots=pivots[1:],
                    selection=selection + (next_beginning,),
                    disconnected_cys=disconnected_cys,
                    connected_cys=connected_cys,
                    neutral_losses_count=neutral_losses_count + is_open,
                    segment_cuts=segment_cuts,
                    fragment_start=fragment_start,
                    modded_residues=modded_residues,
                )

    for b in variant:
        start_new_run(
            segment_cuts={
                s: (variant.segment_bounds(s),) for s in range(variant.segments)
            },
            pivots=(b,),
            breaks_left=allowed_breaks,
            fragment_start=b,
        )

    return result


def write_matched_fragments(
    precursor_matches: List[Dict],
    tryptides: List[Peptide],
    output_path: str,
    max_allowed_breaks: int,
    error_ppm: float,
):
    print(f"Writing the matched fragments to {output_path}")
    fragment_matches = []

    for precursor_match in tqdm.tqdm(precursor_matches):
        scan: Scan = precursor_match["scan"]
        precursor: Precursor = precursor_match["precursor"]
        total_intensity = sum(scan.fragments_intensity)

        # ADD FOR BSA
        # if precursor.mass > 5500 or precursor.to_dict()["prec_max_mc_count"] > 3:
        #    continue

        targets = []
        for frag_id, (frag_mz, frag_intensity) in enumerate(
            zip(scan.fragments_mz, scan.fragments_intensity)
        ):
            coef = scan.prec_charge * PROTON + precursor.mass
            max_charge = min(
                scan.prec_charge,
                math.trunc(coef / (frag_mz - err_margin(frag_mz, error_ppm=error_ppm))),
                # 3,  # ADD FOR BSA
            )

            for ch in range(1, max_charge + 1):
                targets.append(
                    Target(
                        frag_id,
                        frag_mz * ch - PROTON * ch,
                        frag_mz,
                        frag_intensity,
                        ch,
                        frag_intensity / total_intensity,
                    )
                )

        targets = list(sorted(targets, key=lambda t: t.mass))

        variants = precursor.variants(tryptides)
        for variant in variants:
            fragments = _fragments_matching_targets(
                targets,
                variant,
                allowed_breaks=max_allowed_breaks,
                error_ppm=error_ppm,
            )

            base_info: Dict[str, Any] = {
                "scan": scan,
                "precursor": precursor,
                "variant": variant,
                "variant_count": len(variants),
            }

            if not fragments:
                fragment_matches.append(base_info | {"fragment": None})

            for f in fragments:
                fragment_matches.append(base_info | {"fragment": f})

    with open(output_path, "wb") as f:
        pickle.dump(fragment_matches, f)

    return fragment_matches


if __name__ == "__main__":
    import argparse

    args = argparse.ArgumentParser(
        description="Save fragment matches for given precursors"
    )

    # Add the arguments
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
        help="measurement kind (AT/RAT)",
    )
    args.add_argument(
        "--error",
        type=int,
        required=True,
        help="allowed measurement error in ppm",
    )
    args.add_argument(
        "--perror",
        type=int,
        required=True,
        help="allowed precursor measurement error in ppm",
    )
    args.add_argument(
        "--breaks",
        type=int,
        required=True,
        help="upper bound of breaks in matched fragments",
    )
    args.add_argument(
        "--segments",
        type=int,
        required=True,
        help="upper bound of segment count in matched precursors",
    )
    args.add_argument(
        "--code",
        type=str,
        default=None,
        required=False,
        help="code to append to the output file name",
    )

    # Execute the parse_args() method
    args = args.parse_args()
    output_template = (
        "../out/fragment_matches/{}_{}_segments={}_breaks={}_error={}ppm{}.pickle"
    )
    write_matched_fragments(
        precursor_matches=load_precursor_matches(
            args.protein, args.kind, args.segments, args.perror, args.code
        ),
        tryptides=cleave_protein(args.protein),
        output_path=output_template.format(
            args.protein,
            args.kind,
            args.segments,
            args.breaks,
            args.error,
            "" if args.code is None else f"_{args.code}",
        ),
        max_allowed_breaks=args.breaks,
        error_ppm=args.error,
    )

# Optimizations, measured on 845–3896
# - selective charge, 150 -> 120
# - caching modified Cys Residue, 120 -> 113
# - shift optimization, 113 -> 111
# - bundling charges to a multitarget search, 111 -> 55
# - bundling fragments to a multitarget search, 55 -> 10
