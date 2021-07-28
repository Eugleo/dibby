import pickle
import pprint
from typing import Tuple, List, Optional

from pyteomics import mass

from measurement import read_mgf, Scan
from protein import trypsin
from src.common import err_margin, combine_modifications, compute_error
from src.peptide import (
    Peptide,
    Mod,
    IAA_ALKYLATION,
    IAA_PAIR_ALKYLATION,
    MET_OXIDATION,
    CYS_BOND,
)


def match_precursors(
    peptides: List[Peptide],
    scan: Scan,
    max_segments: int,
    alkylation_mod: Mod = IAA_ALKYLATION,
    error_ppm: float = 10,
) -> List:
    target_mass = scan.prec_mass
    H2O = mass.calculate_mass(formula="H2O")
    H2 = mass.calculate_mass(formula="H2")
    H = mass.calculate_mass(formula="H")

    result = []

    def go(
        i: int,
        segments_left: int,
        selection: Tuple[int, ...],
        min_mass: float = H2O,
        base_mass: float = H2O,
        max_mass: float = H2O,
        free_cys_count: int = 0,
        waiting_for_cys: bool = False,
    ) -> None:
        has_alkylated_cys = free_cys_count % 2 == 1
        max_other_bonds = free_cys_count // 2
        min_realistic_mass = (
            min_mass
            + alkylation_mod.mass * has_alkylated_cys
            + (CYS_BOND.mass * max_other_bonds)
        )
        lower_bound = min_realistic_mass - err_margin(min_realistic_mass, error_ppm)

        if not waiting_for_cys:
            max_realistic_mass = max_mass + alkylation_mod.mass * free_cys_count
            upper_bound = max_realistic_mass + err_margin(max_realistic_mass, error_ppm)

            if lower_bound <= target_mass <= upper_bound:
                ranges = list(zip(selection[::2], (selection + (i,))[1::2]))
                possible_mods: List[List[Mod]] = []

                for b, e in ranges:
                    for p in peptides[b:e]:
                        for m, count in p.modifications_anywhere:
                            possible_mods += [[m, None]] * count

                for _ in range(max_other_bonds):
                    possible_mods.append([IAA_PAIR_ALKYLATION, CYS_BOND])

                if has_alkylated_cys:
                    # One Cys has to be alkylated, because it can't be in a bond
                    possible_mods.append([IAA_ALKYLATION])

                mod_combinations = combine_modifications(
                    possible_mods,
                    starting_mass=base_mass,
                    target_mass=target_mass,
                    error_ppm=error_ppm,
                )

                if mod_combinations:
                    segments = (
                        "".join(p.sequence for p in peptides[b:e]) for b, e in ranges
                    )
                    sequence = "+".join(segments)

                    for modifications in mod_combinations:
                        total_mass = base_mass + sum(m.mass for m in modifications)
                        joining_bonds = (max_segments - segments_left) - 1
                        other_bonds = sum(m == CYS_BOND for m in modifications)
                        mcs = [(e - b) - 1 for b, e in ranges]

                        result.append(
                            {
                                "scan": scan,
                                "scan_id": scan.id,
                                "scan_mass": scan.prec_mass,
                                "prec_sequence": sequence,
                                "prec_segments": (max_segments - segments_left),
                                "prec_tryptide_ranges": ranges,
                                "prec_residue_ranges": [
                                    (peptides[start].beginning, peptides[end - 1].end)
                                    for start, end in ranges
                                ],
                                "prec_max_mc_count": max(mcs),
                                "prec_mc": mcs,
                                "prec_cys_bond_count": other_bonds + joining_bonds,
                                "prec_mass": total_mass + scan.prec_charge * H,
                                "prec_mz": (total_mass / scan.prec_charge) + H,
                                "prec_error": compute_error(total_mass, target_mass),
                                "prec_mods": [m.description for m in modifications],
                            }
                        )

                    return

        if i == len(peptides) or lower_bound > target_mass:
            # Either we're out of peptides to add
            # Or our mass is too high beyond repair
            return

        if (
            not waiting_for_cys
            and min(segments_left, free_cys_count) > 0
            and selection[-1] != i  # Can't end if we just started
        ):
            # End current run, begin the next one
            for next_beginning in range(i, len(peptides)):
                go(
                    next_beginning,
                    segments_left=segments_left - 1,
                    selection=selection + (i, next_beginning),
                    min_mass=min_mass + (H2O - H2),
                    base_mass=base_mass + (H2O - H2),
                    max_mass=max_mass + (H2O - H2),
                    free_cys_count=free_cys_count - 1,
                    waiting_for_cys=True,
                )

        # Add current peptide
        new_free_cys = peptides[i].count("C") - waiting_for_cys
        go(
            i + 1,
            segments_left=segments_left,
            selection=selection,
            min_mass=min_mass + peptides[i].min_mass,
            base_mass=base_mass + peptides[i].mid_mass,
            max_mass=max_mass + peptides[i].max_mass,
            free_cys_count=free_cys_count + max(new_free_cys, 0),
            waiting_for_cys=new_free_cys < 0,
        )

    for beginning in range(0, len(peptides)):
        go(beginning, segments_left=max_segments - 1, selection=(beginning,))

    return result


def write_matched_precursors(
    protein: str, kind: str, max_segments: int, error_ppm: float, code: Optional[str]
):
    data_path = f"../data/mgf/190318_{protein}_{kind}_50x_05.mgf"
    seq_path = f"../data/fasta/{args.protein}.fasta"
    output_path = (
        "../out/precursor_matches/{}_{}_segments={}_error={}ppm{}.pickle".format(
            protein, kind, max_segments, error_ppm, "" if code is None else f"_{code}"
        )
    )
    protein = [r.sequence for r in fasta.read(seq_path)][0]
    peptides = []
    for b, e in trypsin(protein):
        seq = protein[b:e]
        met_ox = (MET_OXIDATION, sum(aa == "M" for aa in seq))
        mods = {"M": met_ox} if "M" in seq else {}
        peptides.append(Peptide(b, e, seq, modifications=mods))

    print(f"Loading scans from {data_path}")
    data = list(read_mgf(data_path))

    print(f"Saving matched precursors to {output_path}")
    with open(output_path, "wb") as f:
        for scan in tqdm.tqdm(data):
            precursors = match_precursors(
                peptides,
                scan=scan,
                alkylation_mod=IAA_ALKYLATION,
                max_segments=max_segments,
                error_ppm=error_ppm,
            )
            for precursor in precursors:
                pickle.dump(precursor, f)


if __name__ == "__main__":
    import argparse
    import tqdm
    from pyteomics import fasta

    args = argparse.ArgumentParser(description="Save precursor matches for given scans")

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
        help="measurement type (AT/RAT)",
    )
    args.add_argument(
        "--error",
        type=int,
        required=True,
        help="allowed measurement error in ppm",
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
    args = args.parse_args()
    write_matched_precursors(
        protein=args.protein,
        kind=args.kind,
        max_segments=args.segments,
        error_ppm=args.error,
        code=args.code,
    )
