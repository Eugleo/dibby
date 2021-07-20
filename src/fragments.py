import math
import pickle
from typing import Dict, Tuple, Optional, List, Union

import tqdm
from pyteomics import mass
from common import LYS

from precursor import Mod, Peptide, Residue, err_margin, compute_error, within_bounds

import sys

sys.setrecursionlimit(10000)

from src.measurement import PeptideMeasurement
from src.protein import trypsin


def sort_into(x, xs: Tuple):
    return tuple(sorted(xs + (x,)))


class MultiP:
    _segments: List[Peptide]
    _disulfide_bond: Dict[int, int]
    _modifications: Dict[str, Tuple[Mod, int]]
    _alkylation_mass: float
    _residues: List[Residue]

    def __init__(
        self,
        segments: List[Peptide],
        disulfide_bonds: List[Tuple[int, int]],
        modifications: Dict[str, Tuple[Mod, int]],
        alkylation_mass: float = 57.0214,
    ):
        self._segments = sorted(segments, key=lambda s: s.beginning)
        self.segments = len(segments)
        self._disulfide_bond = {}
        for c1, c2 in disulfide_bonds:
            self._disulfide_bond[c1] = c2
            self._disulfide_bond[c2] = c1
        self._modifications = modifications
        self._alkylation_mass = alkylation_mass

        self._residues = []

        for s in self._segments:
            b = s.beginning
            for i in s:
                res = s[i]
                if res.name == "C" and self._disulfide_bond.get(i, None) is None:
                    self._residues.append(
                        Residue("C", [Mod("Cys alkylation", self._alkylation_mass)])
                    )
                else:
                    self._residues.append(res)

    def __getitem__(self, index: int) -> Residue:
        index_trans = 0
        for segment in self._segments:
            if segment.beginning <= index < segment.end:
                return self._residues[index_trans + index - segment.beginning]
            else:
                index_trans += segment.end - segment.beginning

    def __iter__(self):
        for segment in self._segments:
            for i in segment:
                yield i

    def count(self, residue: str) -> int:
        return sum(s.count(residue) for s in self._segments)

    def bond_partner(self, residue: int) -> Optional[int]:
        return self._disulfide_bond.get(residue, None)

    def modification_on(self, residue: str) -> (Mod, int):
        return self._modifications[residue]

    def can_be_modified(self, residue: int) -> bool:
        return self[residue].name in self._modifications

    def segment(self, residue: int) -> int:
        for i, s in enumerate(self._segments):
            if s.beginning <= residue <= s.end:
                return i

    def segment_beginning(self, segment: int) -> int:
        return self._segments[segment].beginning

    def segment_end(self, segment: int) -> int:
        return self._segments[segment].end

    def __str__(self):
        return "+".join(s.seq for s in self._segments)


def combine_modifications_2(
    modifications: List[List[Union[Mod, None]]],
    starting_mass: float,
    target_mass: float,
    ppm_error: float = 10,
) -> List[List[Mod]]:
    result = []

    def go(i, current, selection):
        if i == len(modifications):
            if within_bounds(current, target_mass, ppm_error):
                result.append(selection)
        else:
            for m in modifications[i]:
                if m is None:
                    go(i + 1, current, selection)
                else:
                    go(i + 1, current + m.mass, selection + (m,))

    go(0, current=starting_mass, selection=())
    return list(set(result))


OH = mass.calculate_mass(formula="OH")
PROTON = mass.calculate_mass(formula="H")
H2 = mass.calculate_mass(formula="H2")
H2O = mass.calculate_mass(formula="H2O")
NH3 = mass.calculate_mass(formula="NH3")
SULPHUR = mass.calculate_mass(formula="S")
Y_ION_MOD = PROTON
B_ION_MOD = -PROTON


def fragments(
    target_masses: List[float], peptide: MultiP, allowed_breaks, ppm_error=10
):

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
        max_i_per_segment: Dict[int, int],
        fragment_start: int,
        modded_residues: Dict[str, int],
    ):
        if i > max_end:
            raise AssertionError("This should never happen, i > max_end")

        if target_i == len(target_masses):
            # We've ran out of targets to fulfil, we're done
            return

        # TODO: Add support for negative-mass modifications
        # TODO: Tighten the lower bound
        # - add must-have modifications
        # - properly count Cys modification
        min_possible_mass = (
            current_mass + B_ION_MOD + (len(disconnected_cys)) * (-H2 - SULPHUR)
        )
        lower_bound = min_possible_mass - err_margin(min_possible_mass, ppm_error)
        if lower_bound > target_masses[target_i]:
            # Jump to the next target, current one is beyond reach (too low)
            continue_run(
                i,
                target_i=target_i + 1,
                min_end=min_end,
                max_end=max_end,
                current_mass=current_mass,
                breaks_left=breaks_left,
                pivots=pivots,
                selection=selection,
                disconnected_cys=disconnected_cys,
                connected_cys=connected_cys,
                neutral_losses_count=neutral_losses_count,
                max_i_per_segment=max_i_per_segment,
                fragment_start=fragment_start,
                modded_residues=modded_residues,
            )

            return

        # We can end the run
        if min_end <= i <= max_end:
            # We are ending the run even though we don't need to
            premature_end = i < max_end

            current_segment = peptide.segment(i)
            new_max_i_per_segment = max_i_per_segment.copy()
            new_max_i_per_segment[current_segment] = i

            # End the run
            start_new_run(
                new_max_i_per_segment,
                pivots,
                breaks_left - premature_end,
                fragment_start,
                target_i,
                current_mass + (B_ION_MOD if premature_end else OH),
                selection + (i,),
                modded_residues,
                disconnected_cys,
                connected_cys,
                neutral_losses_count + premature_end,
            )

            if not premature_end:
                # No point in continuing this run, there's nowhere to continue to
                return

        residue = peptide[i]

        if peptide.can_be_modified(i):
            new_modded_residues = modded_residues.copy()
            new_modded_residues.setdefault(residue.name, 0)
            new_modded_residues[residue.name] += 1
        else:
            new_modded_residues = modded_residues

        in_bond = (j := peptide.bond_partner(i)) is not None and residue.name == "C"
        connected = in_bond and j in disconnected_cys
        bond_is_new = not connected and j not in disconnected_cys
        if in_bond and bond_is_new:
            # We create two branches: one where the bond is kept...
            if j > fragment_start:
                # ...else we've already seen this fragment from the symmetric side
                continue_run(
                    i + 1,
                    target_i,
                    min_end,
                    max_end,
                    current_mass + residue.mass - H2,
                    breaks_left,
                    sort_into(j, pivots),
                    selection,
                    disconnected_cys,
                    connected_cys + (i,),
                    neutral_losses_count,
                    max_i_per_segment,
                    fragment_start,
                    modded_residues,
                )

            # ...and one where it's not
            if breaks_left > 0:
                continue_run(
                    i + 1,
                    target_i,
                    min_end,
                    max_end,
                    current_mass + residue.mass,
                    breaks_left - 1,
                    pivots,
                    selection,
                    disconnected_cys + (i,),
                    connected_cys,
                    neutral_losses_count,
                    max_i_per_segment,
                    fragment_start,
                    modded_residues,
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
                if not connected
                else disconnected_cys,
                connected_cys=(connected_cys + (i,)) if connected else connected_cys,
                neutral_losses_count=neutral_losses_count,
                max_i_per_segment=max_i_per_segment,
                fragment_start=fragment_start,
                modded_residues=new_modded_residues,
            )

    def start_new_run(
        max_i_per_segment: Dict[int, int],
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

            final_mass = current_mass

            for res, seen in modded_residues.items():
                peptide_mod, must_have = peptide.modification_on(res)

                minimum_mods = max(must_have - (peptide.count(res) - seen), 0)

                for _ in range(minimum_mods):
                    potential_mods.append([peptide_mod])

                final_mass += minimum_mods * peptide_mod.mass

                # How many can I have
                maximum_mods = min(must_have, seen)
                # Optional mods
                for _ in range(maximum_mods - minimum_mods):
                    potential_mods.append([None, peptide_mod])

            for _ in range(neutral_losses_count):
                # MAYBE: Make this more granular? Or ditch this altogether
                potential_mods.append(
                    [
                        Mod("–H2O neutral loss", -H2O),
                        Mod("–NH3 neutral loss", -NH3),
                        None,
                    ]
                )

            for c in disconnected_cys:
                symmetric = (
                    (j := peptide.bond_partner(c)) is not None
                ) and j in disconnected_cys

                # Symmetry breaking
                if symmetric and c > j:
                    continue

                if symmetric:
                    potential_mods.append([Mod("-SSH + () or -SH + =S", -H2)])
                else:
                    potential_mods.append(
                        [
                            Mod("–SSH", SULPHUR),
                            Mod("– ()", -(SULPHUR + H2)),
                            Mod("=S", -H2),
                            Mod("–SH", 0),
                        ]
                    )

            ranges = list(zip(selection[::2], selection[1::2]))

            seq = []
            for b, e in ranges:
                s = ""
                for k in range(b, e):
                    s += peptide[k].name
                seq.append(s)
            seq = "+".join(seq)

            # TODO: Make the 190 threshold selectively lower
            # It's basically the mass of the next added residue
            valid_targets = [
                i
                for i, t in enumerate(target_masses)
                if target_masses[target_i] <= t <= target_masses[target_i] + 85
            ]

            for target_i in valid_targets:
                combinations = combine_modifications_2(
                    potential_mods,
                    starting_mass=current_mass,
                    target_mass=target_masses[target_i],
                    ppm_error=ppm_error,
                )

                for modifications in combinations:
                    total_mass = current_mass + sum(m.mass for m in modifications)
                    result.append(
                        (
                            target_i,
                            {
                                "seq": seq,
                                "ranges": ranges,
                                "mass": total_mass,
                                "error": compute_error(
                                    total_mass, target_masses[target_i]
                                ),
                                "mods": modifications,
                            },
                        )
                    )

            return

        pivot = pivots[0]
        segment = peptide.segment(pivot)
        current_segment_beginning = peptide.segment_beginning(segment)
        current_segment_max_i = max_i_per_segment[segment]

        beg_start = max(
            current_segment_max_i,
            fragment_start,
        )
        beg_end = pivot

        end_start = pivot + 1
        end_end = peptide.segment_end(segment)

        first_in_segment = max_i_per_segment[segment] == current_segment_beginning
        shift_optim = max_i_per_segment[segment] != pivot and not first_in_segment

        for b in range(beg_start + shift_optim, beg_end + 1):
            is_break = b > beg_start or (
                b == fragment_start and b > current_segment_beginning
            )

            if not is_break or breaks_left > 0:
                continue_run(
                    b,
                    target_i=target_i,
                    min_end=end_start,
                    max_end=end_end,
                    current_mass=current_mass + (Y_ION_MOD if is_break else PROTON),
                    breaks_left=breaks_left - is_break,
                    pivots=pivots[1:],
                    selection=selection + (b,),
                    disconnected_cys=disconnected_cys,
                    connected_cys=connected_cys,
                    neutral_losses_count=neutral_losses_count + is_break,
                    max_i_per_segment=max_i_per_segment,
                    fragment_start=fragment_start,
                    modded_residues=modded_residues,
                )

    for b in peptide:
        pointers = {s: peptide.segment_beginning(s) for s in range(peptide.segments)}
        start_new_run(
            pointers, pivots=(b,), breaks_left=allowed_breaks, fragment_start=b
        )

    return result


def gen_bonds(cysteines: Tuple[Tuple[int, int], ...], bonds: int, segments: int):
    result = []

    def go(
        cysteines_left: Tuple[Tuple[int, int], ...],
        bonds_left: int,
        segment_connections: Tuple[Tuple[int, int], ...],
        current_bonds: Tuple[Tuple[int, int], ...],
    ):
        if bonds_left == 0 and len(segment_connections) == segments:
            result.append(current_bonds)
            return

        if bonds_left == 0 or len(cysteines_left) == 0:
            return

        current_cys, current_seg = cysteines_left[0]
        others = cysteines_left[1:]

        # This Cys isn't in a bond
        go(others, bonds_left, segment_connections, current_bonds)

        # This Cys is in a bond
        for i, (next_cys, next_seg) in enumerate(others):
            if next_seg != current_seg:
                next_seg_connections = segment_connections + ((current_seg, next_seg),)
            else:
                next_seg_connections = segment_connections

            go(
                others[:i] + others[i + 1 :],
                bonds_left - 1,
                next_seg_connections,
                current_bonds + ((current_cys, next_cys),),
            )

    go(
        cysteines_left=cysteines,
        bonds_left=bonds,
        segment_connections=(),
        current_bonds=(),
    )

    return result


import copy


def gen_multip(
    peptides: List[Peptide], d
) -> List[Tuple[MultiP, List[Tuple[int, int]]]]:
    segments: List[Peptide] = []
    for b, e in d["ranges"]:
        segment = None
        for p in peptides[b:e]:
            if segment is None:
                segment = copy.deepcopy(p)
            else:
                segment += p
        segments.append(segment)

    total_bonds = d["cys_bonds"]

    if total_bonds < len(segments) - 1:
        raise ValueError("There's more segments than bonds + 1.")

    cysteines = []
    for i, s in enumerate(segments):
        cysteines += [(c + s.beginning, i) for c in s.cysteines]

    result = []
    mods = {
        "M": (
            Mod(description="Met Oxidation", mass=15.9949),
            sum(m.description == "Met Oxidation" for m in d["mods"]),
        )
    }

    for bonds in gen_bonds(tuple(cysteines), total_bonds, len(segments) - 1):
        result.append((MultiP(segments, bonds, mods), bonds))

    return result


if __name__ == "__main__":
    import time

    precursors_file = "../out/precursor_matches_lys_at_2_inter_bonds.pickle"
    # fragments_file = "../out/fragment_matches_lys_at_2_inter_bonds.pickle"
    fragments_file = "../out/fragments_matches.txt"

    peptides = []
    for b, e in trypsin(LYS):
        seq = LYS[b:e]
        met_ox = (Mod("met_ox", 15.9949), sum(aa == "M" for aa in seq))
        mods = {"M": met_ox} if "M" in seq else {}
        peptides.append(Peptide(b, e, seq, modifications=mods))

    start_time = time.time()

    with open(fragments_file, "w") as of:
        with open(precursors_file, "rb") as f:
            precursor_matches = []
            try:
                while True:
                    precursor_matches.append(pickle.load(f))
            finally:
                for match in tqdm.tqdm(precursor_matches[200:600]):
                    measurement: PeptideMeasurement = match["measurement"]
                    total_intensity = sum(measurement.fragments_intensity)

                    for precursor in match["matches"]:
                        targets = []
                        for frag, intensity in zip(
                            measurement.fragments_mz,
                            measurement.fragments_intensity,
                        ):
                            coef = measurement.charge * PROTON + precursor["mass"]
                            max_charge = min(
                                measurement.charge,
                                math.trunc(
                                    coef / (frag - err_margin(frag, error_ppm=10))
                                ),
                            )

                            for ch in range(1, max_charge + 1):
                                targets.append(
                                    (
                                        {
                                            "fragment_mz": frag,
                                            "charge": ch,
                                            "intensity": intensity,
                                            "total_intensity": total_intensity,
                                            "score": intensity / total_intensity,
                                        },
                                        frag * ch - PROTON * ch,
                                    ),
                                )
                        targets = sorted(targets, key=lambda t: t[1])

                        for multiprotein, bonds in gen_multip(peptides, precursor):
                            matches = fragments(
                                [t for _, t in targets],
                                multiprotein,
                                allowed_breaks=2,
                            )

                            multip_str = str(multiprotein)

                            matches = sorted(
                                matches,
                                # key=lambda m: m[0],
                                key=lambda m: (
                                    targets[m[0]][0]["fragment_mz"],
                                    targets[m[0]][0]["charge"],
                                ),
                            )

                            for i, m in matches:
                                to_print = {
                                    "measurement": measurement,
                                    "multipeptide": multiprotein,
                                    "bonds": bonds,
                                    "match": m,
                                } | targets[i][0]
                                of.write(f"{to_print}")
                                # pickle.dump(to_print, of)

                end_time = time.time()
                print(f"It took {end_time - start_time} seconds")

# Optimizations, measured on 845–3896
# - selective charge, 150 -> 120
# - caching modified Cys Residue, 120 -> 113
# - shift optimization, 113 -> 111
# - bundling charges to a multitarget search, 111 -> 55
# - bundling fragments to a multitarget search, 55 -> 10
