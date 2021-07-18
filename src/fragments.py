import math
import pickle
from typing import Dict, Tuple, Optional, List, Union

from pyteomics import mass
from common import LYS

from precursor import Mod, Peptide, Residue, err_margin, compute_error, within_bounds


# TODO: Optimize using binary search
from src.measurement import PeptideMeasurement
from src.protein import trypsin


def sort_into(x, xs: Tuple):
    return tuple(sorted(xs + (x,)))


class MultiP:
    _segments: List[Peptide]
    _disulfide_bond: Dict[int, int]
    _modifications: Dict[str, Tuple[Mod, int]]
    _alkylation_mass: float

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

    def __getitem__(self, index: int) -> Residue:
        for segment in self._segments:
            residue = segment[index]
            if residue is not None:
                if residue.name == "C" and self.bond_partner(index) is None:
                    return Residue(
                        residue.name, [Mod("Cys alkylation", self._alkylation_mass)]
                    )
                return residue

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


def fragments(target_mass, peptide: MultiP, allowed_breaks, ppm_error=10):

    result = []

    def go_run(
        i: int,
        min_end: int,
        max_end: int,
        current_mass: float,
        breaks_left: int,
        pivots: Tuple[int, ...],
        selection: Tuple[int, ...],
        broken_cysteines: Tuple[int, ...],
        unbroken_cysteines: Tuple[int, ...],
        neutral_losses_count: int,
        max_i_per_segment: Dict[int, int],
        fragment_start: int,
        modded_residues: Dict[str, int],
    ):
        # TODO: Add support for negative-mass modifications
        # TODO: Tighten the lower bound
        # - add must-have modifications
        # - properly count Cys modification
        min_possible_mass = (
            current_mass + B_ION_MOD + (len(broken_cysteines)) * (-H2 - SULPHUR)
        )
        lower_bound = min_possible_mass - err_margin(min_possible_mass, ppm_error)
        if lower_bound > target_mass:
            # Too heavy, beyond repair, end the whole branch
            return

        if i >= max_end:
            if i > max_end:
                raise AssertionError("This should never happen, i > max_end")

            # We can't grow any longer, end the run
            go(
                i,
                max_i_per_segment,
                current_mass + OH,
                breaks_left,
                broken_cysteines,
                unbroken_cysteines,
                neutral_losses_count,
                pivots,
                selection + (i,),
                fragment_start,
                modded_residues,
            )
            return
        else:
            residue = peptide[i]

            # This residue is (was) part of a disulfide bond
            if residue.name == "C" and ((j := peptide.bond_partner(i)) is not None):
                if j in broken_cysteines:
                    # This Cys (i) has a broken partner, so it has to be broken, too
                    go_run(
                        i + 1,
                        min_end,
                        max_end,
                        current_mass + residue.mass,
                        breaks_left,  # Already added when we were breaking j
                        pivots,
                        selection,
                        broken_cysteines + (i,),
                        unbroken_cysteines,
                        neutral_losses_count,
                        max_i_per_segment,
                        fragment_start,
                        modded_residues,
                    )
                elif j in unbroken_cysteines:
                    # We have already seen this bond
                    # We can't break it, and neither can we jump through it
                    # So, just add the Cys (i) and go on

                    go_run(
                        i + 1,
                        min_end,
                        max_end,
                        current_mass + residue.mass,
                        breaks_left,
                        pivots,
                        selection,
                        broken_cysteines,
                        unbroken_cysteines + (i,),
                        neutral_losses_count,
                        max_i_per_segment,
                        fragment_start,
                        modded_residues,
                    )
                else:
                    # We haven't seen this Cys (i) nor its bond partner yet
                    if breaks_left > 0:
                        # Break the bond
                        go_run(
                            i + 1,
                            min_end,
                            max_end,
                            current_mass + residue.mass,
                            breaks_left - 1,
                            pivots,
                            selection,
                            broken_cysteines + (i,),
                            unbroken_cysteines,
                            neutral_losses_count,
                            max_i_per_segment,
                            fragment_start,
                            modded_residues,
                        )

                    if j > fragment_start:
                        # Keep the bond, add new run
                        go_run(
                            i + 1,
                            min_end,
                            max_end,
                            # Subtract H2 for the bond
                            current_mass + residue.mass - H2,
                            breaks_left,
                            sort_into(j, pivots),
                            selection,
                            broken_cysteines,
                            unbroken_cysteines + (i,),
                            neutral_losses_count,
                            max_i_per_segment,
                            fragment_start,
                            modded_residues,
                        )

            else:
                # Add current residue, continue the run
                if peptide.can_be_modified(i):
                    new_modded_residues = modded_residues.copy()
                    new_modded_residues.setdefault(residue.name, 0)
                    new_modded_residues[residue.name] += 1
                else:
                    new_modded_residues = modded_residues

                go_run(
                    i + 1,
                    min_end,
                    max_end,
                    current_mass + residue.mass,
                    breaks_left,
                    pivots,
                    selection,
                    broken_cysteines,
                    unbroken_cysteines,
                    neutral_losses_count,
                    max_i_per_segment,
                    fragment_start,
                    new_modded_residues,
                )

                # Break this run and end it
                if i >= min_end and breaks_left > 0:
                    # End the run
                    go(
                        i,
                        max_i_per_segment,
                        current_mass + B_ION_MOD,
                        breaks_left - 1,
                        broken_cysteines,
                        unbroken_cysteines,
                        neutral_losses_count + 1,
                        pivots,
                        selection + (i,),
                        fragment_start,
                        modded_residues,
                    )
                    return

    def go(
        max_i: int,
        max_i_per_segment: Dict[int, int],
        current_mass: float,
        breaks_left: int,
        broken_cysteines: Tuple[int, ...],
        unbroken_cysteines,
        neutral_losses_count: int,
        pivots: Tuple[int, ...],
        selection: Tuple[int, ...],
        fragment_start: int,
        modded_residues: Dict[str, int],
    ):

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

            for c in broken_cysteines:
                symmetric = (
                    (j := peptide.bond_partner(c)) is not None
                ) and j in broken_cysteines

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

            combinations = combine_modifications_2(
                potential_mods,
                starting_mass=current_mass,
                target_mass=target_mass,
                ppm_error=ppm_error,
            )

            for modifications in combinations:
                total_mass = current_mass + sum(m.mass for m in modifications)
                result.append(
                    {
                        "seq": seq,
                        "ranges": ranges,
                        "mass": total_mass,
                        "error": compute_error(total_mass, target_mass),
                        "mods": modifications,
                    }
                )

            return

        segment = peptide.segment(max_i)
        new_max_i_per_segment = max_i_per_segment.copy()
        new_max_i_per_segment[segment] = max_i

        pivot = pivots[0]
        current_segment = peptide.segment(pivot)
        current_segment_max_i = new_max_i_per_segment[current_segment]

        beg_start = max(
            peptide.segment_beginning(segment), current_segment_max_i, fragment_start
        )
        beg_end = pivot

        end_start = pivot + 1
        end_end = peptide.segment_end(current_segment)

        at_segment_start = current_segment_max_i == peptide.segment_beginning(segment)
        shift_optim = (
            not at_segment_start
            and current_segment_max_i != pivot
            and not beg_start == fragment_start
        )

        # TODO: Add back shift optim
        for b in range(beg_start + 0, beg_end + 1):
            is_break = b > beg_start or (
                b == fragment_start and b > peptide.segment_beginning(current_segment)
            )

            if not is_break:
                go_run(
                    b,
                    end_start,
                    end_end,
                    current_mass + PROTON,
                    breaks_left,
                    pivots[1:],
                    selection + (b,),
                    broken_cysteines,
                    unbroken_cysteines,
                    neutral_losses_count,
                    new_max_i_per_segment,
                    fragment_start,
                    modded_residues,
                )

            if is_break and breaks_left > 0:
                go_run(
                    b,
                    end_start,
                    end_end,
                    current_mass + Y_ION_MOD,
                    breaks_left - 1,
                    pivots[1:],
                    selection + (b,),
                    broken_cysteines,
                    unbroken_cysteines,
                    neutral_losses_count + 1,
                    new_max_i_per_segment,
                    fragment_start,
                    modded_residues,
                )

    for b in peptide:
        pointers = {s: peptide.segment_beginning(s) for s in range(peptide.segments)}
        go(
            min(pointers.values()),
            pointers,
            0,
            allowed_breaks,
            (),
            (),
            0,
            (b,),
            (),
            b,
            modded_residues={},
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
            counter = 0
            while counter < 600:
                counter += 1

                if counter < 200:
                    continue

                match = pickle.load(f)
                measurement: PeptideMeasurement = match["measurement"]

                total_intensity = sum(measurement.fragments_intensity)

                if counter % 50 == 0:
                    print("Scan", measurement.scan, "count", counter)

                for precursor in match["matches"]:
                    for multiprotein, bonds in gen_multip(peptides, precursor):
                        multip_str = str(multiprotein)

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
                                matches = fragments(
                                    frag * ch - PROTON * ch,
                                    multiprotein,
                                    allowed_breaks=2,
                                )
                                for m in matches:
                                    to_print = {
                                        "scan": measurement.scan,
                                        "precursor_sequence": multip_str,
                                        "bonds": bonds,
                                        "fragment_mz": frag,
                                        "intensity": intensity,
                                        "match": m,
                                        "charge": ch,
                                    }
                                    of.write(f"{to_print}\n")
    end_time = time.time()
    print(f"It took {end_time - start_time} seconds")
