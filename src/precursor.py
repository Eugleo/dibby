import dataclasses
import pickle
from collections import Counter
from typing import Tuple, List, Dict, Union, Set, Iterator, Optional
from enum import Enum, auto
import time

from protein import trypsin
from measurement import read_mgf, PeptideMeasurement
from pyteomics import mass
from common import LYS, BSA


@dataclasses.dataclass
class Mod:
    description: str
    mass: float

    def __hash__(self):
        return (self.description, self.mass).__hash__()


@dataclasses.dataclass
class Residue:
    name: str
    mass: float

    def __init__(self, name: str, modifications=None):
        if modifications is None:
            modifications = []

        self.name = name
        self.mass = (
            mass.calculate_mass(sequence=name)
            - mass.calculate_mass(formula="H2O")
            + sum(m.mass for m in modifications)
        )


@dataclasses.dataclass
class Range:
    minimum: float
    middelum: float
    maximum: float

    def add_constant(self, c: float):
        return Range(self.minimum + c, self.middelum + c, self.maximum + c)

    def from_constant(c: float):
        return Range(c, c, c)

    def __add__(self, other):
        return Range(
            self.minimum + other.minimum,
            self.middelum + other.middelum,
            self.maximum + other.maximum,
        )


class Peptide:
    beginning: int
    end: int
    seq: str
    mass_range: Range

    _modifications: Dict[str, Tuple[Mod, int]]
    _residues: List[Residue]
    _residue_counts: Dict[str, int] = None
    _mass = None
    _minmass = None
    _maxmass = None

    def __init__(
        self,
        beginning: int,
        end: int,
        seq: str,
        modifications: Dict[str, Tuple[Mod, int]],
    ):
        self.beginning = beginning
        self.end = end
        self.seq = seq
        self._residues = [Residue(resname) for resname in seq]
        self._modifications = modifications

        zwitterion_mass = mass.calculate_mass(
            sequence=self.seq, ion_type="M", charge=0
        ) - mass.calculate_mass(formula="H2O")

        pos, neg = 0, 0
        for m, count in modifications.values():
            if m.mass < 0:
                neg += m.mass * count
            else:
                pos += m.mass * count

        self.mass_range = Range(
            zwitterion_mass + neg, zwitterion_mass, zwitterion_mass + pos
        )

    def __getitem__(self, index: int):
        if self.beginning <= index < self.end:
            return self._residues[index - self.beginning]
        return None

    def __iter__(self):
        return range(self.beginning, self.end).__iter__()

    def __add__(self, other):
        if other.beginning != self.end:
            raise ValueError(
                f"Peptides can only be added when they are contiguous. Got {(self.beginning, self.end)} + {other.beginning, other.end} instead."
            )

        merged_mods = self._modifications
        for target, (mod, c2) in other._modifications.items():
            mod2, c1 = merged_mods.setdefault(target, (mod, 0))
            if mod != mod2:
                raise ValueError(
                    f"Peptides can only be added when they have compatible modifications. These two differ at {target}"
                )
            merged_mods[target] = mod, c1 + c2

        return Peptide(self.beginning, other.end, self.seq + other.seq, merged_mods)

    def count(self, amino_acid):
        if self._residue_counts is None:
            self._residue_counts = Counter(self.seq)
        return self._residue_counts[amino_acid]

    @property
    def modifications_anywhere(self) -> Iterator[Tuple[Mod, int]]:
        return (x for x in self._modifications.values())

    @property
    def cysteines(self) -> Iterator[int]:
        return (i for i, res in enumerate(self._residues) if res.name == "C")

    def __repr__(self):
        return f"Peptide(beginning={self.beginning}, end={self.end}, seq={self.seq}, modifications={self._modifications})"


def within_bounds(reference_mass, measured_mass, error_ppm: float = 10):
    return abs(reference_mass - measured_mass) <= err_margin(reference_mass, error_ppm)


def err_margin(reference_mass, error_ppm: float = 10):
    return (error_ppm / 1e6) * reference_mass


def compute_error(reference_mass, measured_mass):
    return 1e6 * abs(measured_mass - reference_mass) / reference_mass


# Pass None when you want to allow to skip a mod
def combine_modifications(
    modifications: List[List[Union[Mod, None]]],
    starting_mass: float,
    target_mass: float,
    error_ppm: float = 10,
) -> List[List[Mod]]:
    result = []

    def go(i: int, current_mass: float, selection: Tuple[Mod, ...] = ()):
        if i == len(modifications):
            if within_bounds(current_mass, target_mass, error_ppm):
                result.append(selection)
        else:
            for m in modifications[i]:
                if m is None:
                    go(i + 1, current_mass, selection)
                else:
                    go(i + 1, current_mass + m.mass, selection + (m,))

    go(0, current_mass=starting_mass)
    return list(set(result))


def match_precursors(
    peptides: List[Peptide],
    measurement: PeptideMeasurement,
    max_segments: int,
    alkylation_mass: float = 57.0214,
    error_ppm: int = 10,
) -> List[str]:
    target_mass = measurement.peptide_mass_estimate
    h2o = mass.calculate_mass(formula="H2O")
    h2 = mass.calculate_mass(formula="H2")

    result = []

    def go(
        i: int,
        segments_left: int,
        selection: Tuple[int, ...],
        current_mass: Range = Range.from_constant(h2o),
        free_cys_count: int = 0,
        waiting_for_cys: bool = False,
    ) -> None:
        has_alkylated_cys = free_cys_count % 2 == 1
        min_mass = current_mass.minimum + alkylation_mass * has_alkylated_cys
        lower_bound = min_mass - err_margin(min_mass, error_ppm)

        if not waiting_for_cys:
            max_mass = current_mass.maximum + alkylation_mass * free_cys_count
            upper_bound = max_mass + err_margin(max_mass, error_ppm)

            if lower_bound <= target_mass <= upper_bound:
                ranges = list(zip(selection[::2], (selection + (i,))[1::2]))
                possible_mods: List[List[Mod]] = []

                for b, e in ranges:
                    for p in peptides[b:e]:
                        for m, count in p.modifications_anywhere:
                            possible_mods += [
                                [Mod(m.description, m.mass), None]
                            ] * count

                max_other_bonds = free_cys_count // 2
                for _ in range(max_other_bonds):
                    possible_mods.append(
                        [Mod("cys_pair_alk", alkylation_mass * 2), None]
                    )

                if has_alkylated_cys:
                    # One Cys has to be alkylated, because it can't be in a bond
                    possible_mods.append([Mod("Alkylated Cys", alkylation_mass)])

                mod_combinations = combine_modifications(
                    possible_mods,
                    starting_mass=current_mass.middelum,
                    target_mass=target_mass,
                    error_ppm=error_ppm,
                )

                if mod_combinations:
                    segments = (
                        "".join(p.seq for p in peptides[b:e]) for b, e in ranges
                    )
                    seq = "+".join(segments)

                    for modifications in mod_combinations:
                        total_mass = current_mass.middelum + sum(
                            m.mass for m in modifications
                        )

                        alkylated_pairs = sum(
                            m.description == "Alkylated Cys Pair" for m in modifications
                        )
                        intra_bonds = max_other_bonds - alkylated_pairs
                        inter_bonds = (max_segments - segments_left) - 1

                        result.append(
                            {
                                "sequence": seq,
                                "ranges": ranges,
                                "cys_bonds": intra_bonds + inter_bonds,
                                "mass": total_mass,
                                "error": compute_error(total_mass, target_mass),
                                "mods": modifications,
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
            for beginning in range(i, len(peptides)):
                go(
                    beginning,
                    segments_left=segments_left - 1,
                    selection=selection + (i, beginning),
                    current_mass=current_mass.add_constant(h2o - h2),
                    free_cys_count=free_cys_count - 1,
                    waiting_for_cys=True,
                )

        # Add current peptide
        new_free_cys = peptides[i].count("C") - waiting_for_cys
        go(
            i + 1,
            segments_left=segments_left,
            selection=selection,
            current_mass=current_mass + peptides[i].mass_range,
            free_cys_count=free_cys_count + max(new_free_cys, 0),
            waiting_for_cys=new_free_cys < 0,
        )

    for beginning in range(0, len(peptides)):
        go(beginning, segments_left=max_segments - 1, selection=(beginning,))

    return result


if __name__ == "__main__":
    measurements = {m.scan: m for m in read_mgf("../data/mgf/190318_LYS_AT_50x_05.mgf")}

    peptides = []
    for b, e in trypsin(LYS):
        seq = LYS[b:e]
        met_ox = (Mod("met_ox", 15.9949), sum(aa == "M" for aa in seq))
        mods = {"M": met_ox} if "M" in seq else {}
        peptides.append(Peptide(b, e, seq, modifications=mods))

    FILE_PATH = "../out/precursor_matches_lys_at_2_inter_bonds.pickle"

    import tqdm

    start_time = time.time()
    with open(FILE_PATH, "wb") as f:
        for scan, measurement in tqdm.tqdm(measurements.items()):
            matches = match_precursors(
                peptides,
                measurement,
                alkylation_mass=57.0214,
                max_segments=3,
                error_ppm=10,
            )
            if matches:
                pickle.dump({"measurement": measurement, "matches": matches}, f)

    end_time = time.time()

    print(f"It took {end_time - start_time} seconds")
