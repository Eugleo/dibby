import dataclasses
import pickle
from collections import Counter
from typing import Tuple, List, Dict, Union, Iterator
import time

from protein import trypsin
from measurement import read_mgf, PeptideMeasurement
from pyteomics import mass


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
    min_mass: float
    mid_mass: float
    max_mass: float

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

        self.min_mass = zwitterion_mass + neg
        self.mid_mass = zwitterion_mass
        self.max_mass = zwitterion_mass + pos

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
                if m is None:  # None = no modification
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
) -> List:
    target_mass = measurement.peptide_mass_estimate
    h2o = mass.calculate_mass(formula="H2O")
    h2 = mass.calculate_mass(formula="H2")

    result = []

    def go(
        i: int,
        segments_left: int,
        selection: Tuple[int, ...],
        min_mass: float = h2o,
        base_mass: float = h2o,
        max_mass: float = h2o,
        free_cys_count: int = 0,
        waiting_for_cys: bool = False,
    ) -> None:
        has_alkylated_cys = free_cys_count % 2 == 1
        min_realistic_mass = min_mass + alkylation_mass * has_alkylated_cys
        lower_bound = min_realistic_mass - err_margin(min_realistic_mass, error_ppm)

        if not waiting_for_cys:
            max_realistic_mass = max_mass + alkylation_mass * free_cys_count
            upper_bound = max_realistic_mass + err_margin(max_realistic_mass, error_ppm)

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
                        [Mod("Alkylated Cys Pair", alkylation_mass * 2), None]
                    )

                if has_alkylated_cys:
                    # One Cys has to be alkylated, because it can't be in a bond
                    possible_mods.append([Mod("Alkylated Cys", alkylation_mass)])

                mod_combinations = combine_modifications(
                    possible_mods,
                    starting_mass=base_mass,
                    target_mass=target_mass,
                    error_ppm=error_ppm,
                )

                if mod_combinations:
                    segments = (
                        "".join(p.seq for p in peptides[b:e]) for b, e in ranges
                    )
                    seq = "+".join(segments)

                    for modifications in mod_combinations:
                        total_mass = base_mass + sum(m.mass for m in modifications)

                        alkylated_pairs = sum(
                            m.description == "Alkylated Cys Pair" for m in modifications
                        )
                        joining_bonds = (max_segments - segments_left) - 1
                        other_bonds = max_other_bonds - alkylated_pairs

                        result.append(
                            {
                                "sequence": seq,
                                "ranges": ranges,
                                "missed_cleavages": max(e - b for b, e in ranges) - 1,
                                "cys_bonds": other_bonds + joining_bonds,
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
                    min_mass=min_mass + (h2o - h2),
                    base_mass=base_mass + (h2o - h2),
                    max_mass=max_mass + (h2o - h2),
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
        "--type",
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

    # Execute the parse_args() method
    args = args.parse_args()

    data_path = f"../data/mgf/190318_{args.protein}_{args.type}_50x_05.mgf"
    seq_path = f"../data/fasta/{args.protein}.fasta"
    output_path = f"../out/precursor_matches/{args.protein}_{args.type}_segments={args.segments}_error={args.error}ppm.pickle"

    measurements = {m.scan: m for m in read_mgf(data_path)}
    protein = [r.sequence for r in fasta.read(seq_path)][0]

    peptides = []
    for b, e in trypsin(protein):
        seq = protein[b:e]
        met_ox = (Mod("Met Oxidation", 15.9949), sum(aa == "M" for aa in seq))
        mods = {"M": met_ox} if "M" in seq else {}
        peptides.append(Peptide(b, e, seq, modifications=mods))

    start_time = time.time()
    with open(output_path, "wb") as f:
        for scan, measurement in tqdm.tqdm(measurements.items()):
            matches = match_precursors(
                peptides,
                measurement,
                alkylation_mass=57.0214,
                max_segments=args.segments,
                error_ppm=args.error,
            )
            for m in matches:
                pickle.dump({"measurement": measurement} | m, f)

    end_time = time.time()

    print(f"It took {end_time - start_time} seconds")
