from collections import Counter
from dataclasses import dataclass
from typing import Dict, Tuple, List, Iterator
from pyteomics.mass import calculate_mass

from src.model.modification import Modification


def trypsin(protein):
    last = 0
    result = []
    for i in range(len(protein) - 1):
        if protein[i] in ["K", "R"] and protein[i + 1] != "P":
            result.append((last, i + 1))
            last = i + 1
    result.append((last, len(protein)))
    return result


def pepsin(protein):
    last = 0
    result = []
    for i in range(len(protein) - 1):
        if protein[i] in ["F", "L", "W", "Y", "A", "E", "Q"]:
            result.append((last, i + 1))
            last = i + 1
    result.append((last, len(protein)))
    return result


@dataclass(init=False, frozen=True)
class Residue:
    name: str
    mass: float

    def __init__(self, name: str, modifications=None):
        if modifications is None:
            modifications = []

        object.__setattr__(self, "name", name)
        m = (
            calculate_mass(sequence=name)
            - calculate_mass(formula="H2O")
            + sum(m.mass for m in modifications)
        )
        object.__setattr__(self, "mass", m)


class Peptide:
    beginning: int
    end: int
    sequence: str
    min_mass: float
    mid_mass: float
    max_mass: float

    modifications: Dict[str, Tuple[Modification, int]]
    _residues: List[Residue]
    _residue_counts: Dict[str, int] = None
    _mass = None
    _minmass = None
    _maxmass = None

    def __init__(
        self,
        beginning: int,
        end: int,
        sequence: str,
        modifications: Dict[str, Tuple[Modification, int]],
    ):
        self.beginning = beginning
        self.end = end
        self.sequence = sequence
        self._residues = [Residue(resname) for resname in sequence]
        self.modifications = modifications

        zwitterion_mass = calculate_mass(
            sequence=self.sequence, ion_type="M", charge=0
        ) - calculate_mass(formula="H2O")

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
                "Peptides aren't contiguous. Got {} + {} instead.".format(
                    (self.beginning, self.end), other.beginning, other.end
                )
            )

        merged_mods = self.modifications
        for target, (mod, c2) in other.modifications.items():
            mod2, c1 = merged_mods.setdefault(target, (mod, 0))
            if mod != mod2:
                msg = "Peptides have incompatible modifications. These two differ at {}"
                raise ValueError(msg.format(target))
            merged_mods[target] = mod, c1 + c2

        return Peptide(
            self.beginning, other.end, self.sequence + other.sequence, merged_mods
        )

    def count(self, amino_acid):
        if self._residue_counts is None:
            self._residue_counts = Counter(self.sequence)
        return self._residue_counts[amino_acid]

    @property
    def modifications_anywhere(self) -> Iterator[Tuple[Modification, int]]:
        return (x for x in self.modifications.values())

    @property
    def cysteines(self) -> Iterator[int]:
        return (i for i, res in enumerate(self._residues) if res.name == "C")

    def __repr__(self):
        return "Peptide(beginning={}, end={}, seq={}, modifications={})".format(
            self.beginning, self.end, self.sequence, self.modifications
        )
