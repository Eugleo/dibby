from dataclasses import dataclass
from typing import List, Optional, Tuple

from pyteomics.mass import calculate_mass

from src.utilities.constants import H2, SULPHUR
from src.utilities.error import within_bounds


# TODO: Shouldn't they be more positive? We're missing the H2 by default, no?
@dataclass(frozen=True, unsafe_hash=True)
class Modification:
    description: str
    mass: float


def combine_modifications(
    modifications: List[List[Optional[Modification]]],
    starting_mass: float,
    target_mass: float,
    error_ppm: float = 10,
) -> List[List[Modification]]:
    result = []

    def go(i: int, current_mass: float, selection: Tuple[Modification, ...] = ()):
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


IAA_ALKYLATION = Modification("Cys Alkylation (IAA)", 57.0214)
IAA_PAIR_ALKYLATION = Modification("2x Cys Alkylation (IAA)", IAA_ALKYLATION.mass * 2)
CYS_BOND = Modification("Disulphide Bond (–H2)", -calculate_mass(formula="H2"))
MET_OXIDATION = Modification("Met Oxidation", 15.9949)


def SEVERED_CYS_BOND_BOTH(i: int, j: int):
    return Modification(f"Cys {i}, Cys {j}: R–SSH + ()–R or R–SH + S=R", -H2)


def SEVERED_CYS_BOND_1(residue: int):
    return Modification(f"Cys {residue}: R–SH", 0)


SEVERED_CYS_BOND_MIN_MASS = -(SULPHUR + H2)
SEVERED_CYS_BOND_MAX_MASS = SULPHUR


def SEVERED_CYS_BOND_2(residue: int):
    return Modification(f"Cys {residue}: R–()", -(SULPHUR + H2))


def SEVERED_CYS_BOND_3(residue: int):
    return Modification(f"Cys {residue}: R=S", -H2)


def SEVERED_CYS_BOND_4(residue: int):
    return Modification(f"Cys {residue}: R–SSH", SULPHUR)
