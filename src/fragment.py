from dataclasses import dataclass
from typing import List, Tuple

from src.modification import Modification


@dataclass(frozen=True)
class Fragment:
    id: int
    sequence: str
    residue_ranges: List[Tuple[int, int]]

    intensity: float
    intensity_ratio: float

    mass: float
    target_mass: float
    mz: float
    target_mz: float
    charge: int
    break_count: int
    error_ppm: float
    modifications: List[Modification]
    connected_bonds: List[Tuple[int, int]]
    disconnected_cys: List[int]

    def to_dict(self):
        return {
            "frag_id": self.id,
            "frag_sequence": self.sequence,
            "frag_residue_ranges": self.residue_ranges,
            "frag_charge": self.charge,
            "frag_mass": self.mass,
            "frag_mz": self.mz,
            "frag_break_count": self.break_count,
            "frag_error_ppm": self.error_ppm,
            "frag_mods": [m.description for m in self.modifications],
            "frag_connected_bonds": self.connected_bonds,
            "frag_disconnected_cys": self.disconnected_cys,
            "frag_interesting_disconnected_cys": [
                c
                for c in self.disconnected_cys
                # TODO:Update to constant Mod
                if f"{c}: –SH" not in self.modifications
            ],
            "frag_intensity": self.intensity,
            "frag_intensity_ratio": self.intensity_ratio,
            "target_mass": self.target_mass,
            "target_mz": self.target_mz,
        }
