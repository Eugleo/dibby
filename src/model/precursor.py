import dataclasses
from typing import List, Tuple

from src.model.peptide import Peptide
from src.model.modification import Modification, MET_OXIDATION
from src.model.variant import Variant


@dataclasses.dataclass(frozen=True)
class Precursor:
    sequence: str
    mass: float
    mz: float
    segments: List[Tuple[int, int]]
    residue_ranges: List[Tuple[int, int]]
    cys_bond_count: int
    alkylation_count: int
    modifications: List[Modification]
    error_ppm: float

    def __hash__(self):
        return (
            self.mass,
            self.mz,
            tuple(self.residue_ranges),
            tuple(self.segments),
            self.cys_bond_count,
            self.alkylation_count,
            self.error_ppm,
            tuple(self.modifications),
        ).__hash__()

    def to_dict(self):
        mcs = [(e - b) - 1 for b, e in self.segments]
        return {
            "prec_sequence": self.sequence,
            "prec_segment_count": len(self.segments),
            "prec_tryptide_ranges": self.segments,
            "prec_residue_ranges": self.residue_ranges,
            "prec_max_mc_count": max(mcs),
            "prec_mc": mcs,
            "prec_cys_bond_count": self.cys_bond_count,
            "prec_mass": self.mass,
            "prec_mz": self.mz,
            "prec_error": self.error_ppm,
            "prec_alkylation_count": self.alkylation_count,
            "prec_mods": [m.description for m in self.modifications],
        }

    def variants(self, tryptides: List[Peptide]) -> List[Variant]:
        segments: List[Peptide] = []
        for sb, se in self.segments:
            segment = tryptides[sb]
            for p in tryptides[sb + 1 : se]:
                segment += p
            segments.append(segment)
        total_bonds = self.cys_bond_count

        if total_bonds < len(segments) - 1:
            raise ValueError("There's more segments than bonds + 1.")

        cysteines = []
        for i, s in enumerate(segments):
            cysteines += [(c + s.beginning, i) for c in s.cysteines]

        result = []
        # TODO: Add support for other types of modifications
        met_oxs = sum(m == MET_OXIDATION for m in self.modifications)
        mods = {"M": (MET_OXIDATION, met_oxs)}

        for bonds in Precursor._gen_bonds(
            tuple(cysteines), total_bonds, len(segments) - 1
        ):
            result.append(Variant(segments, bonds, mods))

        return result

    @staticmethod
    def _gen_bonds(cysteines: Tuple[Tuple[int, int], ...], bonds: int, segments: int):
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
                    next_seg_connections = segment_connections + (
                        (current_seg, next_seg),
                    )
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
