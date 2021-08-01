from typing import List, Dict, Tuple, Optional

from src.model.peptide import Peptide, Residue
from src.model.modification import Modification, IAA_ALKYLATION


class Variant:
    _segments: List[Peptide]
    bonds: List[Tuple[int, int]]
    _disulfide_bond: Dict[int, int]
    _modifications: Dict[str, Tuple[Modification, int]]
    _alkylation_mass: float
    _residues: List[Residue]

    def __init__(
        self,
        segments: List[Peptide],
        disulfide_bonds: List[Tuple[int, int]],
        modifications: Dict[str, Tuple[Modification, int]],
        alkylation_mass: float = 57.0214,
    ):
        self._segments = sorted(segments, key=lambda segment: segment.beginning)
        self.segments = len(segments)
        self.bonds = list(disulfide_bonds)
        self._disulfide_bond = {}
        for c1, c2 in disulfide_bonds:
            self._disulfide_bond[c1] = c2
            self._disulfide_bond[c2] = c1
        self._modifications = modifications
        self._alkylation_mass = alkylation_mass

        self._residues = []

        for s in self._segments:
            for i in s:
                res = s[i]
                if res.name == "C" and self._disulfide_bond.get(i, None) is None:
                    self._residues.append(Residue("C", [IAA_ALKYLATION]))
                else:
                    self._residues.append(res)

    def to_dict(self):
        return {"var_bonds": self.bonds}

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

    def modification_on(self, residue: str) -> (Modification, int):
        return self._modifications[residue]

    def can_be_modified(self, residue: int) -> bool:
        return self[residue].name in self._modifications

    def segment(self, residue: int) -> int:
        for i, s in enumerate(self._segments):
            if s.beginning <= residue < s.end:
                return i

    def segment_bounds(self, segment: int) -> Tuple[int, int]:
        return self.segment_beginning(segment), self.segment_end(segment)

    def segment_beginning(self, segment: int) -> int:
        return self._segments[segment].beginning

    def segment_end(self, segment: int) -> int:
        return self._segments[segment].end

    def __str__(self):
        return "+".join(s.sequence for s in self._segments)
