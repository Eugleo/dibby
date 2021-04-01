import re
from typing import List, Set

import pepfrag as p
import numpy as np


def trypsin(protein):
    last = 0
    result = []
    for i in range(len(protein) - 1):
        if protein[i] in ["K", "R"] and protein[i + 1] != "P":
            result.append((last, i + 1))
            last = i + 1
    result.append((last, len(protein)))
    return result


class Fragment:
    def __init__(self, mz, name=None, intensity=None):
        self.mz = mz
        self.name = name

        # Neutral loss isn't always loss, sometimes it is addition
        m = re.match(r"\[?(([byM])(\d*)([-+]([^]]+))?)\]?\[(\d?)\+\]", name)
        if m:
            self.type = m.group(2)
            self.index = int(m.group(3)) if m.group(3) != "" else None
            self.neutral_loss = m.group(5)
            self.charge = int(m.group(6)) if m.group(6) != "" else 1
        else:
            raise ValueError(f"Incorrect fragment name: {name}")

        self.intensity = intensity

    def __str__(self):
        out = {
            "type": self.type,
            "index": self.index,
            "charge": self.charge,
            "neutral_loss": self.neutral_loss,
        }
        return f"<{self.__class__.__name__} {out}>"

    def __repr__(self):
        out = {
            "type": self.type,
            "index": self.index,
            "charge": self.charge,
            "neutral_loss": self.neutral_loss,
        }
        return f"<{self.__class__.__name__} {out}>"


class Modification:
    def __init__(self, name, mass, amino_acid):
        self.name = name
        self.mass = mass
        self.amino_acid = amino_acid

    def gen_modsites(self, protein_seq):
        return [
            p.ModSite(self.mass, i, self.name)
            for i, aa in enumerate(protein_seq)
            if aa == self.amino_acid
        ]


ALKYLATION_IAA = Modification("cysteine alkylation", 57.0214, "C")
ALKYLATION_EM = Modification("cysteine alkylation", 125.0477, "C")
OXIDATION_MET = Modification("methionine exodation", 15.9949, "M")


class Peptide(p.Peptide):
    def __init__(
        self, protein_seq, beginning, end, charge, modifications: List[Modification]
    ):
        if beginning >= end:
            raise AssertionError("Beginning should be strictly less than end")
        self.charge = charge
        self.seq = protein_seq[beginning:end]
        self.beginning = beginning
        self.end = end
        self.mods = [
            item for mod in modifications for item in mod.gen_modsites(self.seq)
        ]
        self.modstr = ", ".join(m.name for m in modifications)

        self.mass_type = p.MassType.mono
        self.radical = False
        self._fragments = None
        self._noncysteine_fragments = None
        self._fragment_masses = None
        self._noncysteine_fragment_masses = None

    def __contains__(self, item):
        return item in self.seq

    def __len__(self):
        return len(self.seq)

    @property
    def noncysteine_fragments(self):
        if self._noncysteine_fragments is None:
            if "C" not in self.seq:
                self._noncysteine_fragments = self.fragment()
            else:
                cysteines = np.nonzero(np.array(list(self.seq)) == "C")
                fc, lc = np.min(cysteines), np.max(cysteines)
                frags = self.fragment()
                self._noncysteine_fragments = [
                    f
                    for f in frags
                    if (f.type == "b" and f.index <= fc)
                    or (f.type == "y" and f.index > lc)
                ]
        return self._noncysteine_fragments

    def fragment(self) -> List[Fragment]:
        if self._fragments is None:
            self._fragments = [
                Fragment(mz, name=name)
                for mz, name, _ in p.Peptide.fragment(
                    self,
                    ion_types={
                        p.IonType.b: [],
                        p.IonType.y: [],
                    },
                )
            ]
            # self._fragments = [f for f in self._fragments if f.charge == 1]
        return self._fragments

    @property
    def noncysteine_fragment_masses(self):
        if self._noncysteine_fragment_masses is None:
            frags = self.noncysteine_fragments
            self._noncysteine_fragment_masses = np.fromiter(
                (f.mz for f in frags), np.float32, count=len(frags)
            )
        return self._noncysteine_fragment_masses

    @property
    def fragment_masses(self):
        if self._fragment_masses is None:
            frags = self.fragment()
            self._fragment_masses = np.fromiter(
                (f.mz for f in frags), np.float32, count=len(frags)
            )
        return self._fragment_masses


class Protein:
    def __init__(self, sequence, charge=2):
        self.sequence = sequence
        self.charge = charge

    def peptides(self, enzyme, maxskip=2) -> Set[Peptide]:
        peptides = enzyme(self.sequence)

        return {
            Peptide(
                self.sequence,
                beginning=peptides[i][0],
                end=peptides[i + s][1],
                charge=self.charge,
                modifications=m,
            )
            for s in range(maxskip)
            for i in range(len(peptides) - s)
            for m in [[ALKYLATION_IAA], [ALKYLATION_IAA, OXIDATION_MET]]
        }

    def __str__(self):
        return self.sequence
