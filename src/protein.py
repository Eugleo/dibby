import itertools
import re
from typing import List, Set
from pyteomics import mass

import pepfrag as pf
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
        m = re.match(r"\[?(([byM])(\d*)([-+]([^]]+))?)]?\[(\d?)\+]", name)
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

    def everywhere(self, protein_seq):
        return [
            pf.ModSite(self.mass, i, self.name)
            for i, aa in enumerate(protein_seq)
            if aa == self.amino_acid
        ]

    def somewhere(self, protein_seq):
        locs = [i for i, aa in enumerate(protein_seq) if aa == self.amino_acid]
        return [
            [pf.ModSite(self.mass, i, self.name) for i in c]
            for n in range(len(locs) + 1)
            for c in itertools.combinations(locs, n)
        ]


ALKYLATION_IAA = Modification("iaa_c", 57.0214, "C")
ALKYLATION_EM = Modification("em_c", 125.0477, "C")
OXIDATION_MET = Modification("ox_m", 15.9949, "M")


class Peptide:
    _default_ion_types = {
        pf.IonType.b: ["NH3", "H2O"],
        pf.IonType.y: ["NH3", "H2O"],
    }

    def __init__(self, peps):
        self.peptides = [
            pf.Peptide(pep["seq"], pep["charge"], pep["mods"]) for pep in peps
        ]

        self._fragments = None
        self._noncysteine_fragments = None
        self._fragment_masses = None
        self._noncysteine_fragment_masses = None
        self._total_mass = None
        self._modstr = None

    def __contains__(self, item):
        return any(item in p.seq for p in self.peptides)

    def __str__(self):
        return "+".join(p.seq for p in self.peptides)

    @property
    def modstr(self):
        if self._modstr is None:
            self._modstr = str([(m.mod, m.site) for p in self.peptides for m in p.mods])
        return self._modstr

    @property
    def total_mass(self):
        if self._total_mass is None:
            mod = (len(self.peptides) - 1) * mass.calculate_mass(formula="H2")
            self._total_mass = sum(p.mz for p in self.peptides) - mod
        return self._total_mass

    @property
    def ismerged(self):
        return len(self.peptides) > 1

    @property
    def noncysteine_fragments(self):
        if self._noncysteine_fragments is None:
            if all("C" not in p.seq for p in self.peptides):
                self._noncysteine_fragments = self.fragment()
            else:
                self._noncysteine_fragments = []
                for pep in self.peptides:
                    cysteines = [i for i, aa in enumerate(pep.seq) if aa == "C"]
                    fc, lc = np.min(cysteines), np.max(cysteines)
                    for mz, name, _ in pep.fragment(self._default_ion_types):
                        f = Fragment(mz, name=name)
                        if (f.type == "b" and f.index <= fc) or (
                            f.type == "y" and f.index > lc
                        ):
                            self._noncysteine_fragments.append(f)
        return self._noncysteine_fragments

    def fragment(
        self,
        ion_types=None,
    ) -> List[Fragment]:
        if ion_types is None:
            ion_types = self._default_ion_types
        if self._fragments is None:
            self._fragments = [
                Fragment(mz, name=name)
                for pep in self.peptides
                for mz, name, _ in pf.Peptide.fragment(
                    pep,
                    ion_types=ion_types,
                )
            ]
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

    def digest(self, enzyme, maxskip=2) -> Set[Peptide]:
        peptides_initial = enzyme(self.sequence)
        peptides_skip = {
            (peptides_initial[i][0], peptides_initial[i + s][1])
            for s in range(maxskip + 1)
            for i in range(len(peptides_initial) - s)
        }

        cysteines = [i for i, aa in enumerate(self.sequence) if aa == "C"]

        peptides = set()
        peptides_by_c = {}
        for c in cysteines:
            for b, e in peptides_skip:
                if b <= c < e:
                    seq = self.sequence[b:e]
                    for mods in OXIDATION_MET.somewhere(seq):
                        pep = {
                            "seq": seq,
                            "charge": 2,
                            "mods": ALKYLATION_IAA.everywhere(seq) + mods,
                        }
                        peptides_by_c.setdefault(c, []).append(pep)
                        peptides.add(Peptide([pep]))

        for beg, end in itertools.combinations(cysteines, 2):
            for p1 in peptides_by_c[beg]:
                for p2 in peptides_by_c[end]:
                    peptides.add(Peptide([p1, p2]))

        return peptides

    def __str__(self):
        return self.sequence
