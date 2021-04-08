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
            pf.ModSite(self.mass, i + 1, self.name)
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
        pf.IonType.precursor: ["H2O", "NH3", "CO2"],
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
        self._total_mz = None
        self._total_charge = None
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
    def total_mz(self):
        if self._total_mz is None:
            mod = (len(self.peptides) - 1) * mass.calculate_mass(formula="H2")
            total_mass = sum(p.mass for p in self.peptides) - mod
            self._total_mz = (total_mass / self.total_charge) + mass.calculate_mass(
                formula="H"
            )
        return self._total_mz

    @property
    def total_charge(self):
        if self._total_charge is None:
            self._total_charge = sum(p.charge for p in self.peptides)
        return self._total_charge

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

    def digest(self, enzyme, maxskip=2, charges={2, 3, 4, 5, 6}) -> Set[Peptide]:
        peptides_initial = enzyme(self.sequence)
        peptides_skip = {
            (peptides_initial[i][0], peptides_initial[i + s][1])
            for s in range(maxskip + 1)
            for i in range(len(peptides_initial) - s)
        }

        def variants(peptide):
            seq = self.sequence[peptide[0] : peptide[1]]
            alkylations = ALKYLATION_IAA.everywhere(seq)
            return [
                {"charge": ch, "seq": seq, "mods": alkylations + mods}
                for ch in charges
                for mods in OXIDATION_MET.somewhere(seq)
            ]

        cysteines = [i for i, aa in enumerate(self.sequence) if aa == "C"]
        peptide_variants = {p: variants(p) for p in peptides_skip}
        peptides_having_cysteine = {
            c: [p for p in peptides_skip if p[0] <= c < p[1]] for c in cysteines
        }

        peptides_joined = [
            [v] for pepvars in peptide_variants.values() for v in pepvars
        ]
        for beg, end in itertools.combinations(cysteines, 2):
            for p1 in peptides_having_cysteine[beg]:
                for p2 in peptides_having_cysteine[end]:
                    for v1 in peptide_variants[p1]:
                        for v2 in peptide_variants[p2]:
                            if v1["charge"] + v2["charge"] in charges:
                                peptides_joined.append([v1, v2])

        return [Peptide(pepvars) for pepvars in peptides_joined]

    def __str__(self):
        return self.sequence
