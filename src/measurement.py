import protein
import os
import re

from pyteomics import mgf, mass
import numpy as np
from protein import Peptide

PROTON_MASS = mass.calculate_mass(formula="H").conjugate()


class PeptideMeasurement:
    def __init__(
        self,
        scan,
        time,
        charge,
        peptide_mz,
        peptide_intensity,
        peptide_mass_estimate,
        fragments_mz,
        fragments_intensity,
    ):
        self.scan = scan
        self.time = time
        self.charge = charge
        self.peptide_mz = peptide_mz
        self.peptide_intensity = peptide_intensity
        self.peptide_mass_estimate = peptide_mass_estimate
        nonzero = np.nonzero(fragments_intensity)[0]
        indices = np.argsort(fragments_mz[nonzero])
        self.fragments_mz = fragments_mz[nonzero][indices]
        self.fragments_intensity = fragments_intensity[nonzero][indices]

    def to_dicts(self):
        for mz, intenzity in zip(self.fragments_mz, self.fragments_intensity):
            yield {
                "scan": self.scan,
                "time": self.time,
                "charge": self.charge,
                "peptide_mz": self.peptide_mz,
                "peptide_intensity": self.peptide_intensity,
                "peptide_mass_estimate": self.peptide_mass_estimate,
                "fragment_mz": mz,
                "fragment_intensity": intenzity,
            }

    def contains(self, mz, tolerance=0.001):
        check = (self.fragments_mz >= mz - tolerance) & (
            self.fragments_mz <= mz + tolerance
        )
        return len(np.where(check)) > 0

    def score_match(self, peptide: Peptide, tolerance=0.02):
        frags_mz = peptide.fragment_masses(skip_cysteine=True)
        if len(frags_mz) == 0:
            return np.nan
        indices = np.searchsorted(self.fragments_mz, frags_mz, side="right")
        valid_indices = indices < len(self.fragments_mz)
        count = np.count_nonzero(
            (np.abs(self.fragments_mz[indices[valid_indices]] - frags_mz[valid_indices]) <= tolerance)
        )
        return count / len(frags_mz)


def read_mgf(path):
    """
    returns (scan ID, time, charge, mz, mass estimate)
    """
    with mgf.read(path) as reader:
        for i in reader:
            scan = int(re.match(".* scan=([0-9]+)", i["params"]["title"])[1])
            time = i["params"]["rtinseconds"]
            chargelist = i["params"]["charge"]
            if len(chargelist) > 1:
                raise AssertionError("ChargeList length>1 unsupported")
            charge = int(chargelist[0])

            peptide_mz = i["params"]["pepmass"][0]
            peptide_intensity = i["params"]["pepmass"][1]
            peptide_mass_estimate = peptide_mz * charge - charge * PROTON_MASS

            fragments_mz = i["m/z array"]
            fragments_intensity = i["intensity array"]

            yield PeptideMeasurement(
                scan,
                time,
                charge,
                peptide_mz,
                peptide_intensity,
                peptide_mass_estimate,
                fragments_mz,
                fragments_intensity,
            )
