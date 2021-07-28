import protein
import os
import re

from pyteomics import mgf, mass, mzid
import numpy as np
from protein import Peptide

PROTON_MASS = mass.calculate_mass(formula="H").conjugate()


def binary_multiplier(tolerance_ppm=10):
    def multiplier(reference, measurement):
        bound = (tolerance_ppm / 1e6) * reference
        return np.where(np.abs(reference - measurement) <= bound, 1, 0)

    return multiplier


def linear_decay_multiplier(soft_ppm, hard_ppm):
    def multiplier(reference, measurement):
        soft = (soft_ppm / 1e6) * reference
        hard = (hard_ppm / 1e6) * reference
        dif = np.abs(reference - measurement)
        return np.where(
            dif <= soft, 1, 1 - (np.minimum(dif, hard) - soft) / (hard - soft)
        )

    return multiplier


class Scan:
    def __init__(
        self,
        nth_in_order,
        id,
        time,
        charge,
        prec_mz,
        prec_intensity,
        prec_mass,
        fragments_mz,
        fragments_intensity,
    ):
        self.nth_in_order = nth_in_order
        self.id = id
        self.time = time
        self.prec_charge = charge
        self.prec_mz = prec_mz
        self.prec_intensity = prec_intensity
        self.prec_mass = prec_mass
        nonzero = fragments_intensity >= 0.01 * np.max(fragments_intensity)
        indices = np.argsort(fragments_mz[nonzero])

        self.fragments_mz = fragments_mz[nonzero][indices]
        self.fragments_intensity = fragments_intensity[nonzero][indices]
        self.total_intensity = np.sum(self.fragments_intensity)

        self.masses_padded = np.pad(
            self.fragments_mz, 1, constant_values=(-np.inf, np.inf)
        )
        self.intensities_padded = np.pad(self.fragments_intensity, 1, constant_values=0)

    def to_dicts(self):
        for mz, intenzity in zip(self.fragments_mz, self.fragments_intensity):
            yield {
                "scan": self.id,
                "time": self.time,
                "charge": self.prec_charge,
                "peptide_mz": self.prec_mz,
                "peptide_intensity": self.prec_intensity,
                "peptide_mass_estimate": self.prec_mass,
                "fragment_mz": mz,
                "fragment_intensity": intenzity,
            }

    def contains(self, mz, tolerance=0.001):
        check = (self.fragments_mz >= mz - tolerance) & (
            self.fragments_mz <= mz + tolerance
        )
        return len(np.where(check)) > 0

    def score_match(
        self, peptide: Peptide, multiplier=binary_multiplier(tolerance_ppm=10)
    ):
        if peptide.ismerged:
            generated = peptide.fragment_masses
        else:
            generated = peptide.fragment_masses

        indices = np.searchsorted(self.masses_padded, generated)
        left = multiplier(generated, self.masses_padded[indices - 1])
        right = multiplier(generated, self.masses_padded[indices])
        matched_intensities = np.maximum(
            left * self.intensities_padded[indices - 1],
            right * self.intensities_padded[indices],
        )

        return np.sum(matched_intensities) / self.total_intensity


def read_mgf(path):
    """
    returns (scan ID, time, charge, mz, mass estimate)
    """
    with mgf.read(path) as reader:
        for id, i in enumerate(reader):
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

            yield Scan(
                id,
                scan,
                time,
                charge,
                peptide_mz,
                peptide_intensity,
                peptide_mass_estimate,
                fragments_mz,
                fragments_intensity,
            )


def read_mzid(path):
    return mzid.DataFrame(path)
