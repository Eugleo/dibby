import re

from pyteomics import mgf, mass, mzid
import numpy as np

from src.constants import PROTON


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

    def __str__(self):
        return (
            f"Scan(id={self.id}, nth={self.nth_in_order}, precursor={self.prec_mass})"
        )

    def to_dict(self):
        return {
            "scan_id": self.id,
            "scan_nth_in_order": self.nth_in_order,
            "scan_time": self.time,
            "scan_total_intensity": self.total_intensity,
            "prec_charge": self.prec_charge,
        }


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
            peptide_mass_estimate = peptide_mz * charge - charge * PROTON

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
