import os
import re
from glob import glob
from typing import Iterable

from pyteomics import mgf, mass
import csv
from measurement import Measurement

PROTON_MASS = mass.calculate_mass(formula="H").conjugate()


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

            yield Measurement(
                scan,
                time,
                charge,
                peptide_mz,
                peptide_intensity,
                peptide_mass_estimate,
                fragments_mz,
                fragments_intensity,
            )


def to_csv(data: Iterable[Measurement], outpath):
    with open(outpath, mode="w") as outfile:
        cols = [
            "scan",
            "time",
            "charge",
            "peptide_mz",
            "peptide_intensity",
            "peptide_mass_estimate",
            "fragment_mz",
            "fragment_intensity",
        ]
        writer = csv.DictWriter(outfile, cols)
        writer.writeheader()

        for d in data:
            writer.writerows(d.to_dicts())


if __name__ == "__main__":
    for path in glob("../data/mgf/*.mgf"):
        data = read_mgf(path)
        pref, _ = os.path.splitext(path)
        _, name = os.path.split(pref)
        if not os.path.exists("../data/csv"):
            os.mkdir("../data/csv")
        to_csv(data, "../data/csv/" + name + ".csv")
