import itertools
import os
import pickle

from protein import Protein, trypsin, pepsin
import measurement
from measurement import read_mgf

from multiprocessing import Pool
import tqdm


def from_cache(source, cache):
    if os.path.exists(cache):
        with open(cache, "rb") as f:
            print(f"Loading {cache} from cache")
            return pickle.load(f)
    else:
        things = source()
        with open(cache, "wb") as f:
            print(f"Computed and saving to cache at {cache}")
            pickle.dump(things, f)
        return things


multiplier = measurement.binary_multiplier(tolerance_ppm=10)

DIR = "../data/mgf"
FILE = "190318_LYS_AT_50x_05"
MGF_DATA = f"{DIR}/{FILE}.mgf"
PICKLE_DATA = f"{DIR}/{FILE}.pickle"
measurements = from_cache(lambda: list(read_mgf(MGF_DATA)), PICKLE_DATA)


def get_item(t):
    pepid, peptide, threshold = t
    result = []
    for i, m in enumerate(measurements):
        if (
            peptide.total_charge == m.charge
            and abs(peptide.total_mz - m.peptide_mz) <= threshold
        ):
            score = m.score_match(peptide, multiplier)
            if score > 0:
                result.append(
                    {
                        "peptide": peptide,
                        "peptide_id": pepid,
                        "measurement": m,
                        "measurement_id": i,
                        "measurement_scan": m.scan,
                        "score": score,
                        "prec_error_ppm": 1e6
                        * (peptide.total_mz - m.peptide_mz)
                        / peptide.total_mz,
                    }
                )
    return result


if __name__ == "__main__":
    MYST = "KKLTKEGAAALCKMKHLADKVAKERSQELKDRTQNFAGYIEFELYRIDYWLEKLNGPKGRKDGYAKLSDSDIEKVKEIFNKAKDGITKQLPEAKKAGEEAGKLHTEVKKAAENARGQDLDDDTAKSTGLYRVLNWYCITKEERHNATPNCDGIQFRKHYLSVNRSAIDCSSTSYEENYDWSANALQVALNSWEDVKPKKLESAGSDKNCNIGQSSESHPCTMTEEWQTPYKETVEKLRELEDAYQRGKKAHDAMLGYANTAYAVNTKVEQEKPLTTGLEVLFQGPSAEPEA"
    OVA = "GSIGAASMEFCFDVFKELKVHHANENIFYCPIAIMSALAMVYLGAKDSTRTQINKVVRFDKLPGFGDSIEAQCGTSVNVHSSLRDILNQITKPNDVYSFSLASRLYAEERYPILPEYLQCVKELYRGGLEPINFQTAADQARELINSWVESQTNGIIRNVLQPSSVDSQTAMVLVNAIVFKGLWEKAFKDEDTQAMPFRVTEQESKPVQMMYQIGLFRVASMASEKMKILELPFASGTMSMLVLLPDEVSGLEQLESIINFEKLTEWTSSNVMEERKIKVYLPRMKMEEKYNLTSVLMAMGITDVFSSSANLSGISSAESLKISQAVHAAHAEINEAGREVVGSAEAGVDAASVSEEFRADHPFLFCIKHIATNAVLFFGRCVSP"
    LYS = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"

    PROTEIN = (LYS, "lysin")
    ENZYME = (trypsin, "trypsin")
    CHARGES = ({1, 2, 3, 4, 5, 6}, "1–6")
    MISSED_CLEAVAGES = 7

    protein = Protein(LYS)
    peptides = from_cache(
        lambda: list(
            protein.digest(ENZYME[0], charges=CHARGES[0], maxskip=MISSED_CLEAVAGES)
        ),
        f"../data/peptides/{PROTEIN[1]}_{ENZYME[1]}_ch{CHARGES[1]}_mc{MISSED_CLEAVAGES}.pickle",
    )

    #    with open(f"../out/{FILE}_scores.pickle", "wb") as f:
    #        for t in tqdm.tqdm(enumerate(peptides), total=len(peptides)):
    #            r = get_item(t)
    #            if r:
    #                pickle.dump(r, f)

    with Pool() as pool:
        with open(f"../out/{FILE}_scores.pickle", "wb") as f:
            for r in tqdm.tqdm(
                pool.imap_unordered(
                    get_item,
                    (
                        (i, pep, (10000 / 1e6) * pep.total_mz)
                        for i, pep in enumerate(peptides)
                    ),
                    chunksize=100,
                ),
                total=len(peptides),
            ):
                if r:
                    pickle.dump(r, f)
