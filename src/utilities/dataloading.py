import pickle

from pyteomics import fasta

from src.model.modification import MET_OXIDATION
from src.model.peptide import trypsin, Peptide


def load_precursor_matches(protein, kind, segments, error, code, debug=False):
    data_path = (
        "../out/precursor_matches/{}_{}_segments={}_error={}ppm{}.pickle".format(
            protein,
            kind,
            segments,
            error,
            "" if code is None else f"_{code}",
        )
    )

    print(f"Loading precursors from {data_path}...")

    with open(data_path, "rb") as f:
        return pickle.load(f)


def load_fragment_matches(
    protein: str,
    kind: str,
    segments: int,
    breaks: int,
    error: float,
    code: str,
):
    base = "../out/fragment_matches/{}_{}_segments={}_breaks={}_error={}ppm{}.pickle"
    data_path = base.format(
        protein,
        kind,
        segments,
        breaks,
        error,
        "" if code is None else f"_{code}",
    )

    print(f"Loading fragments from {data_path}")

    with open(data_path, "rb") as f:
        return pickle.load(f)


def load_protein(protein: str):
    seq_path = f"../data/fasta/{protein}.fasta"
    return [r.sequence for r in fasta.read(seq_path)][0]


def cleave_protein(protein: str):
    protein = load_protein(protein)
    peptides = []
    for b, e in trypsin(protein):
        seq = protein[b:e]
        met_ox = MET_OXIDATION, sum(aa == "M" for aa in seq)
        mods = {"M": met_ox} if "M" in seq else {}
        peptides.append(Peptide(b, e, seq, modifications=mods))
    return peptides
