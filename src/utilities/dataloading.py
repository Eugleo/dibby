import pickle

from pyteomics import fasta

from src.model.modification import MET_OXIDATION
from src.model.peptide import trypsin, Peptide


def load_precursor_matches(protein, kind, segments, error, code):
    data_path = (
        "../out/precursor_matches/{}_{}_segments={}_error={}ppm{}.pickle".format(
            protein,
            kind,
            segments,
            error,
            "" if code is None else f"_{code}",
        )
    )

    with open(data_path, "rb") as f:
        precursor_matches = []
        try:
            while True:
                precursor_matches.append(pickle.load(f))
        finally:
            return precursor_matches


def load_fragment_matches(
    protein: str, kind: str, segments: int, breaks: int, error: float, code: str
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

    with open(data_path, "rb") as f:
        fragment_matches = []
        try:
            while True:
                fragment_matches.append(pickle.load(f))
        finally:
            return fragment_matches


def cleave_protein(protein: str):
    seq_path = f"../data/fasta/{protein}.fasta"
    protein = [r.sequence for r in fasta.read(seq_path)][0]
    peptides = []
    for b, e in trypsin(protein):
        seq = protein[b:e]
        met_ox = MET_OXIDATION, sum(aa == "M" for aa in seq)
        mods = {"M": met_ox} if "M" in seq else {}
        peptides.append(Peptide(b, e, seq, modifications=mods))
    return peptides
