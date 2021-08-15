import re
import sys
from random import randrange, randint, choices, shuffle
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
from pepfrag import ModSite, IonType, pepfrag
from pyteomics.mass import calculate_mass

from src.fragment_matching import write_matched_fragments
from src.model.fragment import Fragment
from src.model.modification import IAA_ALKYLATION, CYS_BOND
from src.model.peptide import Peptide
from src.model.precursor import Precursor
from src.model.scan import Scan
from src.precursor_matching import write_matched_precursors
from src.utilities.constants import PROTON
from src.utilities.dataloading import cleave_protein


def intersects(t, u):
    x, y = t
    a, b = u
    return not (x >= b or y <= a)


def remove_peptide_duplicates(xs):
    return list(dict(tp) for tp in set(tuple(p.items()) for p in (xs)))


def connected_cys_count(prec):
    return sum(res == "C" for res in prec.sequence) - prec.alkylation_count


def generate_simple_peptides(
    tryptides: List[Peptide],
    cys_bond_tryptides,
    base_count=10_000,
    max_mc_count=5,
):
    peptide_without_bond_cys: List[Dict] = []
    peptide_with_bond_cys: List[Dict] = []
    for _ in range(0, base_count):
        b = randrange(0, len(tryptides) - 1)
        e = randrange(b + 1, min(len(tryptides), b + 2 + max_mc_count))
        if b < e:
            charge = randint(1, 5)
            sequence = "".join(t.sequence for t in tryptides[b:e])
            alkylations = sum(res == "C" for res in sequence)
            cys_overlap = [i for i in cys_bond_tryptides if i in range(b, e)]

            if cys_overlap:
                alkylations -= len(cys_overlap)
            mass = calculate_mass(sequence) + alkylations * IAA_ALKYLATION.mass

            prec: Dict = {
                "charge": charge,
                "precursor": Precursor(
                    sequence=sequence,
                    mass=mass,
                    mz=mass / charge + PROTON,
                    segments=[(b, e)],
                    residue_ranges=[(tryptides[b].beginning, tryptides[e - 1].end)],
                    cys_bond_count=0,
                    alkylation_count=alkylations,
                    modifications=[],
                    error_ppm=0,
                ),
            }

            if cys_overlap:
                peptide_with_bond_cys.append(prec)
            else:
                peptide_without_bond_cys.append(prec)

    return (
        remove_peptide_duplicates(peptide_without_bond_cys),
        remove_peptide_duplicates(peptide_with_bond_cys),
    )


def generate_dipeptides(peptides: List[Dict], max_charge=5):
    dipeptides = []
    for i, s in enumerate(peptides):
        prec: Precursor = s["precursor"]
        for t in peptides[i:]:
            qrec: Precursor = t["precursor"]
            if not intersects(prec.segments[0], qrec.segments[0]):
                charge = randint(1, max_charge)
                ps = sorted([prec, qrec], key=lambda p: p.segments[0][0])
                mass = prec.mass + qrec.mass + CYS_BOND.mass
                joined = Precursor(
                    sequence=ps[0].sequence + "+" + ps[1].sequence,
                    mass=mass,
                    mz=mass / charge + PROTON,
                    segments=ps[0].segments + ps[1].segments,
                    residue_ranges=ps[0].residue_ranges + ps[1].residue_ranges,
                    cys_bond_count=1,
                    alkylation_count=prec.alkylation_count + qrec.alkylation_count,
                    modifications=ps[0].modifications + ps[1].modifications,
                    error_ppm=0,
                )
                dipeptides.append({"charge": charge, "precursor": joined})
    return remove_peptide_duplicates(dipeptides)


def generate_unipeptides(peptides: List[Dict]):
    unipeptides = []
    for s in peptides:
        p: Precursor = s["precursor"]
        if connected_cys_count(p) == 2:
            charge = s["charge"]
            precursor = Precursor(
                p.sequence,
                p.mass + CYS_BOND.mass,
                (p.mass + CYS_BOND.mass) / charge + PROTON,
                p.segments,
                p.residue_ranges,
                p.cys_bond_count,
                p.alkylation_count,
                p.modifications,
                p.error_ppm,
            )
            unipeptides.append({"charge": charge, "precursor": precursor})
    return remove_peptide_duplicates(unipeptides)


def valid_frags(frags, cys, length):
    def ok(frag):
        if "b" in frag[1]:
            return frag[2] > cys
        else:
            return frag[2] >= (length - cys)

    return [f for f in frags if ok(f)]


def charge_from_code(code):
    match = re.match(r".*\[(\d+)?\+]$", code)
    if match.group(1) is None:
        return 1
    else:
        return int(match.group(1))


def split_on_simple_frags(seq, frags, cysteines):
    b, e = seq
    safe = []
    unsafe = []
    for f in frags:
        mass, code, i = f
        if "b" in code:
            if not any(b <= c < b + i for c in cysteines):
                safe.append(f)
                continue
        else:
            if not any(e - i <= c < e for c in cysteines):
                safe.append(f)
                continue
        unsafe.append(f)

    return safe, unsafe


def simple_fragment(
    id, sequence, residue_range, charge, mz, break_count, intensity_ratio, intensity=10
):
    return Fragment(
        id=id,
        sequence=sequence,
        residue_ranges=residue_range,
        intensity=intensity,
        intensity_ratio=intensity_ratio,
        target_mass=(mz - PROTON) * charge,
        mass=(mz - PROTON) * charge,
        target_mz=mz,
        mz=mz,
        charge=charge,
        break_count=break_count,
        error_ppm=0,
        modifications=[IAA_ALKYLATION for res in sequence if res == "C"],
        connected_bonds=[],
        disconnected_cys=[],
    )


def fragment_sequence(seq, frag, residue_range):
    _, code, i = frag
    sequence = seq[:i] if "b" in code else seq[-i:]
    b, e = residue_range
    frag_residue_range = (b, b + i) if "b" in code else (e - i, e)

    return sequence, frag_residue_range


def simple_frags_to_fragments(frags, prec_sequence, prec_residue_range, precursor):
    fragments = []
    for id, frag in enumerate(frags):
        mz, code, i = frag
        frag_charge = charge_from_code(code)
        frag_sequence, frag_residue_range = fragment_sequence(
            prec_sequence, frag, prec_residue_range
        )
        fragment = simple_fragment(
            id=id,
            sequence=frag_sequence,
            residue_range=[frag_residue_range],
            charge=frag_charge,
            mz=mz,
            break_count=int(prec_residue_range != frag_residue_range),
            intensity_ratio=1 / len(frags),
        )

        fragments.append(
            {"fragment": fragment, "precursor": precursor, "var_bonds": []}
        )
    return fragments


def safe_choose_n(fragments, n=50):
    return list(sorted(list(set(choices(fragments, k=n)))))


def pepfrag_fragments(
    sequence: str,
    residue_range: Tuple[int, int],
    charge: int,
    ion_types,
    bond_cys_res: List[int],
    count=50,
):
    frags = pepfrag.Peptide(
        sequence,
        charge=charge,
        modifications=[
            ModSite(IAA_ALKYLATION.mass, ri + 1, IAA_ALKYLATION.description)
            for ri, (ai, res) in enumerate(zip(sequence, range(*residue_range)))
            if res == "C" and ai not in bond_cys_res
        ],
    ).fragment(ion_types=ion_types)
    return safe_choose_n(frags, count)


def fragment_simple_peptide(
    peptide: Dict,
    bond_cys_res: List[int],
    count=50,
    ion_types=None,
):
    if ion_types is None:
        ion_types = {IonType.y: [], IonType.b: [], IonType.precursor: []}

    precursor: Precursor = peptide["precursor"]
    sequence = precursor.sequence
    residue_range = precursor.residue_ranges[0]

    # if connected_cys_count(precursor) == 0:
    frags = pepfrag_fragments(
        sequence=precursor.sequence,
        residue_range=residue_range,
        charge=peptide["charge"],
        bond_cys_res=bond_cys_res,
        ion_types=ion_types,
        count=count,
    )
    return simple_frags_to_fragments(frags, sequence, residue_range, precursor)


def fragment_dipeptide(
    peptide: Dict, bond_cys_res: List[int], ion_types=None, count=50
):
    if ion_types is None:
        ion_types = {IonType.y: [], IonType.b: [], IonType.precursor: []}

    max_charge = peptide["charge"]
    precursor: Precursor = peptide["precursor"]
    ps, qs = precursor.sequence.split("+")
    prr, qrr = precursor.residue_ranges

    result = []
    building_fragments = []
    for sequence, residue_range in [(ps, prr), (qs, qrr)]:
        frags = pepfrag_fragments(
            sequence=sequence,
            residue_range=residue_range,
            charge=1,
            ion_types=ion_types,
            bond_cys_res=bond_cys_res,
            count=count,
        )

        simple_frags, cys_frags = split_on_simple_frags(
            residue_range, frags, bond_cys_res
        )
        result += simple_frags_to_fragments(
            simple_frags, sequence, residue_range, precursor
        )
        shuffle(cys_frags)
        building_fragments.append(
            [
                fr["fragment"]
                for fr in simple_frags_to_fragments(
                    cys_frags, sequence, residue_range, precursor
                )
            ]
        )

    for i, (pf, qf) in enumerate(
        choices(list(zip(building_fragments[0], building_fragments[1])), k=count)
    ):
        total_charge = randint(1, max_charge)
        total_mass = pf.mz + qf.mz + CYS_BOND.mass - 2 * PROTON

        if "C" not in pf.sequence or "C" not in qf.sequence:
            continue

        fragment = Fragment(
            id=0,
            sequence=pf.sequence + "+" + qf.sequence,
            residue_ranges=pf.residue_ranges + qf.residue_ranges,
            intensity=10,
            intensity_ratio=1,
            mass=total_mass,
            target_mass=total_mass,
            mz=total_mass / total_charge + PROTON,
            target_mz=total_mass / total_charge + PROTON,
            charge=total_charge,
            break_count=pf.break_count + qf.break_count,
            error_ppm=0,
            modifications=qf.modifications + pf.modifications,
            connected_bonds=tuple([(72, 119)]),
            disconnected_cys=tuple([]),
        )

        result.append(
            {"fragment": fragment, "precursor": precursor, "var_bonds": [(72, 119)]}
        )
    return result


def fragment_unipeptide(
    peptide: Dict, bond_cys_res: List[int], ion_types=None, count=50
):
    if ion_types is None:
        ion_types = {IonType.y: [], IonType.b: [], IonType.precursor: []}

    max_charge = peptide["charge"]
    precursor: Precursor = peptide["precursor"]
    sequence = precursor.sequence
    residue_range = precursor.residue_ranges[0]

    frags = pepfrag_fragments(
        sequence=sequence,
        residue_range=residue_range,
        charge=1,
        ion_types=ion_types,
        bond_cys_res=bond_cys_res,
        count=count,
    )

    simple_frags, cys_frags = split_on_simple_frags(residue_range, frags, bond_cys_res)

    result = []

    b_ions, y_ions = [], []
    for frag in cys_frags:
        if "b" in frag[1]:
            b_ions.append(frag)
        else:
            y_ions.append(frag)

    fragments = []
    for ions in [b_ions, y_ions]:
        fragments.append(
            [
                fr["fragment"]
                for fr in simple_frags_to_fragments(
                    b_ions, sequence, residue_range, precursor
                )
            ]
        )

    for i, (pf, qf) in enumerate(
        choices(list(zip(fragments[0], fragments[1])), k=count)
    ):
        if "C" not in pf.sequence or "C" not in qf.sequence:
            continue
        total_charge = randint(1, max_charge)

        pr, qr = pf.residue_ranges[0], qf.residue_ranges[0]
        if intersects(pr, qr):
            continue

        total_mass = pf.mz + qf.mz + CYS_BOND.mass - 2 * PROTON

        fragment = Fragment(
            id=i,
            sequence=pf.sequence + "+" + qf.sequence,
            residue_ranges=pf.residue_ranges + qf.residue_ranges,
            intensity=10,
            intensity_ratio=1,
            mass=total_mass,
            target_mass=total_mass,
            mz=total_mass / total_charge + PROTON,
            target_mz=total_mass / total_charge + PROTON,
            charge=total_charge,
            break_count=2 if pr[1] != qr[0] else 1,
            error_ppm=0,
            modifications=qf.modifications + pf.modifications,
            connected_bonds=tuple([(72, 119)]),
            disconnected_cys=tuple([]),
        )

        result.append(
            {"fragment": fragment, "precursor": precursor, "var_bonds": [(72, 119)]}
        )
    return result


def generate_fragments(peptide: Dict, **kwargs):
    precursor: Precursor = peptide["precursor"]
    if precursor.cys_bond_count == 0:
        return fragment_simple_peptide(peptide, **kwargs)
    elif len(precursor.segments) == 2:
        return fragment_dipeptide(peptide, **kwargs)
    else:
        return fragment_unipeptide(peptide, **kwargs)


if __name__ == "__main__":
    import argparse

    args = argparse.ArgumentParser(description="Generate precursors and fragments")

    args.add_argument(
        "--protein",
        type=str,
        required=True,
        help="protein code (usually three letters)",
    )
    args.add_argument(
        "--kind",
        type=str,
        choices=["AT", "RAT"],
        required=True,
        help="measurement type (AT/RAT)",
    )
    args.add_argument(
        "--prec_error",
        type=int,
        required=True,
        help="allowed precursor error in ppm",
    )
    args.add_argument(
        "--frag_error",
        type=int,
        required=True,
        help="allowed fragment error in ppm",
    )
    args.add_argument(
        "--prec_segments",
        type=int,
        required=True,
        help="upper bound of segment count in matched precursors",
    )
    args.add_argument(
        "--frag_breaks",
        type=int,
        required=True,
        help="upper bound of break count in matched fragments",
    )
    args = args.parse_args()

    tryptides = cleave_protein(args.protein)

    print(f"Generating precursors...")
    precursors, cys_peptides = generate_simple_peptides(tryptides, [7, 10])
    if args.kind == "AT":
        precursors += generate_dipeptides(cys_peptides)
        precursors += generate_unipeptides(cys_peptides)

    print(f"In total there's {len(precursors)} precursors.")

    sys.exit()

    scans: List[Scan] = []
    fragment_records = []
    precursor_records = []

    for i, peptide in enumerate(precursors):
        precursor: Precursor = peptide["precursor"]
        fragments = generate_fragments(peptide, bond_cys_res=[72, 119])
        fragment_objects: List[Fragment] = [f["fragment"] for f in fragments]
        scan = Scan(
            nth_in_order=i,
            id=i,
            time=i,
            charge=peptide["charge"],
            prec_mz=precursor.mz,
            prec_intensity=100,
            prec_mass=precursor.mass,
            fragments_mz=np.array(sorted([f.mz for f in fragment_objects])),
            fragments_intensity=np.array([f.intensity for f in fragment_objects]),
            threshold=0,
        )
        scans.append(scan)

        precursor_records.append(scan.to_dict() | precursor.to_dict())

        fragment_records += [
            scan.to_dict()
            | fr["precursor"].to_dict()
            | {"var_bonds": fr["var_bonds"]}
            | fr["fragment"].to_dict()
            for fr in fragments
        ]

    precursor_path = (
        "../out/precursor_matches/{}_{}_segments={}_error={}ppm.pickle".format(
            args.protein, args.kind, args.prec_segments, args.prec_error
        )
    )
    precursor_matches = write_matched_precursors(
        tryptides,
        scans,
        precursor_path,
        max_segments=args.prec_segments,
        error_ppm=args.prec_error,
    )

    precursor_match_records = []
    for pm in precursor_matches:
        precursor_match_records.append(pm["scan"].to_dict() | pm["precursor"].to_dict())
    prec_df = pd.DataFrame(precursor_match_records)
    precursor_csv_path = (
        "../out/csv/precursor_matches_{}_{}_segments={}_error={}ppm.pickle".format(
            args.protein, args.kind, args.prec_segments, args.prec_error
        )
    )
    print(f"Saving precursor csv to {precursor_csv_path}")
    prec_df.to_csv(precursor_csv_path, index=False)

    fragment_path = (
        "../out/fragment_matches/{}_{}_segments={}_breaks={}_error={}ppm.pickle".format(
            args.protein,
            args.kind,
            args.prec_segments,
            args.frag_breaks,
            args.frag_error,
        )
    )
    print(f"Computing fragments...")
    fragment_matches = write_matched_fragments(
        precursor_matches=precursor_matches,
        tryptides=tryptides,
        output_path=fragment_path,
        max_allowed_breaks=args.frag_breaks,
        error_ppm=args.frag_error,
    )
    fragment_match_records = []
    for fm in fragment_matches:
        fragment_match_records.append(
            fm["scan"].to_dict()
            | fm["precursor"].to_dict()
            | fm["variant"].to_dict()
            | (fm["fragment"].to_dict() if fm["fragment"] is not None else {})
            | {"prec_variant_count": fm["variant_count"]}
        )

    frag_df = pd.DataFrame(fragment_match_records)
    fragment_csv_path = "../out/fragment_matches/fragment_matches_{}_{}_segments={}_breaks={}_error={}ppm.pickle".format(
        args.protein,
        args.kind,
        args.prec_segments,
        args.frag_breaks,
        args.frag_error,
    )
    print(f"Saving fragments csv to {fragment_csv_path}")
    frag_df.to_csv(fragment_csv_path, index=False)
