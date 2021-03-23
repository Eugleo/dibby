#!/usr/bin/env python3

import re
from pyteomics import mgf, mass

from Bio.SeqUtils.ProtParam import ProteinAnalysis

__proton_weight__ = mass.calculate_mass(formula="H").conjugate()


class MZ:
    def __init__(self, scan, time, intensity, charge, mz, est_mass, pepmass):
        self.scan = scan
        self.time = time
        self.intensity = intensity
        self.charge = charge
        self.mz = mz
        self.mass = est_mass
        self.pepmass = pepmass

    def __str__(self):
        return (
            f"<MZ scan={self.scan} time={self.time} intensity={self.intensity} charge={self.charge}"
            + " mz={self.mz} mass_estimate={self.mass} pepmass={self.pepmass}>"
        )


def read_mgf(fn):
    """
    returns (scan ID, time, charge, mz, mass estimate)
    """
    with mgf.read(fn) as reader:
        for i in reader:
            # new version
            # scan=int(re.match(".* scans: \"([0-9]+)\"", i['params']['title'])[1])
            # old version
            scan = int(re.match(".* scan=([0-9]+)", i["params"]["title"])[1])

            time = i["params"]["rtinseconds"]
            chargelist = i["params"]["charge"]
            if len(chargelist) > 1:
                raise AssertionError("ChargeList length>1 unsupported")
            charge = int(chargelist[0])
            # print(i['params'].keys())
            mz = i["params"]["pepmass"][0]
            intensity = i["params"]["pepmass"][1]
            # print(mz)
            est_mass = mz * charge - charge * __proton_weight__
            yield MZ(scan, time, intensity, charge, mz, est_mass, mz)


def estimate_peptide(seq, monoiso):
    # return mass.calculate_mass(sequence=seq).conjugate()
    return ProteinAnalysis(seq, monoisotopic=monoiso).molecular_weight()


def estimate_unlink_peptides(seq1, seq2):
    return estimate_peptide(seq1) + estimate_peptide(seq2)


def sulf_link_mass_decrease():
    return mass.calculate_mass(formula="H2").conjugate()


def estimate_sulf_link_peptides(seq1, seq2, nlinks=1):
    return estimate_unlink_peptides(seq1, seq2) - nlinks * sulf_link_mass_decrease()


def find_matches(ms, mass, eps=0.02):
    for i in ms:
        if abs(i.mass - mass) < eps:
            yield (i.scan, i.mz, i.charge, i.mass, i.time, i.intensity, i.pepmass)


def break_peptide(pep, maxseq=3, mis=False):
    # Poiss rozdělení možná
    parts = []
    last = 0
    for i in range(len(pep) - 1):
        if pep[i] in ["K", "R"] and pep[i + 1] != "P":
            parts.append((last, i + 1))
            last = i + 1
    parts.append((last, len(pep)))

    res = set()
    for l in range(1, maxseq + 1):
        for b in range(len(parts) - l + 1):
            begin = parts[b][0]
            end = parts[b + l - 1][1]
            res.add((begin, end))
            if mis and end < len(pep):
                end += 1
                res.add((begin, end))

    return sorted([([(b, e)], []) for (b, e) in res])


def subsets(x):
    if not x:
        return [[]]
    else:
        xs = subsets(x[1:])
        return [[x[0]] + a for a in xs] + xs


def connection_valid(rs1, rs2):
    rs = sorted(rs1 + rs2)
    laste = -1
    for (b, e) in rs:
        if laste > b:
            return False
        laste = e
    return True


def is_in_ranges(pos, rs):
    return any(b <= pos and pos < e for (b, e) in rs)


def connect_sulfide_bridge(fs, lb, le):
    if lb >= le:
        raise ValueError("lb < le, please")
    lb -= 1
    le -= 1

    res = []
    for (rngs, sbs) in fs:
        if is_in_ranges(lb, rngs):  # we must connect this with others
            for (rngs2, sbs2) in fs:
                if not is_in_ranges(le, rngs2):
                    continue
                elif (rngs, sbs) == (rngs2, sbs2):
                    res.append((rngs, sorted(sbs + [(lb, le)])))
                elif connection_valid(rngs, rngs2):
                    res.append(
                        (
                            sorted(list(set(rngs + rngs2))),
                            sorted(sbs + sbs2 + [(lb, le)]),
                        )
                    )
                else:
                    pass  # the connection could not have happened
        elif is_in_ranges(le, rngs):
            continue  # ends are handled elsewhere
        else:  # nothing happened, just carry on
            res.append((rngs, sbs))

    return sorted(res)


def incrboths(x, fst=1, snd=1):
    return [(a + fst, b + snd) for (a, b) in x]


# Můstky někdy váží více / méně, třeba to mít označené
def frag_to_desc(pep, rngs, sbs, monoiso):
    subpeps = [pep[b:e] for (b, e) in rngs]
    mass = (
        sum(estimate_peptide(subpep, monoiso) for subpep in subpeps)
        - len(sbs) * sulf_link_mass_decrease()
    )
    return {
        "ranges": incrboths(rngs, 1, 0),
        "seqs": subpeps,
        "mass": mass,
        "sulf_bridges": incrboths(sbs),
    }


def connect_sulfide_bridges(fs, bs):
    for (b, e) in bs:
        fs = connect_sulfide_bridge(fs, b, e)
    return fs


def fragments_to_mols(pep, fs, monoiso):
    return [frag_to_desc(pep, rngs, sbs, monoiso) for (rngs, sbs) in fs]


def search_mols(ms, mols, eps):
    for m in mols:
        r = dict(m)
        r["matches"] = list(find_matches(ms, m["mass"], eps))
        r["match_count"] = len(r["matches"])
        # if(r['match_count']==1): print(r['matches'])
        yield r


def search_peptide_fragments(
    p, ms, bridges, maxseq=3, miscleavage=False, eps=0.1, monoiso=True
):
    frags = break_peptide(p, maxseq=maxseq, mis=miscleavage)
    # Spojení více můstky?
    frags = connect_sulfide_bridges(frags, bridges)

    return search_mols(ms, fragments_to_mols(p, frags, monoiso), eps)


def search_peptide_fragments_to_tsv(out, *args, **kwargs):
    d = search_peptide_fragments(*args, **kwargs)

    with open(out, "w") as of:
        of.write(
            "seq_pos\tbridges\test_mass\tseq\tmatch_count\tscan_id\tmass_comp\tcharge\tpepmass\tscan_time\tpeak_intensity\n"
        )
        for i in d:
            pfx = "\t".join(
                [
                    ",".join([str(b) + "--" + str(e) for (b, e) in i["ranges"]]),
                    ",".join([str(b) + "--" + str(e) for (b, e) in i["sulf_bridges"]]),
                    str(i["mass"]),
                    ",".join(i["seqs"]),
                    str(i["match_count"]),
                ]
            )
            for (sc, mz, z, mass, time, intensity, pepmass) in i["matches"]:
                of.write(
                    "\t".join(
                        [
                            pfx,
                            str(sc),
                            str(mass),
                            str(z),
                            str(pepmass),
                            str(time),
                            str(intensity),
                        ]
                    )
                    + "\n"
                )


def gen_bridges_between(aa, seq):
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] == aa and seq[j] == aa:
                yield (i + 1, j + 1)


def bridges_between(aa, seq):
    return list(gen_bridges_between(aa, seq))


def test_peptide():
    return "DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"


# eps může být menší
data = list(read_mgf("190318_LYS_AT_50x_05.mgf"))
search_peptide_fragments_to_tsv(
    out="AT_bridges.tsv",
    p=test_peptide(),
    ms=data,
    # bridges=bridges_between('C',test_peptide()), #all bridges (too slow)
    bridges=[(6, 127), (30, 115), (64, 80), (76, 94)],
    eps=0.02,
    monoiso=True,
)
search_peptide_fragments_to_tsv(
    out="AT_nobridge.tsv", p=test_peptide(), ms=data, bridges=[], eps=0.02, monoiso=True
)

data = list(read_mgf("190318_LYS_RAT_50x_05.mgf"))
search_peptide_fragments_to_tsv(
    out="RAT_bridges.tsv",
    p=test_peptide(),
    ms=data,
    bridges=[(6, 127), (30, 115), (64, 80), (76, 94)],
    eps=0.02,
    monoiso=True,
)
search_peptide_fragments_to_tsv(
    out="RAT_nobridge.tsv",
    p=test_peptide(),
    ms=data,
    bridges=[],
    eps=0.02,
    monoiso=True,
)

# TODO
# projít scholar
# bibliografie
# k čemu to bude, pokud to uděláme

# NICE TO HAVE
# automatické srovnání controlu s něčím (Bayes + likelihood?)
# kde mohou být můstky?, pokud je nemáme dané
#   heuristiky

# BP
# které můstky tam jsou, které ne, a proč (report)
# check, jestli se peptid bude dělit na MSMS fragmenty (pnost?)
# proč je to složité (matika), nerozdhodnutelnost (příklady), obrázek

# SUPER NICE TO HAVE (practically required)
# jakákali heuristika (prožezávání atp)

# MAYBE
# genetický přístup?
#   pstor fragmentů, skóre = v měření se nám objevil náš fragment

# CO OD PETRY
# přesná chemie
# data označená + neoznačená
# k čemu to bude, až to uděláme
# jak poznáme, že tam předtím byl můstek
# mění se distribuce spekter když něco nalepíme na cysteniy
