import itertools

import math

from protein import Protein, trypsin
from measurement import read_mgf
import numpy as np
import pandas as pd

OVA = "GSIGAASMEFCFDVFKELKVHHANENIFYCPIAIMSALAMVYLGAKDSTRTQINKVVRFDKLPGFGDSIEAQCGTSVNVHSSLRDILNQITKPNDVYSFSLASRLYAEERYPILPEYLQCVKELYRGGLEPINFQTAADQARELINSWVESQTNGIIRNVLQPSSVDSQTAMVLVNAIVFKGLWEKAFKDEDTQAMPFRVTEQESKPVQMMYQIGLFRVASMASEKMKILELPFASGTMSMLVLLPDEVSGLEQLESIINFEKLTEWTSSNVMEERKIKVYLPRMKMEEKYNLTSVLMAMGITDVFSSSANLSGISSAESLKISQAVHAAHAEINEAGREVVGSAEAGVDAASVSEEFRADHPFLFCIKHIATNAVLFFGRCVSP"
LYS = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
BSA = "DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"

protein = Protein(LYS)
peptides = protein.peptides(trypsin)
mzs = {pep: pep.mz for pep in peptides}
has_cysteine = {pep: "C" in pep for pep in peptides}

# Protein Pilot
# Paragon
# AA tagy ze spektra

err_ppm = 10
result = []
for i, m in enumerate(read_mgf("../data/mgf/190318_LYS_AT_50x_05.mgf")):
    if i % 500 == 0:
        print(f"Done: {i}")

    for pep in peptides:
        if abs(mzs[pep] - m.peptide_mz) > err_ppm * mzs[pep]:
            score = m.score_match(pep, err_ppm=err_ppm, skip_cysteine=False)
        elif m.peptide_mz > mzs[pep] and has_cysteine[pep]:
            score = m.score_match(pep, err_ppm=err_ppm, skip_cysteine=True)
        else:
            score = 0

        result.append(
            {
                "peptide": pep.seq,
                "score": score,
                "mod": pep.modstr,
                "scan": m.scan,
                "measurement": i,
                # "silico_frags": ";".join(str(m) for m in pep.fragment_masses()),
                # "frags": ";".join(str(f) for f in m.fragments_mz),
                # "frags_intensity": ";".join(str(f) for f in m.fragments_intensity)
            }
        )

#df = pd.DataFrame(result)
#df.to_csv("../out/lys_peptide_matches_new_score.csv", index=False)

# Nechat pokud
# má vysokou pnost
# nebo obsahuje cystein a je kolem něj díra, a jinak má dobré pokrytí
# nebo je krátký a obsahuje cystein, pak na pokrytí nezáleží

# a u těchto se podívat na to, co v tom spektru přebývá?

