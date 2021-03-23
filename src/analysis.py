import itertools

import math

from protein import Protein, trypsin
from measurement import read_mgf
import numpy as np
import pandas as pd

LYS = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
BSA = "DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"

protein = Protein(LYS)
peptides = list(protein.peptides(trypsin))
mzs = {pep: pep.mz for pep in peptides}

tolerance = 0.1
result = []
count = 0
for m in read_mgf("../data/mgf/190318_BSA_AT_50x_05.mgf"):
    count += 1
    if count % 100 == 0:
        print(f"Done: {count}")

    res = []
    for pep in peptides:
        if mzs[pep] > m.peptide_mz + tolerance:
            res.append(0)
        else:
            res.append(m.score_match(pep, tolerance=0.02))

    result.append(res)

df = pd.DataFrame(result, columns=[f"{n}" for n in range(np.shape(result)[1])])
df.to_csv("../out/lys_peptide_matches.csv")

# Nechat pokud
# má vysokou pnost
# nebo obsahuje cystein a je kolem něj díra, a jinak má dobré pokrytí
# nebo je krátký a obsahuje cystein, pak na pokrytí nezáleží

# a u těchto se podívat na to, co v tom spektru přebývá?
