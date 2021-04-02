from protein import Protein, trypsin
from measurement import read_mgf
import pandas as pd

OVA = "GSIGAASMEFCFDVFKELKVHHANENIFYCPIAIMSALAMVYLGAKDSTRTQINKVVRFDKLPGFGDSIEAQCGTSVNVHSSLRDILNQITKPNDVYSFSLASRLYAEERYPILPEYLQCVKELYRGGLEPINFQTAADQARELINSWVESQTNGIIRNVLQPSSVDSQTAMVLVNAIVFKGLWEKAFKDEDTQAMPFRVTEQESKPVQMMYQIGLFRVASMASEKMKILELPFASGTMSMLVLLPDEVSGLEQLESIINFEKLTEWTSSNVMEERKIKVYLPRMKMEEKYNLTSVLMAMGITDVFSSSANLSGISSAESLKISQAVHAAHAEINEAGREVVGSAEAGVDAASVSEEFRADHPFLFCIKHIATNAVLFFGRCVSP"
LYS = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
BSA = "DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"

# LYS můstky
# VFGRCELAAA + WIRGCRL
# GNWVCAAKFE + WRNRCKGTDV
# SRWWCNDGRT + CNIPCSALLS
# SRNLCNIPCS + ASVNCAKKIV


protein = Protein(LYS)
peptides = protein.digest(trypsin)

# Protein Pilot
# Paragon
# AA tagy ze spektra

soft_err_ppm = 10
hard_err_ppm = 20
result = []
peps_with_threshold = [(pep, (hard_err_ppm / 1e6) * pep.total_mass) for pep in peptides]
for i, m in enumerate(read_mgf("../data/mgf/190318_LYS_AT_50x_05.mgf")):
    if i % 500 == 0:
        print(f"Done: {i}")

    for pep, t in peps_with_threshold:
        if abs(pep.total_mass - m.peptide_mz) <= t:
            score = m.score_match(
                pep,
                soft_err_ppm=soft_err_ppm,
                hard_err_ppm=hard_err_ppm,
            )

            result.append(
                {
                    "peptide": str(pep),
                    "is_valid": pep.ismerged or len(pep.peptides[0].seq) > 5,
                    "score": score,
                    "mod": pep.modstr,
                    "scan": m.scan,
                    "mass_error": pep.total_mass - m.peptide_mz,
                    "measurement": i,
                    # "silico_frags": ";".join(str(m) for m in pep.fragment_masses()),
                    # "frags": ";".join(str(f) for f in m.fragments_mz),
                    # "frags_intensity": ";".join(str(f) for f in m.fragments_intensity)
                }
            )
df = pd.DataFrame(result)
df.to_csv("../out/lys_peptide_matches_ultra_new_score.csv", index=False)

# Nechat pokud
# má vysokou pnost
# nebo obsahuje cystein a je kolem něj díra, a jinak má dobré pokrytí
# nebo je krátký a obsahuje cystein, pak na pokrytí nezáleží

# a u těchto se podívat na to, co v tom spektru přebývá?
