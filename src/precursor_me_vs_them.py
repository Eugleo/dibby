# precursor_matches_df = pd.DataFrame(precursor_matches)
# precursor_matches_df["scan_nth"] = [
#     scan.nth_in_order for scan in precursor_matches_df["scan"]
# ]
# precursor_matches_df = precursor_matches_df[
#     [
#         "scan_id",
#         "scan_nth",
#         "prec_mz",
#         "prec_sequence",
#         "prec_mass",
#         "prec_error",
#         "prec_mods",
#         "prec_cys_bond_count",
#         "prec_mc",
#     ]
# ]
#
# precursor_matches_df
#
# #%%
#
# from pymzid.read_mzid import Mzid
# from fragments import compute_error
#
# mgf_id = Mzid("../data/mgf/190318_LYS_AT_50x_05.mzid")
# mgf_id.read_psm()
# msgf_matches_df = mgf_id.psm_df
#
# msgf_matches_df.head(30)
#
# #%%
#
# msgf_matches_df["spectrum_id"] = [
#     int(s.removeprefix("index=")) for s in msgf_matches_df["spectrum_id"]
# ]
# msgf_matches_df["pep_id"] = [s.removeprefix("Pep_") for s in msgf_matches_df["pep_id"]]
# msgf_matches_df["their_error_ppm"] = [
#     compute_error(float(reference), float(measured))
#     for reference, measured in zip(msgf_matches_df["calc_mz"], msgf_matches_df["mz"])
# ]
#
# msgf_matches_df = msgf_matches_df[
#     ["spectrum_id", "calc_mz", "pep_id", "their_error_ppm"]
# ]
# msgf_matches_df = msgf_matches_df.rename(
#     columns={"pep_id": "their_sequence", "calc_mz": "their_mz"}
# )
# msgf_matches_df = msgf_matches_df.set_index("spectrum_id")
# msgf_matches_df
#
# #%%
#
# precursor_matches_df.join(msgf_matches_df, how="outer").to_csv(
#     "../out/my_vs_their_matches_rat.csv"
# )
