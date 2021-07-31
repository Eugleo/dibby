library(tidyverse)

ideal <-
  read_delim(
    "out/csv/precursor_matches_TEST_OVA_segments=3_error=50ppm_ideal.csv",
    delim = ","
  ) %>%
  select(-scan_nth_in_order, -scan_time, -scan_total_intensity)

df <- read_delim(
  "out/csv/precursor_matches_TEST_OVA_segments=3_error=50ppm.csv",
  delim = ","
) %>%
  select(-scan_nth_in_order, -scan_time, -scan_total_intensity)

identical <- function(xs, ys) {
  !(is.na(xs) | is.na(ys)) & xs == ys
}

prec <-
  df %>%
  full_join(ideal, by = c("scan_id"), suffix = c("", "_ideal")) %>%
  mutate(
    prec_identical =
      identical(prec_sequence, prec_sequence_ideal) &
        identical(prec_segment_count, prec_segment_count_ideal) &
        identical(prec_tryptide_ranges, prec_tryptide_ranges_ideal) &
        identical(prec_residue_ranges, prec_residue_ranges_ideal) &
        identical(prec_max_mc_count, prec_max_mc_count_ideal) &
        identical(prec_mc, prec_mc) &
        identical(prec_cys_bond_count, prec_cys_bond_count_ideal) &
        identical(prec_alkylation_count, prec_alkylation_count_ideal) &
        identical(prec_mods, prec_mods_ideal),
    prec_big_error = !identical(prec_sequence, prec_sequence_ideal)
  )

# Běžná chyba špatných věcí
df %>%
  ggplot(aes(prec_error)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~identical, scales = "free_y")

# False positives
prec %>%
  filter(
    prec_big_error,
    is.na(prec_error) | prec_error < 5,
  ) %>%
  View()

# 858 true scanů
# 1,283 matchů, z toho
# - 858 TP
# - 425 FP, a z nich
#   - 365 má MC < 5, a znich
#   - 260 má chybu < 5ppm, a z nich
#   - 20 nejsou jen shifty správného prekurzoru, a z nich
#   - zbydou 4 unikátní sekvence, když pominu jejich náboj

# FP když odstraním věci s posunem breaku uprostřed
prec %>%
  mutate(
    prec_break_removed = str_remove(prec_sequence, "\\+"),
  ) %>%
  filter(
    !prec_identical,
    prec_max_mc_count < 5,
    is.na(prec_error) | prec_error < 5,
    # str_remove(prec_sequence_ideal, "\\+") != prec_break_removed,
  )
  # distinct(prec_break_removed, prec_sequence_ideal, .keep_all = TRUE)

# FN
df %>%
  filter(big_error, is.na(prec_error))

# Scany s 0 TP
df %>%
  filter(prec_cys_bond_count_ideal == 0) %>%
  group_by(scan_id) %>%
  mutate(no_match = sum(identical) == 0) %>%
  filter(no_match)


# Fragments ---------------------------------------------------------------

ideal <-
  read_delim(
    "out/csv/fragment_matches_TEST_OVA_segments=3_error=50ppm_ideal.csv",
    delim = ","
  ) %>%
  select(
    -scan_nth_in_order, -scan_time, -scan_total_intensity,
    -frag_intensity, -frag_intensity_ratio, -target_mass,
  )

df <- read_delim(
  "out/csv/fragment_matches_TEST_OVA_segments=3_error=50ppm.csv",
  delim = ","
) %>%
  select(
    -scan_nth_in_order, -scan_time, -scan_total_intensity,
    -frag_intensity, -frag_intensity_ratio, -target_mass,
  )
frags <-
  df %>%
  full_join(ideal, by = c("scan_id", "target_mz"), suffix = c("", "_ideal")) %>%
  mutate(
    prec_identical =
      identical(prec_sequence, prec_sequence_ideal) &
        identical(prec_segment_count, prec_segment_count_ideal) &
        identical(prec_tryptide_ranges, prec_tryptide_ranges_ideal) &
        identical(prec_residue_ranges, prec_residue_ranges_ideal) &
        identical(prec_max_mc_count, prec_max_mc_count_ideal) &
        identical(prec_mc, prec_mc) &
        identical(prec_cys_bond_count, prec_cys_bond_count_ideal) &
        identical(prec_alkylation_count, prec_alkylation_count_ideal) &
        identical(prec_mods, prec_mods_ideal),
    frag_identical =
      identical(frag_sequence, frag_sequence_ideal) &
        identical(frag_residue_ranges, frag_residue_ranges_ideal) &
        identical(frag_charge, frag_charge_ideal) &
        identical(frag_break_count, frag_break_count_ideal) &
        identical(frag_connected_bonds, frag_connected_bonds_ideal) &
        identical(frag_disconnected_cys, frag_disconnected_cys_ideal) &
        identical(frag_interesting_disconnected_cys, frag_interesting_disconnected_cys_ideal) &
        identical(prec_mods, prec_mods_ideal),
    frag_big_error = !identical(frag_sequence, frag_sequence_ideal)
  )

# Převážná většina má chybu > 5
# Zbyde cca 900, z nich polovina jsou sekvenční duplikáty
# Zbytečně tam hážu několik breaků — JAK Z TOHO?
frags %>%
  filter(
    frag_big_error,
    # prec_cys_bond_count_ideal == 0,
    # frag_charge <= 3,
    prec_error < 5,
    frag_error_ppm < 5,
    if_else(prec_cys_bond_count == 0, frag_break_count < 2, TRUE),
  ) %>%
  distinct(frag_sequence, frag_sequence_ideal, .keep_all = TRUE) %>%
  select(
    scan_id, prec_identical, prec_sequence, prec_sequence_ideal, frag_sequence,
    frag_sequence_ideal,
    frag_charge, frag_charge_ideal, target_mz, frag_mz, frag_mz_ideal,
    frag_error_ppm
  ) %>%
  View()

# 36,117 TP
# 21,660 FP, ale z nich
# - jen 18,365 s prec chybou < 5, a znich
# - jen 3,389 s frag chybou < 5, a z nich
# - jen 177 takových, který nejsou z bonded prec a mají < 2 breaky
# - ...a tam je 71 unikátních kombinací

# FN
frags %>%
  # filter(prec_cys_bond_count_ideal == 0) %>%
  filter(is.na(prec_sequence)) %>%
  View()
