library(tidyverse)

df <-
  read_delim(
    "out/csv/fragment_matches_LIP_AT_segments=3_breaks=2_error=15ppm.csv",
    delim = ";"
  )

df %>% ggplot(aes(var_is_good)) + geom_bar() + scale_y_log10()

OVA_BONDS <- c("(72, 119)")
LYS_BONDS <- c("(5, 126)", "(29, 114)", "(75, 93)", "(63, 79)")
BSA_BONDS <-
  c('(52, 61)',
    '(74, 90)',
    '(89, 100)',
    '(122, 167)',
    '(166, 175)',
    '(198, 244)',
    '(243, 251)',
    '(263, 277)',
    '(276, 287)',
    '(314, 359)',
    '(358, 367)',
    '(390, 436)',
    '(435, 446)',
    '(459, 475)',
    '(474, 485)',
    '(512, 557)',
    '(556, 565)')
LYS_BONDS_NAMES <- c("var_has_5_126", "var_has_29_114", "var_has_75_93", "var_has_63_79")

prec_is_good <-
  df %>%
  filter(!is.na(var_is_good)) %>%
  group_by(scan_id, prec_sequence) %>%
  summarise(prec_is_good = sum(var_is_good) > 0)
df <- df %>% left_join(prec_is_good)

scan_is_good <-
  df %>%
  filter(!is.na(prec_is_good)) %>%
  group_by(scan_id) %>%
  summarise(scan_is_good = sum(prec_is_good) > 0)
df <- df %>% left_join(scan_is_good)

# Dobrých scanů je zhruba čtvrtina oproti těm špatným
df %>%
  distinct(scan_id, scan_is_good) %>%
  ggplot(aes(scan_is_good)) +
  geom_bar()

# Dobrých prekurzorů už je méně
df %>%
  distinct(scan_id, prec_sequence, prec_is_good) %>%
  ggplot(aes(prec_is_good)) +
  geom_bar()

# ...a dobrých variant ještě méně
df %>%
  distinct(scan_id, prec_sequence, var_bonds, var_is_good) %>%
  ggplot(aes(var_is_good)) +
  geom_bar()

plot_diff <- function(df, ...) {
  df %>%
    pivot_longer(
      c(...),
      names_to = "variable",
      values_to = "value"
    ) %>%
    ggplot(aes(var_is_good, value)) +
    geom_boxplot(aes(color=prec_is_good, fill=scan_is_good)) +
    facet_wrap(vars(variable), ncol = 4, scales = "free_y") +
    scale_color_manual(values=c("#215264", "#26734A")) +
    scale_fill_manual(values = c("#C16E45", "#D3B27C"))
}

var_intensity_ratios <-
  df %>%
  group_by(scan_id, prec_sequence, var_bonds) %>%
  distinct(frag_id, .keep_all = TRUE) %>%
  summarise(var_intensity_ratio =
              sum(frag_intensity) / first(scan_total_intensity)
  ) %>%
  ungroup()

df <-
  df %>%
  left_join(var_intensity_ratios) %>%
  mutate(
    prec_mod_count = lengths(str_split(prec_mods, " \\+ ")),
    frag_mod_count = lengths(str_split(frag_mods, " \\+ ")),
    frag_disconnected_cys_count =
      str_split(frag_disconnected_cys, " \\+ ") %>%
      map_dbl(~sum(!is.na(.x))),
    var_has_bonds =
      str_split(var_bonds, " \\+ ") %>%
      map(~ as_tibble_row(setNames(LYS_BONDS %in% .x, LYS_BONDS_NAMES)))
  ) %>%
  unnest(var_has_bonds) %>%
  group_by(scan_id, prec_sequence, var_bonds) %>%
  mutate(
    across(
      c(frag_mass, frag_charge, frag_error_ppm, frag_mod_count, frag_break_count, frag_intensity_ratio, frag_disconnected_cys_count),
      list(median = median)
    )
  )

df %>%
  plot_diff(
    frag_mass, frag_mass_median,
    frag_charge, frag_charge_median,
    frag_error_ppm, frag_error_ppm_median,
    frag_mod_count, frag_mod_count_median,
    frag_break_count, frag_break_count_median,
    frag_intensity_ratio, frag_intensity_ratio_median,
    frag_disconnected_cys_count, frag_disconnected_cys_count_median
  )

df %>%
  plot_diff(
    prec_charge, prec_segment_count, prec_max_mc_count, prec_cys_bond_count,
    prec_mass, prec_error, prec_alkylation_count, prec_variant_count,
    prec_mod_count, var_intensity_ratio
  )

df %>%
  group_by(scan_id, target_mz) %>%
  summarise(n = n_distinct(frag_sequence)) %>%
  filter(n < 25) %>%
  ggplot(aes(n)) +
  geom_bar()

df %>%
  pivot_longer(c(var_has_5_126, var_has_29_114, var_has_63_79, var_has_75_93)) %>%
  ggplot(aes(var_is_good)) +
  geom_bar() +
  facet_wrap(~ name)

# Chybí nám data o posledních dvou můstcích
df %>% filter(var_has_5_126) %>% ggplot(aes(var_is_good)) + geom_bar()
df %>% filter(var_has_29_114) %>% ggplot(aes(var_is_good)) + geom_bar()
df %>% filter(var_has_63_79) %>% ggplot(aes(var_is_good)) + geom_bar()
df %>% filter(var_has_75_93) %>% ggplot(aes(var_is_good)) + geom_bar()

df %>%
  distinct(scan_id, prec_sequence, var_bonds, .keep_all = TRUE) %>%
  filter(var_is_good & prec_is_good)

# df %>%
#   separate_rows(connected_bonds, sep="\\+") %>%
#   filter(!is.na(connected_bonds)) %>%
#   ggplot(aes(connected_bonds)) +
#   geom_bar() +
#   theme(legend.position="none") +
#   coord_flip()
#
# df %>%
#   separate_rows(precursor_bonds, sep="\\+") %>%
#   filter(!is.na(precursor_bonds)) %>%
#   ggplot(aes(precursor_bonds)) +
#   geom_bar() +
#   theme(legend.position="none") +
#   coord_flip()

# Old analysis ------------------------------------------------------------

# Počet variant a tryptide probability
# df %>%
#   group_by(scan_id) %>%
#   summarise(var_count = n_distinct(variant_id), var_prob = first(variant_probability)) %>%
#   ungroup() %>%
#   filter(var_count < 150) %>%
#   ggplot(aes(var_count, var_prob)) +
#   geom_boxplot(aes(group=var_count)) +
#   scale_y_log10()

# Počet segmentů a tryptide probability
# df %>%
#   group_by(scan, variant_id) %>%
#   summarise(seg_count = lengths(str_split(variant_seq, "\\+")), var_prob = first(variant_probability)) %>%
#   ungroup() %>%
#   ggplot(aes(seg_count, var_prob)) +
#   geom_boxplot(aes(group=seg_count)) +
#   scale_y_log10()

# Počet segmentů
# df %>%
#   group_by(scan, variant_id) %>%
#   summarise(seg_count = lengths(str_split(variant_seq, "\\+")), var_prob = first(variant_probability)) %>%
#   ungroup() %>%
#   ggplot(aes(seg_count)) +
#   geom_bar()

# Počet namatchovaných fargmentů a probability
# df %>%
#   group_by(scan, variant_id) %>%
#   summarise(match_count = n_distinct(target_mz), var_prob = first(variant_probability)) %>%
#   ungroup() %>%
#   ggplot(aes(match_count, var_prob)) +
#   geom_point() +
#   scale_y_log10()

# % namatchované intenzity a probability
# df %>%
#   group_by(scan, variant_id) %>%
#   distinct(target_mz, .keep_all = TRUE) %>%
#   summarise(
#     match_percent = sum(intensity) / first(total_intensity),
#     var_prob = first(variant_probability)
#   ) %>%
#   ungroup() %>%
#   ggplot(aes(match_percent, var_prob)) +
#   geom_point() +
#   scale_y_log10()

# Jak vypadá pnost obsahující určité můstky
# df %>%
#   group_by(scan_id, prec_sequence, var_bonds) %>%
#   distinct(target_mz, .keep_all = TRUE) %>%
#   summarise(
#     match_percent = sum(frag_intensity) / first(scan_total_intensity),
#     var_prob = first(variant_probability)
#   ) %>%
#   ungroup() %>%
#   ggplot(aes(match_percent, var_prob)) +
#   geom_point()+
#   scale_y_log10()

# Kolik prekurzorů se většinou matchne na jeden scan
df %>%
  filter(scan > 4000) %>%
  group_by(scan) %>%
  summarise(n = n_distinct(precursor_seq)) %>%
  ungroup() %>%
  ggplot(aes(scan, n, color = n > 1)) +
  geom_col()

# Kolik prekurzorů se většinou matchne na jeden scan (boxplot)
df %>%
  filter(scan > 4000) %>%
  group_by(scan) %>%
  summarise(n = n_distinct(precursor_seq)) %>%
  ungroup() %>%
  mutate(gr = "T") %>%
  ggplot(aes(gr, n)) +
  geom_boxplot()



# SD of scores of different matches in scans
df %>%
  filter(scan > 2000) %>%
  group_by(scan, precursor_seq, precursor_bonds) %>%
  summarise(score = sum(intensity) / first(total_intensity)) %>%
  count(scan, score) %>%
  ungroup() %>%
  group_by(scan) %>%
  mutate(deviation = sd(score)) %>%
  filter(n > 1) %>%
  ggplot(aes(scan, deviation)) +
  geom_col(aes(color = scan)) +
  theme(legend.position="none")


# Matching precursor vs MSGF+ ---------------------------------------------

# Na LYS AT

my_vs_their <- read_csv("out/my_vs_their_matches.csv")

# Které našli oni a my ne
# Nechybí nám žádný prekurzor (jen ty s chybou > 200ppm)
my_vs_their %>%
  filter(
    is.na(my_sequence) | my_sequence != their_sequence
  ) %>%
  arrange(spectrum_id) %>%
  View()

# Které oni našli a my ne, když nám vk jednomu spektru stačí jeden společný match
my_vs_their %>%
  mutate(
    their_sequence_nomod = str_remove_all(their_sequence, "\\+\\d+")
  ) %>%
  group_by(spectrum_id) %>%
  mutate(
    no_match =
      all(is.na(my_sequence)) ||
      all(their_sequence_nomod != my_sequence)
  ) %>%
  ungroup() %>%
  filter(no_match) %>%
  View()

# Které oni našli a my ne, když porovnáme můj a jejich top match v daném scanu
my_vs_their %>%
  mutate(
    their_sequence_nomod = str_remove_all(their_sequence, "\\+\\d+")
  ) %>%
  group_by(spectrum_id) %>%
  arrange(my_error_ppm) %>%
  mutate(
    no_match =
      all(is.na(my_sequence)) ||
      first(their_sequence_nomod) != first(my_sequence)
  ) %>%
  ungroup() %>%
  filter(no_match) %>%
  View()


# Které jsme my našli a oni ne
# Nadbývá nám několik prekurzorů
my_vs_their %>%
  mutate(
    their_sequence_nomod = str_remove_all(their_sequence, "\\+\\d+")
  ) %>%
  group_by(spectrum_id) %>%
  arrange(my_error_ppm) %>%
  mutate(
    no_match =
      all(is.na(their_sequence)) ||
      all(their_sequence_nomod != my_sequence)
  ) %>%
  ungroup() %>%
  filter(no_match & cys_bonds == 0) %>%
  View()
