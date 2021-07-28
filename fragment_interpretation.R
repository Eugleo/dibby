library(tidyverse)

to_pairs <- function(xs) {
  (xs %>% str_split(","))[[1]] %>% parse_double()
}

df <-
  read_delim(
    "out/csv/fragment_matches_OVA_AT_segments=3_breaks=2_error=15ppm.csv",
    delim = ";"
  )

OVA_BONDS <- c("(72, 119)")
LYS_BONDS <- c("(5, 126)", "(29, 114)", "(75, 93)", "(63, 79)")

good_vars <-
  df %>%
  group_by(prec_sequence, var_bonds) %>%
  separate_rows(var_bonds, sep="\\ + ") %>%
  summarise(
    is_good = case_when(
      all(is.na(var_bonds)) ~ NA,
      TRUE ~ all(var_bonds %in% OVA_BONDS)
    )
  )

df <- df %>% left_join(good_vars)

# Počet variant a tryptide probability
df %>%
  group_by(scan_id) %>%
  summarise(var_count = n_distinct(variant_id), var_prob = first(variant_probability)) %>%
  ungroup() %>%
  filter(var_count < 150) %>%
  ggplot(aes(var_count, var_prob)) +
  geom_boxplot(aes(group=var_count)) +
  scale_y_log10()

# Počet segmentů a tryptide probability
df %>%
  group_by(scan, variant_id) %>%
  summarise(seg_count = lengths(str_split(variant_seq, "\\+")), var_prob = first(variant_probability)) %>%
  ungroup() %>%
  ggplot(aes(seg_count, var_prob)) +
  geom_boxplot(aes(group=seg_count)) +
  scale_y_log10()

# Počet segmentů
df %>%
  group_by(scan, variant_id) %>%
  summarise(seg_count = lengths(str_split(variant_seq, "\\+")), var_prob = first(variant_probability)) %>%
  ungroup() %>%
  ggplot(aes(seg_count)) +
  geom_bar()

# Počet namatchovaných fargmentů a probability
df %>%
  group_by(scan, variant_id) %>%
  summarise(match_count = n_distinct(target_mz), var_prob = first(variant_probability)) %>%
  ungroup() %>%
  ggplot(aes(match_count, var_prob)) +
  geom_point() +
  scale_y_log10()

# % namatchované intenzity a probability
df %>%
  group_by(scan, variant_id) %>%
  distinct(target_mz, .keep_all = TRUE) %>%
  summarise(
    match_percent = sum(intensity) / first(total_intensity),
    var_prob = first(variant_probability)
  ) %>%
  ungroup() %>%
  ggplot(aes(match_percent, var_prob)) +
  geom_point() +
  scale_y_log10()

# Jak vypadá pnost obsahující určité můstky
df %>%
  group_by(scan, variant_seq, variant_bonds) %>%
  distinct(target_mz, .keep_all = TRUE) %>%
  summarise(
    match_percent = sum(intensity) / first(total_intensity),
    var_prob = first(variant_probability)
  ) %>%
  ungroup() %>%
  ggplot(aes(match_percent, var_prob)) +
  geom_point()+
  scale_y_log10()

# Kolik je prekurzorových matchů s alespoň jednou správnou variantou?
df %>%
  filter(!is.na(is_good)) %>%
  group_by(scan, is_good) %>%
  summarise(n_variants = n_distinct(variant_id)) %>%
  summarise(
    ratio = sum(n_variants[is_good]) / sum(n_variants),
    absolute = sum(n_variants[is_good])
  ) %>%
  ggplot(aes(ratio > 0)) +
  geom_bar()

bad_precursor_scans <-
  df %>%
  filter(!is.na(is_good)) %>%
  group_by(scan, is_good) %>%
  summarise(n_variants = n_distinct(variant_id)) %>%
  summarise(
    ratio = sum(n_variants[is_good]) / sum(n_variants),
    absolute = sum(n_variants[is_good])
  ) %>%
  filter(absolute == 0) %>%
  select(scan)


df %>%
  semi_join(bad_precursor_scans) %>%
 filter(is.na(is_good))


df %>%
  separate_rows(connected_bonds, sep="\\+") %>%
  filter(!is.na(connected_bonds)) %>%
  ggplot(aes(connected_bonds)) +
  geom_bar() +
  theme(legend.position="none") +
  coord_flip()

df %>%
  separate_rows(precursor_bonds, sep="\\+") %>%
  filter(!is.na(precursor_bonds)) %>%
  ggplot(aes(precursor_bonds)) +
  geom_bar() +
  theme(legend.position="none") +
  coord_flip()


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
