library(tidyverse)
# install.packages("gridExtra")
library(gridExtra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

get_paths <- function(names, type, prefix = "190318", suffix = "50x_05") {
  fp <- file.path("..", "data", "csv")
  names %>%
    map(function(name) {
      filename <- paste0(paste(prefix, name, type, suffix, sep = "_"), ".csv")
      list(
        path = file.path(fp, filename),
        name = name,
        type = type
      )
    })
}

PROTS <- c("BSA", "LYS", "LIP", "OVA")

df <-
  c(get_paths(PROTS, "AT"), get_paths(PROTS, "RAT")) %>%
  map(~ read_csv(.x$path) %>% mutate(type = .x$type, protein = .x$name)) %>%
  bind_rows()

na_scans <-
  df_at %>%
  group_by(scan) %>%
  select_if(~ any(is.na(.x))) %>%
  summarise_all(~ sum(is.na(.))) %>%
  filter(peptide_intensity > 0)

cleanup_spectra <- function(df) {
  df %>% filter(
    (protein != "LYS" | total_intensity < 3e8) &
      (protein != "LIP" | total_intensity < 1e8)
  )
}

plots <-
  PROTS %>%
  map(
    function(p) {
      df %>%
        filter(protein == p) %>%
        group_by(protein, type, scan) %>%
        summarise(
          peptide_mz = first(peptide_mz),
          peptide_intensity = first(peptide_intensity)
        ) %>%
        group_by(protein, type, peptide_mz) %>%
        summarise(total_intensity = sum(peptide_intensity)) %>%
        cleanup_spectra() %>%
        ggplot(aes(peptide_mz, total_intensity)) +
        geom_col(width = 2) +
        facet_wrap(~type, ncol = 1, scales = "free_y") +
        labs(title = p)
    }
  )

plots %>% grid.arrange(grobs = .)

df_at %>%
  filter(scan == 2206) %>%
  ggplot(aes(fragment_mz, fragment_intensity)) +
  geom_col(width = 2)


# RAT AT analysis ---------------------------------------------------------

at_pepmass <-
  df %>%
  filter(protein == "BSA") %>%
  filter(type == "AT") %>%
  group_by(scan) %>%
  summarise(peptide_mz = first(peptide_mz)) %>%
  filter(!is.na(peptide_mz)) %>%
  pull(peptide_mz) %>%
  unique() %>%
  sort()


rat_pepmass <-
  df %>%
  filter(protein == "BSA") %>%
  filter(type == "RAT") %>%
  group_by(scan) %>%
  summarise(peptide_mz = first(peptide_mz)) %>%
  filter(!is.na(peptide_mz)) %>%
  pull(peptide_mz) %>%
  unique() %>%
  sort()

findInterval(at_pepmass, rat_pepmass)

count <- 0
i <- 1
j <- 1
imax <- length(at_pepmass)
while (i <= imax) {
  at <- at_pepmass[i]
  rat <- rat_pepmass[j]

  if (near(at, rat, tol = 0.01)) {
    i <- i + 1
  }
  else if (at < rat) {
    count <- count + 1
    i <- i + 1
  }
  else {
    j <- j + 1
  }
}


# Matches -----------------------------------------------------------------

sort_by <- function(xs, fun = str_length) {
  u <- xs %>% unique()
  l <- u %>%
    map_dbl(fun) %>%
    sort(index.return = TRUE)
  u[l$ix]
}

pepmod_lvls <- function(peptides, mods, pepmods) {
  transpose(list(peptides, mods, pepmods)) %>%
    sort_by(~ str_length(.x[[1]])) %>%
    map_chr(~ .x[[3]])
}

isok <- function(peptides) {
  ((peptides %>% map_dbl(str_length)) > 5) #| (peptides %>% map_lgl(~ str_detect(.x, "C")))
}

mch <-
  read_csv("../out/lys_analysis_rat.csv") %>%
  mutate(
    measurement = measurement + 1,
    pepmod = if_else(is.na(mod), peptide, paste(peptide, mod, sep = "_"))
  ) %>%
  mutate(
    pepmod = factor(pepmod, levels = pepmod_lvls(peptide, mod, pepmod))
  )

mch2 <-
  read_csv("../out/lys_analysis_at.csv") %>%
  mutate(
    measurement = measurement + 1,
    pepmod = if_else(is.na(mod), peptide, paste(peptide, mod, sep = "_")),
    pepid = as.numeric(as.factor(pepmod))
  ) %>%
  mutate(
    pepmod = factor(pepmod, levels = pepmod_lvls(peptide, mod, pepmod))
  )

# mch2 <-
#   read_csv("../out/lys_peptide_matches_ultra_new_score.csv") %>%
#   mutate(
#     pepmod = if_else(is.na(mod), peptide, paste(peptide, mod, sep = "_"))
#   ) %>%
#   mutate(
#     pepmod = factor(pepmod, levels = pepmod_lvls(peptide, mod, pepmod))
#   )
#
# mch2 %>%
#   filter(isok(peptide)) %>%
#   pivot_longer(c(score, new_score)) %>%
#   ggplot(aes(value)) +
#   geom_histogram(bins=60) +
#   scale_y_continuous(trans='log10') +
#   scale_x_continuous(limits=c(0, 1)) +
#   facet_wrap(~name, ncol = 1)

scores <-
  mch2 %>%
  filter(isok(peptide)) %>%
  group_by(measurement) %>%
  summarise(score = max(score))

scores %>%
  mutate(score_bin = cut_width(score, width = 0.1, center = 0.1)) %>%
  ggplot(aes(score_bin)) +
  geom_bar() +
  scale_y_continuous(trans = "log10")

# LYS můstky
# VFGRCELAAA + WIRGCRL
# GNWVCAAKFE + WRNRCKGTDV
# SRWWCNDGRT + CNIPCSALLS
# SRNLCNIPCS + ASVNCAKKIV

View(mch2 %>% filter(str_detect(peptide, "\\+")) %>% arrange(desc(score)))

View(
  mch2 %>%
    arrange(desc(score)) %>%
    select(measurement, scan, peptide, score)
)

matches <-
  mch %>%
  arrange(desc(score)) %>%
  group_by(scan) %>%
  mutate(score_rank = row_number(desc(score))) %>%
  filter(score_rank <= 1) %>%
  arrange(desc(score)) %>%
  summarise(pep = peptide)

okok <-
  mch %>%
  filter(isok(peptide)) %>%
  arrange(desc(score))


mch2 %>%
  group_by(measurement) %>%
  mutate(score_rank = row_number(desc(score))) %>%
  ungroup() %>%
  ggplot(aes(pepid, measurement)) +
  geom_tile(aes(fill = score)) +
  guides(x = guide_axis(angle = 90))



c <-
  mch %>%
  filter(scan == 4592) %>%
  arrange(desc(score)) %>%
  filter(peptide == "CELAAAMKR")

silico <-
  tibble(
    type = "silico", mass = str_split(c$silico_frags, ";")[[1]] %>% as.numeric(),
    intensity = 2e06
  )

nature <-
  tibble(
    type = "nature", mass = str_split(c$frags, ";")[[1]] %>% as.numeric(),
    intensity = str_split(c$frags_intensity, ";")[[1]] %>% as.numeric()
  )

silico %>%
  mutate(
    type =
      if_else(
        mass %>% map_lgl(~ min(abs(nature$mass - .x)) < 0.01),
        "silico",
        "missing"
      )
  ) %>%
  bind_rows(nature) %>%
  ggplot(aes(mass, intensity)) +
  geom_col(aes(color = type))


# MZID --------------------------------------------------------------------

t <-
  read_tsv("/Users/eugen/code/bp/bp-code/data/mgf/190318_LYS_RAT_50x_05.tsv") %>%
  mutate(
    best_match =
      gsub("\\+|\\d|,", "", str_extract(Peptide, "(?<=\\.).*(?=\\.)"))
  )

get_matches <- function(df, reference) {
  df %>%
    arrange(desc(score)) %>%
    group_by(measurement) %>%
    summarise(
      best_match = first(peptide),
      scan = first(scan),
      score = first(score),
      pepid = first(pepid)
    ) %>%
    full_join(
      reference %>%
        select(ScanNum, best_match, SpecID),
      by = c("scan" = "ScanNum")
    ) %>%
    mutate(
      is_match =
        case_when(
          is.na(best_match.x) ~ "FN",
          is.na(best_match.y) ~ "FP",
          best_match.x == best_match.y ~ "TP",
          TRUE ~ "TF"
        )
    )
}

get_matches(mch, t) %>%
  ggplot(aes(score)) +
  geom_density(aes(color = is_match))

mch %>%
  arrange(desc(score)) %>%
  group_by(measurement, scan) %>%
  summarise(best_match = first(peptide)) %>%
  full_join(t %>% select(ScanNum, best_match, SpecID), by = c("scan" = "ScanNum")) %>%
  filter(best_match.x == best_match.y) %>%
  View()
