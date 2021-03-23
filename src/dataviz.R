library(tidyverse)
install.packages("gridExtra")
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

mch <-
  read_csv("../out/lys_peptide_matches.csv") %>%
  rename(measurement = X1) %>%
  pivot_longer(-measurement, names_to = "peptide", values_to = "score") %>%
  mutate(peptide = factor(peptide, levels=0:max(as.numeric(peptide))))

mch %>%
  #filter(measurement < 400) %>%
  filter(score > 0.2) %>%
  ggplot(aes(peptide, measurement)) +
  geom_tile(aes(fill = score))
