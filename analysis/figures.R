#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Figures for the M. tuberculosis DST concordance project.
#
# Five figures, each justified by the question rather than by appearance:
#
#   1. Confidence intervals      the imprecision IS the result at n=30
#   2. Confusion matrices        where the discordances actually sit
#   3. Variant count by lineage  reference bias, and the clustering prediction
#   4. Coverage over target genes the evidence behind every wild-type call
#   5. Variant heatmap           katG R463L as a lineage marker, visibly
#
# Deliberately excluded: PCA (the lineage barcode already answers that, with a
# purpose-built method), pathway analysis (unjustifiable across five genes),
# transition/transversion ratios (a QC metric bearing on nothing asked here).
#
# Usage, from the repository root:
#   Rscript analysis/figures.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

theme_set(
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(colour = "grey30", size = 9),
      plot.caption = element_text(colour = "grey45", size = 8, hjust = 0),
      strip.text = element_text(face = "bold")
    )
)

OUT <- "results/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

conc <- read_tsv("docs/concordance_primary.tsv", show_col_types = FALSE)
lin  <- read_tsv("docs/lineage_assignments.tsv", show_col_types = FALSE)
filt <- read_tsv("docs/filter_stats.tsv", show_col_types = FALSE)
vars <- read_tsv("docs/target_variants.tsv", show_col_types = FALSE)
cat_ <- read_tsv("docs/catalogue_matches.tsv", show_col_types = FALSE)
cov  <- read_tsv("docs/coverage_by_gene.tsv", show_col_types = FALSE)

depth <- read_table("docs/depth_summary.txt",
                    col_names = c("sample", "mean_depth"),
                    show_col_types = FALSE)

# Top-level lineage, with the unassigned isolate kept visible rather than
# dropped. Its absence of markers is informative (see findings/lineage).
lin <- lin %>%
  mutate(lineage_top = if_else(lineage == "unassigned",
                               "unassigned",
                               str_extract(lineage, "^[a-zA-Z]+[0-9]+")),
         lineage_top = factor(lineage_top,
                              levels = c("lineage1", "lineage2",
                                         "lineage3", "lineage4",
                                         "unassigned")))

LINCOL <- c(lineage1 = "#8c510a", lineage2 = "#d8b365",
            lineage3 = "#5ab4ac", lineage4 = "#01665e",
            unassigned = "grey70")

# ---------------------------------------------------------------------------
# Figure 1: confidence intervals
# ---------------------------------------------------------------------------
# Wilson score interval. The normal approximation gives a width of zero for
# 13/13, implying certainty from thirteen observations, and can produce bounds
# outside [0,1]. Neither is acceptable when several cells here are perfect.

wilson <- function(x, n, z = 1.959963985) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- x / n
  d <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  half <- (z / d) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(max(0, centre - half), min(1, centre + half))
}

perf <- conc %>%
  count(drug, predicted, phenotype) %>%
  pivot_wider(names_from = c(predicted, phenotype),
              values_from = n, values_fill = 0) %>%
  transmute(
    drug,
    tp = R_R, fp = R_S, fn = S_R, tn = S_S
  ) %>%
  rowwise() %>%
  mutate(
    sens = tp / (tp + fn),
    spec = tn / (tn + fp),
    sens_lo = wilson(tp, tp + fn)[1], sens_hi = wilson(tp, tp + fn)[2],
    spec_lo = wilson(tn, tn + fp)[1], spec_hi = wilson(tn, tn + fp)[2]
  ) %>%
  ungroup()

ci_long <- perf %>%
  select(drug, sens, sens_lo, sens_hi, spec, spec_lo, spec_hi,
         tp, fp, fn, tn) %>%
  pivot_longer(c(sens, spec), names_to = "metric", values_to = "estimate") %>%
  mutate(
    lo = if_else(metric == "sens", sens_lo, spec_lo),
    hi = if_else(metric == "sens", sens_hi, spec_hi),
    k  = if_else(metric == "sens", tp, tn),
    N  = if_else(metric == "sens", tp + fn, tn + fp),
    metric = recode(metric, sens = "Sensitivity", spec = "Specificity"),
    label = sprintf("%s\n%s", drug, metric),
    label = fct_rev(factor(label))
  )

p1 <- ggplot(ci_long, aes(x = estimate, y = label)) +
  geom_vline(xintercept = c(0.25, 0.5, 0.75), colour = "grey92") +
  geom_vline(xintercept = 1, colour = "grey80", linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.14,
                linewidth = 0.7, colour = "grey35",
                orientation = "y") +
  geom_point(size = 3.2, colour = "#01665e") +
  geom_text(aes(x = hi, label = sprintf("  %d/%d", k, N)),
            hjust = 0, size = 3, colour = "grey35") +
  scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(
    title = "Genotypic prediction of resistance: estimates and 95% intervals",
    subtitle = paste("Wilson score intervals.",
                     "Every lower bound falls between 0.59 and 0.77."),
    x = NULL, y = NULL,
    caption = paste(
      "The intervals, not the point estimates, are the result. This sample",
      "cannot distinguish good performance\nfrom perfect performance, and",
      "cannot establish that rifampicin prediction outperforms isoniazid",
      "prediction."
    )
  )

ggsave(file.path(OUT, "fig1_confidence_intervals.png"), p1,
       width = 7.5, height = 3.6, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figure 2: confusion matrices
# ---------------------------------------------------------------------------

cm <- conc %>%
  count(drug, predicted, phenotype) %>%
  complete(drug, predicted = c("R", "S"), phenotype = c("R", "S"),
           fill = list(n = 0)) %>%
  mutate(
    correct = predicted == phenotype,
    predicted = factor(predicted, levels = c("R", "S"),
                       labels = c("Predicted\nresistant",
                                  "Predicted\nsusceptible")),
    phenotype = factor(phenotype, levels = c("R", "S"),
                       labels = c("Phenotype R", "Phenotype S"))
  )

p2 <- ggplot(cm, aes(x = phenotype, y = fct_rev(predicted))) +
  geom_tile(aes(fill = correct), colour = "white", linewidth = 2) +
  geom_text(aes(label = n), size = 6, fontface = "bold",
            colour = "grey15") +
  facet_wrap(~ drug) +
  scale_fill_manual(values = c(`TRUE` = "#c7eae5", `FALSE` = "#f6e8c3"),
                    guide = "none") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(hjust = 1)) +
  labs(title = "Where the discordances sit",
       subtitle = "Primary analysis, quality-ranked phenotypes",
       x = NULL, y = NULL)

ggsave(file.path(OUT, "fig2_confusion_matrices.png"), p2,
       width = 7, height = 3.4, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figure 3: variant count by lineage
# ---------------------------------------------------------------------------
# Variant count against H37Rv is a proxy for phylogenetic distance from the
# reference, which is itself lineage 4. This figure shows the clustering that
# was predicted before lineage was assigned, with the lineage labels that
# explain it.

vc <- filt %>%
  left_join(lin, by = "sample") %>%
  arrange(filtered) %>%
  mutate(sample = factor(sample, levels = sample))

p3 <- ggplot(vc, aes(x = filtered, y = sample, colour = lineage_top)) +
  geom_segment(aes(x = 0, xend = filtered, yend = sample),
               colour = "grey88", linewidth = 0.4) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = LINCOL, name = NULL,
                      labels = c("Lineage 1", "Lineage 2", "Lineage 3",
                                 "Lineage 4", "Unassigned")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(axis.text.y = element_text(size = 7),
        legend.position = "top") +
  labs(
    title = "Filtered variants per isolate, against H37Rv",
    subtitle = paste("H37Rv is lineage 4. Distance from the reference tracks",
                     "phylogenetic distance."),
    x = "Filtered variants", y = NULL,
        caption = paste(
      "Grouping was predicted from these counts before lineage was assigned",
      "(see docs/findings/).\nLineages 2 and 3 overlap: both are roughly",
      "equidistant from lineage 4, so variant count resolves\ndepth, not",
      "direction."
    )
  )

ggsave(file.path(OUT, "fig3_variants_by_lineage.png"), p3,
       width = 9, height = 7.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figure 4: coverage over the target genes
# ---------------------------------------------------------------------------
# Every "no resistance mutation" claim rests on the assumption that a mutation
# would have been seen. This is the evidence for that assumption.

INH_RIF <- c("katG", "rpoB", "inhA", "fabG1", "ahpC")

cov_t <- cov %>%
  left_join(lin, by = "sample") %>%
  mutate(target = gene %in% INH_RIF,
         gene = fct_reorder(gene, mean_depth, .fun = median))

p4 <- ggplot(cov_t, aes(x = mean_depth, y = gene)) +
  geom_vline(xintercept = 30, colour = "#bf812d", linetype = "dashed") +
  geom_boxplot(aes(fill = target), outlier.shape = NA, alpha = 0.5,
               linewidth = 0.3) +
  geom_jitter(height = 0.18, size = 1, alpha = 0.55, colour = "grey25") +
  scale_fill_manual(values = c(`TRUE` = "#5ab4ac", `FALSE` = "grey85"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    annotate("text", x = 3, y = "ethA", label = "ERR8975061\nethA deleted",
           hjust = 0, size = 2.6, colour = "#8c510a", lineheight = 0.9) +
  labs(
    title = "Mean depth over resistance-associated genes",
    subtitle = paste("Each point is one isolate. Teal genes are the five",
                     "relevant to isoniazid and rifampicin."),
    x = "Mean depth", y = NULL,
    caption = paste(
      "Dashed line at 30x. All 30 isolates have 100% of katG, rpoB, inhA,",
      "fabG1 and ahpC covered at 30x or above,\nincluding 200 bp promoter",
      "flanks. Wild-type calls for the two target drugs rest on data, not on",
      "absence of data."
    )
  )

ggsave(file.path(OUT, "fig4_coverage_target_genes.png"), p4,
       width = 7, height = 4.6, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Figure 5: variant heatmap
# ---------------------------------------------------------------------------
# Synonymous variants are excluded: they cannot alter the protein and would
# triple the width for no interpretive gain.

grade_short <- function(x) {
  first <- substr(x, 1, 1)
  case_when(
    is.na(x) | x == "" ~ "not graded",
    first == "1" ~ "1 assoc",
    first == "2" ~ "2 assoc interim",
    first == "3" ~ "3 uncertain",
    first == "4" ~ "4 not assoc interim",
    first == "5" ~ "5 not assoc",
    TRUE ~ "not graded"
  )
}

hm <- vars %>%
  filter(effect != "synonymous") %>%
  left_join(cat_ %>% select(sample, pos, ref, alt, INH_grade, RIF_grade),
            by = c("sample", "pos", "ref", "alt")) %>%
  distinct(sample, notation, gene, .keep_all = TRUE) %>%
  mutate(grade = grade_short(coalesce(na_if(INH_grade, ""),
                                      na_if(RIF_grade, "")))) %>%
  left_join(lin %>% select(sample, lineage_top), by = "sample") %>%
  left_join(conc %>% filter(drug == "INH") %>%
              select(sample, INH = phenotype), by = "sample") %>%
  left_join(conc %>% filter(drug == "RIF") %>%
              select(sample, RIF = phenotype), by = "sample")

samp_order <- hm %>%
  distinct(sample, lineage_top) %>%
  arrange(lineage_top, sample) %>%
  pull(sample)

var_order <- hm %>%
  count(notation, gene) %>%
  arrange(gene, desc(n)) %>%
  pull(notation)

hm <- hm %>%
  mutate(sample = factor(sample, levels = samp_order),
         notation = factor(notation, levels = unique(var_order)))

GRADECOL <- c(
  "1 assoc"             = "#a6611a",
  "2 assoc interim"     = "#dfc27d",
  "3 uncertain"         = "grey75",
  "4 not assoc interim" = "#80cdc1",
  "5 not assoc"         = "#018571",
  "not graded"          = "grey88"
)

p5 <- ggplot(hm, aes(x = notation, y = sample, fill = grade)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = GRADECOL, name = "WHO grading") +
  facet_grid(lineage_top ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6.5),
        panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 7),
        legend.position = "top") +
  labs(
    title = "Non-synonymous and promoter variants in the five target genes",
    subtitle = "Isolates grouped by lineage. Synonymous variants excluded.",
    x = NULL, y = NULL,
    caption = paste(
      "katG R463L (graded 5, not associated with resistance) is present in",
      "every non-lineage-4 isolate and absent\nfrom every lineage-4 isolate.",
      "It is a phylogenetic marker. Presence in a resistance gene does not",
      "imply\na resistance mechanism."
    )
  )

ggsave(file.path(OUT, "fig5_variant_heatmap.png"), p5,
       width = 9, height = 7.5, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------

message("\nwrote:")
for (f in list.files(OUT, full.names = TRUE)) message("  ", f)

message("\nsanity checks:")
message("  isolates in concordance table: ",
        n_distinct(conc$sample))
message("  isolates with a lineage row:   ", n_distinct(lin$sample))
message("  isolates in coverage table:    ", n_distinct(cov$sample))
message("  variants in heatmap:           ", nrow(hm))
