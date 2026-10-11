################################################################################
# Multiplex microsphere immunoassay of goat and sheep sera, Punjab, Pakistan
#
# Cross-sectional study of influenza virus exposure in goats and sheep.
# Serum MFI values (Luminex MAGPIX) are classified as MFI-positive or negative
# per antigen with an expectation-maximization (EM) algorithm (mclust).
#
# Workflow
#   1. Import and clean the raw MFI data
#   2. Exclude sera with MFI >= 8,000 to more than one antigen
#   3. EM classification of MFI values (run once, used by all later sections)
#   4. Positive controls and MFI boxplots with EM cutoffs
#   5. Seroprevalence by district
#   6. Co-reactivity between antigens (radar plots)
#   7. Heatmaps of MFI-positive sera
#   8. MFI of positive sera, goat versus sheep (Kruskal-Wallis, Dunn test)
#   9. Percentage of MFI-positive sera per antigen
#  10. Agreement between multiplex MFI and ELISA, and ELISA seropositivity
#  11. Dose-response of selected sera
#
# Input files expected in `input_dir`:
#   AllSamples.xlsx          raw MFI values, sheets 1 to 6
#   AllSamples_PosCtl.xlsx   raw MFI values of control sera, sheets 1 to 9
#   ELISA_H5N1.xlsx, ELISA_N1.xlsx, ELISA_N2.xlsx   (columns Sample, ELISA_OD)
#   Dose_Response.xlsx
################################################################################


#### 0. Setup ##################################################################

library(mclust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(patchwork)
library(ggh4x)
library(ggnewscale)
library(ggrepel)
library(ggbreak)
library(PMCMRplus)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(grid)
library(fmsb)
library(scales)
library(stringr)

input_dir <- "input"   # folder containing the input files listed above

hosts     <- c("Sheep", "Goat")
districts <- c("Gujranwala", "Kasur", "Lahore", "Sheikhupura")

# Antigens in plotting order. Names encode subtype, strain and year.
antigens <- c("H5N1_Cam_NPH_23", "H5N8_Ast_3212_20",
              "H7N9_AH_1_13",
              "H9N2_HK_3239_08",
              "pdm_Cam_SW64_21", "hpdm_Ca_07_09",
              "H3_Cam_e0826_20", "H3_Thai_CU_P53_12",
              "N1_hpdm_Ca_04_09",
              "N2_CS1_Alb_217_17", "N2_SZN_CHN_JG20_19",
              "Vic_Au_13594_21")

# Antigens that are plotted together on the x axis
virus_group <- c("group1", "group1", "group2", "group3",
                 "group4", "group4", "group5", "group5",
                 "group6", "group7", "group7", "group8")

axis_labels <- setNames(paste0(antigens, "!", virus_group), antigens)

# Number of mixture components fitted by mclust for each antigen (order of
# `antigens`), chosen per host from the BIC.
em_components <- list(
  Sheep = c(3, 3, 3, 2, 1, 1, 1, 4, 2, 1, 3, 3),
  Goat  = c(3, 3, 3, 2, 4, 4, 1, 1, 2, 3, 3, 3)
)

exclusion_cutoff <- 8000   # MFI above which a serum is flagged as high

# Antigen order and colours used for the summary plots
plot_order <- c("H5N1_Cam_NPH_23", "H5N8_Ast_3212_20", "H7N9_AH_1_13",
                "H9N2_HK_3239_08", "hpdm_Ca_07_09", "pdm_Cam_SW64_21",
                "H3_Thai_CU_P53_12", "H3_Cam_e0826_20", "N1_hpdm_Ca_04_09",
                "N2_CS1_Alb_217_17", "N2_SZN_CHN_JG20_19", "Vic_Au_13594_21")
antigen_cols <- setNames(hue_pal()(12), plot_order)


#### 1. Import and clean raw data ##############################################

read_raw_sheet <- function(sheet_idx) {
  d <- read.xlsx(file.path(input_dir, "AllSamples.xlsx"), sheet = sheet_idx)
  colnames(d) <- gsub("[/-]", "_", colnames(d))
  d$Location <- NULL
  d %>%
    # negative or zero MFI values are set to 1
    mutate(across(NP_H3_Cam_e0826_20:hpdm_Ca_07_09, ~ifelse(. < 0, 1, .))) %>%
    mutate(across(NP_H3_Cam_e0826_20:hpdm_Ca_07_09, ~ifelse(. == 0, 1, .))) %>%
    gather(Antigen, MFI, NP_H3_Cam_e0826_20:hpdm_Ca_07_09)
}

all_samples <- do.call(rbind, lapply(1:6, read_raw_sheet))

# Sample names encode host (first letter) and district (second letter)
host_code <- substr(all_samples$Sample, 1, 1)
loc_code  <- substr(all_samples$Sample, 2, 2)
all_samples$Host <- ifelse(host_code == "S", "Sheep",
                           ifelse(host_code == "G", "Goat", host_code))
all_samples$Location <- unname(c(L = "Lahore", G = "Gujranwala",
                                 K = "Kasur", S = "Sheikhupura")[loc_code])
all_samples <- all_samples[host_code != "B", ]   # remove blank wells


#### 2. Exclude sera with MFI >= 8,000 to more than one antigen ################

excluded_samples <- lapply(setNames(hosts, hosts), function(h) {
  high <- all_samples[all_samples$Host == h & all_samples$MFI >= exclusion_cutoff, ]
  n_antigens <- tapply(high$Antigen, high$Sample, function(x) length(unique(x)))
  names(n_antigens)[n_antigens >= 2]
})

clean_samples <- all_samples[!(all_samples$Sample %in% unlist(excluded_samples)) &
                               all_samples$Antigen %in% antigens, ]

# Number of sera per host and district, before and after exclusion
count_sera <- function(df) {
  df %>%
    distinct(Sample, Host, Location) %>%
    count(Host, Location) %>%
    pivot_wider(names_from = Location, values_from = n, values_fill = 0) %>%
    mutate(Total = rowSums(across(where(is.numeric))))
}
sera_before <- count_sera(all_samples)
sera_after  <- count_sera(clean_samples)
sera_before
sera_after


#### 3. EM classification of MFI values (performed once) #######################

# Cutoff rule: the cutoff is the lowest MFI in the highest component. If the
# next lower component reaches above that lowest value (overlap), the lowest
# value is dropped first. A single component means no positives.
em_cutoff <- function(mfi, class) {
  top <- max(class)
  if (top == 1 || length(unique(class)) == 1) {
    return(list(cutoff = max(mfi), top_rows = rep(FALSE, length(mfi))))
  }
  lower <- sort(unique(class), decreasing = TRUE)[2]
  top_rows <- class == top
  if (max(mfi[class == lower]) > min(mfi[top_rows])) {
    top_rows <- top_rows & mfi != min(mfi[top_rows])
  }
  list(cutoff = min(mfi[top_rows]), top_rows = top_rows)
}

# Same rule applied to any numeric vector (used for the ELISA in section 10)
em_threshold <- function(x, G) {
  fit <- Mclust(x, G = G, verbose = FALSE)
  em_cutoff(x, fit$classification)$cutoff
}

classify_antigen <- function(d, G) {
  fit <- Mclust(d$MFI, G = G, verbose = FALSE)
  d$emClass <- fit$classification
  co <- em_cutoff(d$MFI, d$emClass)
  d$Cutoff   <- co$cutoff
  d$TopClass <- co$top_rows   # MFI-positive sera used for district counts
  d
}

em_table <- lapply(setNames(hosts, hosts), function(h) {
  host_df <- clean_samples[clean_samples$Host == h, ]
  d <- do.call(rbind, lapply(seq_along(antigens), function(i) {
    classify_antigen(host_df[host_df$Antigen == antigens[i], ], em_components[[h]][i])
  }))
  d$Gp <- factor(axis_labels[d$Antigen], levels = axis_labels)
  d$Calculation <- d$MFI - d$Cutoff
  d$Shape <- factor(ifelse(d$Calculation > 0, "16",
                           ifelse(d$Calculation < 0, "1", NA)))
  d
})

# Percentage of sera in the highest component (MFI-positive) per antigen
em_summary <- lapply(setNames(hosts, hosts), function(h) {
  do.call(rbind, lapply(antigens, function(a) {
    d <- em_table[[h]][em_table[[h]]$Antigen == a, ]
    n_pos <- sum(d$emClass == max(d$emClass))
    pct   <- round(100 * n_pos / nrow(d), 2)
    if (pct == 100) {   # one component: no positives
      n_pos <- 0
      pct <- 0
    }
    data.frame(Host = h, Antigen = a, Total = nrow(d),
               Positive = n_pos, Percentage = pct)
  }))
})

# MFI distributions with fitted components
em_histograms <- lapply(setNames(hosts, hosts), function(h) {
  lapply(antigens, function(a) {
    d <- em_table[[h]][em_table[[h]]$Antigen == a, ]
    ggplot(d, aes(x = MFI, color = as.factor(emClass))) +
      geom_histogram(aes(y = after_stat(density)), bins = 75) +
      geom_density() +
      scale_color_discrete(labels = seq_along(unique(d$emClass))) +
      ggtitle(paste0(a, "_", h)) +
      guides(color = guide_legend(title = "emClass")) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold.italic", size = 12))
  })
})
wrap_plots(em_histograms$Goat,  ncol = 4)
wrap_plots(em_histograms$Sheep, ncol = 4)

em_summary$Goat
em_summary$Sheep


#### 4. Positive controls and MFI boxplots with EM cutoffs #####################

control_raw <- NULL
for (sheet_idx in 1:9) {
  d <- read.xlsx(file.path(input_dir, "AllSamples_PosCtl.xlsx"), sheet = sheet_idx)
  colnames(d) <- gsub("[/-]", "_", colnames(d))
  d$Location <- NULL
  control_raw <- rbind(control_raw,
                       d[grep("Ferret", d$Sample), ],
                       d[grep("PBST", d$Sample), ],
                       d[grep("V", d$Sample), ])   # human convalescent sera
}
control_raw$NP_H3_Cam_e0826_20 <- NULL

control_long <- control_raw %>%
  mutate(across(H5N8_Ast_3212_20:hpdm_Ca_07_09, ~ifelse(. < 0, 1, .))) %>%
  mutate(across(H5N8_Ast_3212_20:hpdm_Ca_07_09, ~ifelse(. == 0, 1, .))) %>%
  gather(Antigen, MFI, H5N8_Ast_3212_20:hpdm_Ca_07_09)
control_long$Gp <- factor(axis_labels[control_long$Antigen], levels = axis_labels)

lbl <- gsub(" 1:500", "", control_long$Sample)
lbl[grep("_019", lbl)] <- paste0("Human-H1-", lbl[grep("_019", lbl)])
for (p in c("_020", "_146", "_198")) {
  lbl[grep(p, lbl)] <- paste0("Human-H3-", lbl[grep(p, lbl)])
}
control_long$sampleLabel <- lbl
control_long <- control_long[!grepl("_042", control_long$sampleLabel), ]  # HPIV patient

# Human control sera (not shown in the boxplots below)
control_human <- control_long[grepl("Human", control_long$sampleLabel), ]
control_human <- control_human[!grepl("Human-H3-V2_020|Human-H3-V2_146",
                                      control_human$sampleLabel), ]

# Ferret control sera: mean MFI per antigen (PBST wells excluded)
other_labels <- unique(control_long$sampleLabel[!grepl("Human", control_long$sampleLabel)])
control_means <- do.call(rbind, lapply(antigens, function(a) {
  do.call(rbind, lapply(other_labels, function(l) {
    d <- control_long[control_long$Antigen == a & control_long$sampleLabel == l, ]
    data.frame(Sample = d$Sample[1], Gp = as.character(d$Gp[1]),
               MFI = mean(d$MFI), sampleLabel = d$sampleLabel[1])
  }))
}))
control_means <- control_means[!is.na(control_means$Gp) & control_means$Sample != "PBST", ]
control_means$Gp    <- factor(control_means$Gp, levels = axis_labels)
control_means$Shape <- factor("17")

boxplot_em <- function(d, controls, host) {
  cutoffs <- distinct(d, Gp, Cutoff)
  ggplot(d) +
    geom_jitter(aes(x = Gp, y = MFI, colour = Gp, shape = Shape),
                position = position_jitter(width = 0.2), show.legend = FALSE) +
    geom_boxplot(aes(Gp, MFI), outlier.size = NA, alpha = 0, width = 0.2,
                 colour = "gray") +
    new_scale_color() +
    scale_shape_manual(values = c(1, 16, 17)) +
    geom_jitter(data = controls, aes(x = Gp, y = MFI, shape = Shape),
                position = position_jitter(width = 0.2)) +
    geom_errorbar(data = cutoffs, aes(x = Gp, ymin = Cutoff, ymax = Cutoff),
                  linetype = "solid", color = "red", linewidth = 1.2, width = 0.6) +
    scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Virus") +
    scale_y_continuous(breaks = seq(0, 15000, by = 2500)) +
    theme_bw() +
    geom_text_repel(data = controls, aes(Gp, MFI, label = sampleLabel), size = 2.3) +
    theme(legend.position = "right",
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold.italic", size = 12),
          axis.text = element_text(angle = 90, hjust = 1, size = 10)) +
    ggtitle(host)
}

fig_goat  <- boxplot_em(em_table$Goat,  control_means, "Goat")
fig_sheep <- boxplot_em(em_table$Sheep, control_means, "Sheep")
fig_goat / fig_sheep


#### 5. Seroprevalence by district #############################################

# Per host, antigen and district: MFI-positive sera (highest component),
# n tested, and the percentage of the district's sera that are positive.
district_counts <- do.call(rbind, lapply(hosts, function(h) {
  do.call(rbind, lapply(antigens, function(a) {
    d <- em_table[[h]][em_table[[h]]$Antigen == a, ]
    n_total    <- nrow(d)
    n_positive <- if (max(d$emClass) == 1) 0 else sum(d$emClass == max(d$emClass))
    do.call(rbind, lapply(districts, function(dist) {
      n_dist <- sum(d$Location == dist)
      pos    <- sum(d$TopClass & d$Location == dist)
      data.frame(Host = h, Antigen = a, Virus = axis_labels[[a]], District = dist,
                 Total = n_total, Positive = n_positive,
                 PosNr = pos, PosPerc = pos / n_total,
                 PosProp = pos / n_positive,
                 Tested_district = n_dist,
                 Pct_positive_district = 100 * pos / n_dist)
    }))
  }))
}))
district_counts$Virus <- factor(district_counts$Virus, levels = axis_labels)

plot_by_district <- function(host) {
  pd <- district_counts[district_counts$Host == host, ]
  base_theme <- theme(panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.background = element_blank(),
                      legend.title = element_blank(),
                      legend.key = element_blank(),
                      legend.background = element_blank(),
                      axis.line = element_line(),
                      panel.grid.major.x = element_blank(),
                      plot.title = element_text(face = "bold.italic", size = 12),
                      axis.text = element_text(angle = 90, hjust = 1, size = 10),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(color = "black", size = 12))
  totals <- distinct(pd, Virus, Total, Positive)
  title  <- paste0("Multiplex Serology Prevalence Plot, host = ", host)
  p1 <- ggplot(pd) +
    geom_bar(aes(x = Virus, y = PosNr, fill = District), position = "stack", stat = "identity") +
    geom_text(data = totals, aes(x = Virus, y = Total, label = Total), vjust = -1, nudge_y = 0.5) +
    scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Virus") +
    labs(y = "Number of positive samples", title = title) + base_theme
  p2 <- ggplot(pd) +
    geom_bar(aes(x = Virus, y = 100 * PosPerc, fill = District), position = "stack", stat = "identity") +
    geom_text(data = totals, aes(x = Virus, y = 100 * Positive / Total,
                                 label = paste0(round(100 * Positive / Total), "%")),
              vjust = -1, nudge_y = 0.5) +
    scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Virus") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(y = "Percentage of positive samples (%)", title = title) + base_theme
  p3 <- ggplot(pd) +
    geom_bar(aes(x = Virus, y = PosProp, fill = District), position = "stack", stat = "identity") +
    scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "Virus") +
    labs(y = "Proportion of positive samples", title = title) + base_theme
  (p1 / p2 / p3) + plot_layout(guides = "collect")
}
fig_district_sheep <- plot_by_district("Sheep")
fig_district_goat  <- plot_by_district("Goat")
fig_district_goat
fig_district_sheep


#### 6. Co-reactivity between antigens (radar plots) ###########################

# For each antigen (row), the number of its MFI-positive sera that are also
# MFI-positive for each antigen (column).
radar_antigens <- antigens[antigens != "Vic_Au_13594_21"]

coreactivity <- function(host) {
  pos <- em_table[[host]][em_table[[host]]$TopClass, ]
  mat <- matrix(0, length(radar_antigens), length(radar_antigens),
                dimnames = list(radar_antigens, radar_antigens))
  for (a in radar_antigens) {
    for (b in radar_antigens) {
      mat[a, b] <- length(intersect(pos$Sample[pos$Antigen == a],
                                    pos$Sample[pos$Antigen == b]))
    }
  }
  mat
}
coreactivity_mat <- lapply(setNames(hosts, hosts), coreactivity)

radar_chart <- function(data, color, caxislabels, title) {
  radarchart(data, axistype = 1,
             pcol = color, pfcol = alpha(color, 0.2), plwd = 2, plty = 1,
             cglcol = "grey", cglty = 1, cglwd = 0.8, axislabcol = "grey",
             vlcex = 0.7, caxislabels = caxislabels, title = title, seg = 4)
}

draw_radars <- function(host) {
  mat <- coreactivity_mat[[host]]
  top <- max(mat)
  op <- par(mar = rep(0.8, 4), mfrow = c(3, 4))
  for (a in radar_antigens) {
    d <- as.data.frame(rbind(rep(top, ncol(mat)), rep(0, ncol(mat)), mat[a, ]))
    rownames(d) <- c("Max", "Min", a)
    radar_chart(d, antigen_cols[[a]], round(seq(0, top, length.out = 5)), a)
  }
  par(op)
}
draw_radars("Goat")
draw_radars("Sheep")


#### 7. Heatmaps of MFI-positive sera ##########################################

rank_scale <- function(x) {
  r <- rank(x)
  (r - min(r)) / (max(r) - min(r))
}

positive_heatmap <- function(host) {
  d <- em_table[[host]][em_table[[host]]$Calculation >= 0, ]
  mat <- reshape2::acast(d, Sample ~ Antigen, value.var = "MFI",
                         fun.aggregate = mean, fill = 0)
  mat <- as.matrix(apply(as.data.frame(mat), 2, rank_scale))
  col_fun <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                        c("#4A4A4A", viridis(4, option = "D")))
  Heatmap(mat, name = "MFI", col = col_fun,
          cluster_rows = TRUE, cluster_columns = TRUE,
          show_row_names = FALSE, show_column_names = TRUE,
          column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
          heatmap_legend_param = list(title = "Rank(MFI)",
                                      at = c(0, 0.25, 0.5, 0.75, 1),
                                      labels = c("0", "0.25", "0.5", "0.75", "1"),
                                      legend_height = unit(40, "mm"),
                                      title_position = "topcenter"))
}
draw(positive_heatmap("Goat"))
draw(positive_heatmap("Sheep"))


#### 8. MFI of positive sera, goat versus sheep ################################

positives <- do.call(rbind, lapply(em_table, function(d) d[d$Shape %in% "16", ]))
positives$Antigen <- factor(positives$Antigen, levels = plot_order, ordered = TRUE)
positives$Host_Grouping <- paste(positives$Antigen, positives$Host, sep = "_")

fig_positives_combined <- ggplot(positives) +
  geom_boxplot(aes(Antigen, MFI, colour = Antigen, fill = Antigen),
               alpha = 0.5, width = 0.5, linewidth = 2) +
  stat_boxplot(aes(Antigen, MFI, colour = Antigen), geom = "errorbar",
               width = 0.5, linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 15000, by = 2500)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks.x = element_blank())

fig_positives_by_host <- ggplot(positives) +
  geom_boxplot(aes(Antigen, MFI, group = interaction(Antigen, Host), colour = Host),
               position = position_dodge2(preserve = "single", width = 0.8),
               alpha = 1, width = 0.8, linewidth = 0.5) +
  stat_boxplot(aes(Antigen, MFI, group = interaction(Antigen, Host), colour = Host),
               position = position_dodge2(preserve = "single", width = 0.8),
               geom = "errorbar", width = 0.8, linewidth = 1.5) +
  scale_y_continuous(breaks = seq(0, 15000, by = 2500)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1, vjust = 1))

fig_positives_combined / fig_positives_by_host

# Kruskal-Wallis with Dunn post hoc test (FDR-adjusted P values) across all
# antigen and host groups
dunn <- kwAllPairsDunnTest(positives$MFI, factor(positives$Host_Grouping),
                           p.adjust.method = "fdr")
dunn_pvalues <- dunn$p.value
dunn_pvalues


#### 9. Percentage of MFI-positive sera per antigen ############################

percentage_barplot <- function(host, break_range, tick_labels, margin) {
  d <- em_summary[[host]]
  d$Antigen <- factor(d$Antigen, levels = plot_order, ordered = TRUE)
  d$Group   <- virus_group[match(d$Antigen, antigens)]
  ggplot(d, aes(Antigen, y = Percentage, fill = Antigen)) +
    geom_bar(stat = "identity", position = "stack", colour = "black") +
    scale_y_break(break_range, space = 0.4, ticklabels = tick_labels, expand = c(0, 0)) +
    labs(x = "", y = "Percentage (%)") +
    scale_fill_manual(values = antigen_cols) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    theme(panel.background = element_rect(fill = "white"),
          axis.line.y = element_line(linetype = 1, linewidth = 1.5, colour = "black"),
          axis.line.x = element_line(linetype = 1, linewidth = 1.5, colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold", size = 15, colour = "black"),
          axis.title.y = element_text(face = "bold", size = 15, angle = 90, hjust = 0.65),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          legend.position = "none",
          plot.margin = margin)
}
fig_pct_sheep <- percentage_barplot("Sheep", c(55, 70), c(seq(0, 70, 5), 75),
                                    unit(c(1, 10, 1, 10), "mm"))
fig_pct_goat  <- percentage_barplot("Goat", c(35, 60), c(seq(0, 70, 5), 70),
                                    unit(c(1, 1, 1, 10), "mm"))
fig_pct_goat
fig_pct_sheep


#### 10. Multiplex MFI versus ELISA ############################################

# Uses the raw MFI of all sera that have an ELISA result (no exclusion step).
elisa_compare <- function(antigen, elisa_file, plot_title, xlim = NULL, ylim = NULL) {
  elisa <- read.xlsx(file.path(input_dir, elisa_file))
  by_host <- lapply(setNames(c("Goat", "Sheep"), c("Goat", "Sheep")), function(h) {
    d <- all_samples[all_samples$Host == h & all_samples$Antigen == antigen,
                     c("Sample", "MFI")]
    na.omit(merge(d, elisa[, c("Sample", "ELISA_OD")], by = "Sample", all.x = TRUE))
  })
  combined <- rbind(by_host$Goat, by_host$Sheep)

  # Pearson correlation by host and combined
  for (d in list(by_host$Goat, by_host$Sheep, combined)) {
    print(cor.test(d$MFI, d$ELISA_OD, method = "pearson"))
  }

  # Linear regression, raw and log10 scale (combined hosts)
  print(summary(lm(MFI ~ ELISA_OD, data = combined)))
  fit <- summary(lm(log10(MFI) ~ log10(ELISA_OD), data = combined))
  print(fit)
  f_stat <- fit$fstatistic[1]
  p_val  <- pf(f_stat, fit$fstatistic[2], fit$fstatistic[3], lower.tail = FALSE)

  p <- ggplot(combined, aes(x = log10(ELISA_OD), y = log10(MFI))) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, colour = "red") +
    theme_classic() +
    labs(title = plot_title,
         subtitle = sprintf("Adj R² = %.3f | F = %.2f | p = %.3g",
                            fit$adj.r.squared, f_stat, p_val),
         x = "log10(ELISA OD)", y = "log10(MFI)")
  if (!is.null(xlim)) p <- p + scale_x_continuous(limits = xlim)
  if (!is.null(ylim)) p <- p + scale_y_continuous(limits = ylim)
  list(data = by_host, plot = p)
}

h5_elisa <- elisa_compare("H5N1_Cam_NPH_23", "ELISA_H5N1.xlsx", "H5N1_Cam_NPH_23",
                          xlim = c(-0.4, 0.6), ylim = c(2, 4.1))
n1_elisa <- elisa_compare("N1_hpdm_Ca_04_09", "ELISA_N1.xlsx", "N1_hpdm_Ca_04_09")
n2_elisa <- elisa_compare("N2_SZN_CHN_JG20_19", "ELISA_N2.xlsx", "N2_SZN_CHN_JG20_19")

h5_elisa$plot
n1_elisa$plot | n2_elisa$plot

# ELISA seropositivity for H5N1: cutoff from the EM rule applied to ELISA OD
# (up to 3 components for goats, up to 4 for sheep), exact binomial 95% CI
elisa_seropositivity <- function(host, G) {
  d <- h5_elisa$data[[host]]
  cutoff <- em_threshold(d$ELISA_OD, G)
  n_pos  <- sum(d$ELISA_OD >= cutoff)
  ci     <- binom.test(n_pos, nrow(d))
  data.frame(Host = host, Cutoff = cutoff, Positive = n_pos, Tested = nrow(d),
             Percentage = 100 * unname(ci$estimate),
             CI_low = 100 * ci$conf.int[1], CI_high = 100 * ci$conf.int[2])
}
elisa_h5_seropositivity <- rbind(elisa_seropositivity("Goat", 1:3),
                                 elisa_seropositivity("Sheep", 1:4))
elisa_h5_seropositivity

# Exact binomial 95% CI for the multiplex MFI-positive sera (section 3)
mfi_seropositivity <- do.call(rbind, lapply(em_summary, function(d) {
  ci <- lapply(seq_len(nrow(d)), function(i) binom.test(d$Positive[i], d$Total[i]))
  d$CI_low  <- sapply(ci, function(x) 100 * x$conf.int[1])
  d$CI_high <- sapply(ci, function(x) 100 * x$conf.int[2])
  d
}))
mfi_seropositivity


#### 11. Dose-response of selected sera ########################################

dose <- read.xlsx(file.path(input_dir, "Dose_Response.xlsx"))
dose$Host     <- ifelse(substr(dose$Name, 1, 1) == "S", "Sheep", "Goat")
dose$Dilution <- as.numeric(str_extract(dose$Sample, "(?<=1:)\\d+"))

dose_plot <- function(column) {
  ggplot(dose, aes(x = Dilution, y = .data[[column]], colour = Host)) +
    geom_point(size = 3) +
    geom_line(aes(group = Name)) +
    scale_x_log10() +
    theme_classic() +
    labs(x = "Dilution (log scale)", y = "MFI", title = column)
}
dose_columns <- c("H5N8-Ast_3212_20", "H5N1-Cam_NPH_23", "H7N9-AH_1_13",
                  "H9N2-HK_3239_08", "N1-hpdm-Ca_04_09")
dose_plots <- lapply(dose_columns, dose_plot)
wrap_plots(dose_plots, ncol = 3)
