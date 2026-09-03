suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(kableExtra))
args <- commandArgs(TRUE)
type <- if (length(args) == 0) "cor" else match.arg(args[1], c("cor", "diag"))
dir.create("./tables", recursive = TRUE, showWarnings = FALSE)

data_avg_rel_bal <- readRDS(paste0("./score/", "balanced_hier_", type, ".rds")) |>
  pivot_longer(-c("setting", "p", "nr", "nh", "var", "sim", "name"), names_to = "err") |>
  pivot_wider() |>
  mutate(bench=`scr_sa-none`) |>
  pivot_longer(-c("setting", "p", "nr", "nh", "var", "sim", "err", "bench")) |>
  group_by(setting, p, nr, nh, name, err) |>
  summarise(value = exp(mean(log(value)-log(bench))),
            .groups = "drop")

dfbase <- data_avg_rel_bal |>
  filter(name %in% c(paste0("p", 1:100, "-none")) & (!is.na(value))) |>
  group_by(setting, p, nr, nh, err) |>
  summarise(min = min(value),
            name = name[which.min(value)], .groups = "drop") |>
  rename(value = min)

dfuni <- rbind(data_avg_rel_bal |>
                 filter(name %in% c(paste0("p", 1:100, "-shr")) & (!is.na(value))) |>
                 group_by(setting, p, nr, nh, err) |>
                 summarise(min = min(value),
                           name = "base-shr", .groups = "drop") |>
                 rename(value = min),
               dfbase |>
                 mutate(name = sub("none", "shr", name)) |>
                 select(-value) |>
                 left_join(data_avg_rel_bal, by = join_by(setting, p, nr, nh, err, name)) |>
                 mutate(name = "baseA-shr"), 
               dfbase |> 
                 mutate(name = "baseA-none"))

comb <- c("base$^{\\ast}$"="baseA-none", 
          "base$^{\\ast}_{\\text{shr}}$"="baseA-shr", 
          "base$_{\\text{shr}}$"="base-shr", 
          #"ew"="scr_sa-none", 
          "ow$_{\\text{var}}$"="scr_var-none", 
          "ow$_{\\text{cov}}$"="scr_cov-none",
          "src"="src_sa-shr", 
          "scr$_{\\text{ew}}$"="scr_sa-shr", 
          "scr$_{\\text{var}}$"="scr_var-shr", 
          "scr$_{\\text{cov}}$"="scr_cov-shr", 
          #"occ$_{\\text{wls}}$"="occ-wls", 
          "occ$_{\\text{be}}$"="occ-shrbe")

tab_bal <- rbind(data_avg_rel_bal |>
                   filter(name %in% comb),
                 dfuni) |>
  mutate(err = recode(err,  "mae" = "AvgRelMAE", "mse" = "AvgRelMSE"),
         name = recode(name, !!!setNames(names(comb), comb)),
         name = factor(name, names(comb), ordered = TRUE)) |>
  group_by(setting, p, nr, nh, err) |>
  mutate(minv = min(value),
         minv2 = sort(value)[2]) |> 
  ungroup() |>
  mutate(bold = value == minv,
         italic = value == minv2,
         value = ifelse(
           value > 1,
           cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, italic = italic,
                     color = "red"),
           cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, italic = italic)
         )) |>
  arrange(setting, err, name, p, nr) |>
  select(-minv, -bold, -minv2, -italic) |>
  select(setting, err, name, p, nr, value)

data_avg_rel_unb <- readRDS(paste0("./score/", "unbalanced_hier_", type, ".rds")) |>
  pivot_longer(-c("setting", "p", "nr", "nh", "var", "sim", "name"), names_to = "err") |>
  pivot_wider() |>
  rename(`scr_sa-none` = `sa-none`,
         `scr_cov-none` = `cov-none`,
         `scr_var-none` = `var-none`) |>
  mutate(bench=`scr_sa-none`) |>
  pivot_longer(-c("setting", "p", "nr", "nh", "var", "sim", "err", "bench")) |>
  group_by(setting, p, nr, nh, name, err) |>
  summarise(value = exp(mean(log(value)-log(bench))),
            .groups = "drop")

tab_unb <- data_avg_rel_unb |>
  filter(name %in% comb) |>
  mutate(err = recode(err,  "mae" = "AvgRelMAE", "mse" = "AvgRelMSE"),
         name = recode(name, !!!setNames(names(comb), comb)),
         name = factor(name, names(comb), ordered = TRUE)) |>
  group_by(setting, p, nr, nh, err) |>
  mutate(minv = min(value),
         minv2 = sort(value)[2]) |> 
  ungroup() |>
  mutate(bold = value == minv,
         italic = value == minv2,
         value = ifelse(
           value > 1,
           cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, italic = italic,
                     color = "red"),
           cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, italic = italic)
         )) |>
  arrange(setting, err, name, p, nr) |>
  select(-minv, -bold, -minv2, -italic) |>
  select(setting, err, name, p, nr, value)

tab <- rbind(tab_bal |>
               add_column(order = 0),
             tab_unb |>
               add_column(order = 1)) |>
  pivot_wider(names_from = c("order", "name"), values_from = value, names_sep = " ") |>
  mutate(setting = parse_number(setting))

for(pattern_id in c("AvgRelMSE", "AvgRelMAE")){
  tab_tmp <- tab |>
    filter(err == pattern_id) |>
    select(-err)
  options(knitr.kable.NA = '')
  
  rows_space <- seq(1, NROW(tab_tmp), 3)[-1]-1
  rows_line <- seq(1, NROW(tab_tmp), 9)[-1]-1
  rows_space <- rows_space[!(rows_space %in% rows_line)]
  tab_tmp |>
    group_by(setting) |>
    mutate(p = ifelse(duplicated(p), NA, p)) |>
    ungroup() |>
    mutate(setting = ifelse(duplicated(setting), NA, setting)) |>
    kbl(format = "latex", digits = 3, booktabs = TRUE, 
        linesep = "",
        align = paste0(c(rep("c", NCOL(tab_tmp))), collapse = ""),
        col.names = c("\\textbf{Sett.}", "\\textbf{p}", "\\textbf{N}", 
                      paste0("\\rotatebox{90}{", 
                             simplify2array(strsplit(colnames(tab_tmp)[-c(1:3)], " "))[2,], "}")),
        escape = FALSE) |>
    add_header_above(c(" " = 3,
                       "Balanced panel of forecasts"= length(unique(tab_bal$name)),
                       "Unbalanced panel of forecasts"= length(unique(tab_unb$name))), 
                     escape = FALSE, bold = TRUE, line_sep = 0, 
                     line = FALSE) |>
    column_spec(c(3, 3+length(unique(tab_bal$name))), border_right = T) |>
    row_spec(seq(1, NROW(tab_tmp), 9)[-1]-1, 
             extra_latex_after = "\\midrule") |>
    row_spec(rows_space, 
             extra_latex_after = "\\addlinespace[0.25em]") |>
    save_kable(paste0("./tables/p_", pattern_id, "_", type,".tex"), self_contained = FALSE)
}

comb <- c("base$^{\\ast}$"="baseA-none", 
          "base$^{\\ast}_{\\text{shr}}$"="baseA-shr", 
          "base$_{\\text{shr}}$"="base-shr", 
          #"ew"="scr_sa-none", 
          "ow$_{\\text{var}}$"="scr_var-none", 
          "ow$_{\\text{cov}}$"="scr_cov-none",
          "src"="src_sa-shr", 
          "scr$_{\\text{ew}}$"="scr_sa-shr", 
          #"scr$_{\\text{ew-wls}}$"="scr_sa-wls", 
          "scr$_{\\text{var}}$"="scr_var-shr", 
          #"scr$_{\\text{var-wls}}$"="scr_var-wls", 
          "scr$_{\\text{cov}}$"="scr_cov-shr", 
          #"scr$_{\\text{cov-wls}}$"="scr_cov-wls", 
          "occ$_{\\text{wls}}$"="occ-wls", 
          "occ$_{\\text{bv}}$"="occ-shrbv", 
          "occ$_{\\text{shr}}$"="occ-shr", 
          "occ$_{\\text{be}}$"="occ-shrbe")

for(sel in c("balanced", "unbalanced")){
  if(sel == "balanced"){
    data_avg_rel <- rbind(data_avg_rel_bal |>
                            filter(name %in% comb),
                          dfuni)
  }else{
    data_avg_rel <- data_avg_rel_unb |>
      filter(name %in% comb)
  }
  tab <- data_avg_rel |>
    mutate(err = recode(err,  "mae" = "AvgRelMAE", "mse" = "AvgRelMSE"),
           name = recode(name, !!!setNames(names(comb), comb)),
           name = factor(name, names(comb), ordered = TRUE)) |>
    group_by(setting, p, nr, nh, err) |>
    mutate(minv = min(value),
           minv2 = sort(value)[2]) |> 
    ungroup() |>
    mutate(bold = value == minv,
           italic = value == minv2,
           value = ifelse(
             value > 1,
             cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, italic = italic,
                       color = "red"),
             cell_spec(sprintf("%.3f", value), format = "latex", bold = bold, italic = italic)
           )) |>
    arrange(setting, err, name, p, nr) |>
    select(-minv, -bold, -minv2, -italic) |>
    select(setting, err, name, p, nr, value)|>
    pivot_wider(names_from = name, values_from = value) |>
    mutate(setting = parse_number(setting))
  
  for(pattern_id in c("AvgRelMSE", "AvgRelMAE")){
    tab_tmp <- tab |>
      filter(err == pattern_id) |>
      select(-err)
    options(knitr.kable.NA = '')
    
    rows_space <- seq(1, NROW(tab_tmp), 3)[-1]-1
    rows_line <- seq(1, NROW(tab_tmp), 9)[-1]-1
    rows_space <- rows_space[!(rows_space %in% rows_line)]
    tab_tmp |>
      group_by(setting) |>
      mutate(p = ifelse(duplicated(p), NA, p)) |>
      ungroup() |>
      mutate(setting = ifelse(duplicated(setting), NA, setting)) |>
      kbl(format = "latex", digits = 3, booktabs = TRUE, 
          linesep = "",
          align = paste0(c(rep("c", NCOL(tab_tmp))), collapse = ""),
          col.names = c("\\textbf{Sett.}", "\\textbf{p}", "\\textbf{N}", 
                        paste0("\\rotatebox{90}{",colnames(tab_tmp)[-c(1:3)], "}")),
          escape = FALSE) |>
      add_header_above(c("", "", "", "\\\\textbf{Approach}"= NCOL(tab_tmp)-3), 
                       escape = FALSE, bold = TRUE, line_sep = 0) |>
      column_spec(3, border_right = T) |>
      row_spec(seq(1, NROW(tab_tmp), 9)[-1]-1, 
               extra_latex_after = "\\midrule") |>
      row_spec(rows_space, 
               extra_latex_after = "\\addlinespace[0.25em]") |>
      save_kable(paste0("./tables/oa_", sel, "_", pattern_id, "_", type,".tex"), self_contained = FALSE)
  }
}
