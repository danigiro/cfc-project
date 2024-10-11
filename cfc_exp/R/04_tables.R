suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(kableExtra))
dir.create(file.path(".", "tables"), recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(TRUE)

if(length(args) == 0){
  name <- "energy"
}else{
  name <- args[1]
}

bench <- "snaive-none-base"
data <- readRDS(file.path(".", "fc", name, paste0(name, "_score.rds")))
if(name %in% c("energy")){
  comb <- c("stlf"="stlf-none-base", 
            "arima"="arima-none-base", 
            "tbats"="tbats-none-base",
            "ew"="tb+st+ar-sa-base", 
            "ow$_{\\text{var}}$"="tb+st+ar-var-base", 
            "ow$_{\\text{cov}}$"="tb+st+ar-cov0-base",
            "stlf$_{\\text{shr}}$"="stlf-cs-shr", 
            "arima$_{\\text{shr}}$"="arima-cs-shr", 
            "tbats$_{\\text{shr}}$"="tbats-cs-shr", 
            "src"="tb+st+ar-src_sa-shr", 
            "scr$_{\\text{ew}}$"="tb+st+ar-scr_sa-shr", 
            "scr$_{\\text{var}}$"="tb+st+ar-scr_var-shr", 
            "scr$_{\\text{cov}}$"="tb+st+ar-scr_cov0-shr", 
            "occ"="tb+st+ar-occ-shrbe")
  hid <- setNames(c(1:7, 0), c(1:7, "1:7")) 
  nts <- c(23, 8, 15)
  ncomb <- c(3, 3, 3, 5)
  
  comb_oa <- c("stlf"="stlf-none-base", 
               "arima"="arima-none-base", 
               "tbats"="tbats-none-base",
               "ew"="tb+st+ar-sa-base", 
               "ow$_{\\text{var}}$"="tb+st+ar-var-base", 
               "ow$_{\\text{cov}}$"="tb+st+ar-cov0-base",
               "stlf$_{\\text{shr}}$"="stlf-cs-shr", 
               "arima$_{\\text{shr}}$"="arima-cs-shr", 
               "tbats$_{\\text{shr}}$"="tbats-cs-shr", 
               "src"="tb+st+ar-src_sa-shr", 
               "scr$_{\\text{ew}}$"="tb+st+ar-scr_sa-shr", 
               "scr$_{\\text{var}}$"="tb+st+ar-scr_var-shr", 
               "scr$_{\\text{cov}}$"="tb+st+ar-scr_cov0-shr", 
               "occ$_{\\text{bv}}$"="tb+st+ar-occ-shrbv", 
               "occ$_{\\text{shr}}$"="tb+st+ar-occ-shr", 
               "occ"="tb+st+ar-occ-shrbe")
  ncomb_oa <- c(3, 3, 3, 7)
}

df <- full_join(
  data |>
    filter(name %in% comb, h %in% hid) |>
    pivot_longer(cols = c(mse, mae), names_to = "err", values_to = "value"),
  data |>
    filter(name == bench, h %in% hid) |> 
    select(-name) |>
    pivot_longer(cols = c(mse, mae), names_to = "err", values_to = "bench"),
  by = join_by(h, var, bts, err)
)

df <- df |>
  full_join(
    mutate(df, bts = -1), 
    by = join_by(h, var, name, bts, err, value, bench)
  ) |>
  group_by(h, name, err, bts) |>
  summarise(value = exp(mean(log(value/bench))), .groups = "drop") |>
  mutate(
    name = recode(name, !!!setNames(names(comb), comb)),
    name = factor(name, names(comb), ordered = TRUE),
    h = recode(h, !!!setNames(names(hid), hid)),
    h = factor(h, names(hid), ordered = TRUE)
  ) |>
  arrange(err, bts, h, name) |>
  filter(bts == -1) |>
  select(-bts) |>
  group_by(h, err) |>
  mutate(minv = min(value),
         minv2 = sort(value)[2]) |> 
  ungroup() |>
  mutate(bold = value == minv,
         italic = value == minv2,
         value = ifelse(
           value > 1,
           cell_spec(sprintf("%.4f", value), format = "latex", bold = bold, italic = italic,
                     color = "red"),
           cell_spec(sprintf("%.4f", value), format = "latex", bold = bold, italic = italic)
         )) |>
  arrange(err, h, name) |>
  select(-minv, -bold, -minv2, -italic) |>
  pivot_wider(names_from = h, values_from = value) |>
  select(-err)

df |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "lcccccccc",
      col.names = c("\\multicolumn{1}{l|}{\\textbf{Approach}}", colnames(df)[-1]),
      escape = FALSE) |>
  add_header_above(c("", "Forecast horizon" = length(hid)), 
                   escape = TRUE, bold = TRUE, line_sep = 0) |>
  pack_rows(index = setNames(rep(length(comb), 2), c(paste0("\\textbf{AvgRelMAE} - All ", 
                                                            nts[1], " time series"), 
                                                     paste0("\\textbf{AvgRelMSE} - All ", 
                                                            nts[1], " time series"))), 
            latex_gap_space = "0.3em", 
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  
  row_spec(seq(0, nrow(df), length(comb))[-c(1, nrow(df)/length(comb) +1)], 
           extra_latex_after = "\\midrule") |>
  column_spec (1, border_right = T) |>
  pack_rows(index = rep(setNames(ncomb, c("Base (incoherent forecasts)",
                                          "Combination (incoherent forecasts)",
                                          "Single model reconciliation",
                                          "Coherent combination")), 2), 
            latex_gap_space = "0em",
            indent = FALSE, escape = FALSE, latex_align = "l", bold = FALSE, italic = TRUE) |>
  save_kable(paste0("./tables/p_", name, ".tex"), self_contained = FALSE)

for(errid in c("mse", "mae")){
  
  df <- full_join(
    data |>
      filter(name %in% comb_oa, h %in% hid) |>
      pivot_longer(cols = c(mse, mae), names_to = "err", values_to = "value"),
    data |>
      filter(name == bench, h %in% hid) |> 
      select(-name) |>
      pivot_longer(cols = c(mse, mae), names_to = "err", values_to = "bench"),
    by = join_by(h, var, bts, err)
  )
  
  dfk <- df |>
    full_join(
      mutate(df, bts = -1), 
      by = join_by(h, var, name, bts, err, value, bench)
    ) |>
    group_by(h, name, err, bts) |>
    summarise(value = exp(mean(log(value/bench))), .groups = "drop") |>
    mutate(
      name = recode(name, !!!setNames(names(comb_oa), comb_oa)),
      name = factor(name, names(comb_oa), ordered = TRUE),
      h = recode(h, !!!setNames(names(hid), hid)),
      h = factor(h, names(hid), ordered = TRUE)
    ) |>
    arrange(err, bts, h, name) |>
    filter(err == errid) |>
    select(-err) |>
    group_by(h, bts) |>
    mutate(minv = min(value),
           minv2 = sort(value)[2]) |> 
    ungroup() |>
    mutate(bold = value == minv,
           italic = value == minv2,
           value = ifelse(
             value > 1,
             cell_spec(sprintf("%.4f", value), format = "latex", bold = bold, italic = italic,
                       color = "red"),
             cell_spec(sprintf("%.4f", value), format = "latex", bold = bold, italic = italic)
           )) |>
    arrange(bts, h, name) |>
    select(-minv, -bold, -minv2, -italic) |>
    pivot_wider(names_from = h, values_from = value) |>
    select( -bts)
  
  dfk |>
    kbl(format = "latex", digits = 3, booktabs = TRUE, 
        linesep = "",
        align = "lcccccccc",
        col.names = c("\\multicolumn{1}{l|}{\\textbf{Approach}}", colnames(dfk)[-1]),
        escape = FALSE) |>
    add_header_above(c("", "Forecast horizon" = length(hid)), 
                     escape = TRUE, bold = TRUE, line_sep = 0) |>
    pack_rows(index = setNames(rep(length(comb_oa), 3), c(paste0("All ", nts[1], " time series"), 
                                                       paste0(nts[2], " upper time series"), 
                                                       paste0(nts[3], " bottom time series"))), 
              latex_gap_space = "0.3em", 
              indent = FALSE, escape = FALSE, latex_align = "c", bold = TRUE) |>
    
    row_spec(seq(0, nrow(dfk), length(comb_oa))[-c(1, nrow(dfk)/length(comb_oa) +1)], 
             extra_latex_after = "\\midrule") |>
    column_spec (1, border_right = T) |>
    pack_rows(index = rep(setNames(ncomb_oa, c("Base (incoherent forecasts)",
                                            "Combination (incoherent forecasts)",
                                            "Single model reconciliation",
                                            "Coherent combination")), 3), 
              latex_gap_space = "0em",
              indent = FALSE, escape = FALSE, latex_align = "l", bold = FALSE, italic = TRUE) |> 
    save_kable(paste0("./tables/oa_", name, "_", errid, ".tex"), self_contained = FALSE)
}

