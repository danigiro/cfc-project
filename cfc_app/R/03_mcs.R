suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(kableExtra))
source("./R/mcs_fun.R")
dir.create(file.path(".", "tables"), recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(TRUE)

if(length(args) == 0){
  name <- "energy"
}else{
  name <- args[1]
}

bench <- "tb+st+ar-sa-base"
data <- readRDS(file.path(".", "fc", name, paste0(name, "_score.rds")))
if(name %in% c("energy")){
  comb <- c("tbats"="tbats-none-base",
            "tbats$_{\\text{shr}}$"="tbats-cs-shr", 
            "ew"="tb+st+ar-sa-base", 
            "ow$_{\\text{var}}$"="tb+st+ar-var-base", 
            "ow$_{\\text{cov}}$"="tb+st+ar-cov0-base",
            "src"="tb+st+ar-src_sa-shr", 
            "scr$_{\\text{ew}}$"="tb+st+ar-scr_sa-shr", 
            "scr$_{\\text{var}}$"="tb+st+ar-scr_var-shr", 
            "scr$_{\\text{cov}}$"="tb+st+ar-scr_cov0-shr", 
            #"occ$_{\\text{wls}}$"="tb+st+ar-occ-wls", 
            "occ$_{\\text{be}}$"="tb+st+ar-occ-shrbe")
  hid <- setNames(c(1:7, 0), c(1:7, "1:7")) 
  nts <- c(23, 8, 15)
  ncomb <- c(2, 3, 5)
  
  comb_oa <- c("stlf"="stlf-none-base", 
               "arima"="arima-none-base", 
               "tbats"="tbats-none-base",
               "stlf$_{\\text{shr}}$"="stlf-cs-shr", 
               "arima$_{\\text{shr}}$"="arima-cs-shr", 
               "tbats$_{\\text{shr}}$"="tbats-cs-shr", 
               "ew"="tb+st+ar-sa-base", 
               "ow$_{\\text{var}}$"="tb+st+ar-var-base", 
               "ow$_{\\text{cov}}$"="tb+st+ar-cov0-base",
               "src"="tb+st+ar-src_sa-shr", 
               "scr$_{\\text{ew}}$"="tb+st+ar-scr_sa-shr", 
               "scr$_{\\text{var}}$"="tb+st+ar-scr_var-shr", 
               "scr$_{\\text{cov}}$"="tb+st+ar-scr_cov0-shr", 
               "occ$_{\\text{bv}}$"="tb+st+ar-occ-shrbv", 
               "occ$_{\\text{shr}}$"="tb+st+ar-occ-shr", 
               "occ$_{\\text{wls}}$"="tb+st+ar-occ-wls", 
               "occ$_{\\text{be}}$"="tb+st+ar-occ-shrbe")
  ncomb_oa <- c(3, 3, 3, 8)
}

# Model confidence set with tibble
MCSfun <- function(x){
  set.seed(146282)
  tmp <- estMCS(x, test = "t.max", l = 7, B = 10000)
  out <- as_tibble(tmp, rownames = "name") |>
    rename(pval = "MCS p-val",
           name = name) |>
    select(name, pval)
  out
}

data <- readRDS(file.path(".", "fc", name, paste0(name, "_mcb.rds")))
var_bts <- readRDS(file.path(".", "fc", name, paste0(name, "_score.rds"))) |>
  select(var, bts) |>
  unique() |>
  filter(bts == 1) |>
  pull(var) |>
  as.character()

dfmcs <- data |>
  filter(name %in% comb, h %in% c(0, 1)) |>
  pivot_longer(-any_of(c("h", "var", "name", "bts", "ite")), names_to = "err")

dfmcs <- bind_rows(dfmcs |>
                     mutate(type = as.numeric(var %in% var_bts)),
                   dfmcs |>
                     mutate(type = -1)) |>
  select(any_of(c("h", "name", "var", "value", "ite", "err", "type"))) |>
  group_by(h, err, type, var) |>
  nest() |>
  mutate(data = map(data, pivot_wider, names_from = name),
         data = map(data, function(x) select(x, any_of(unname(comb)))),
         data = map(data, MCSfun)) |>
  unnest(cols = c(data)) |>
  ungroup() |>
  group_by(name, err, type, h) |>
  summarise("0.05" = 100*sum(pval>0.05)/length(pval),
            "0.10" = 100*sum(pval>0.10)/length(pval),
            "0.15" = 100*sum(pval>0.15)/length(pval),
            "0.20" = 100*sum(pval>0.20)/length(pval),
            "0.25" = 100*sum(pval>0.25)/length(pval), .groups = "drop") |>
  pivot_longer(-c("name", "err", "type", "h"), names_to = "pval") |>
  mutate(pval = as.numeric(pval))

tab <- dfmcs |>
  mutate(
    h = recode(h, !!!setNames(names(hid), hid)),
    h = factor(h, names(hid), ordered = TRUE),
    name = recode(name, !!!setNames(names(comb), comb)),
    name = factor(name, names(comb), ordered = TRUE)) |>
  filter(type == -1,
         pval %in% c(0.05, 0.10, 0.20)) |>
  group_by(pval, h, type, err) |>
  mutate(ming = sort(unique(value), decreasing = TRUE)[1],
         ming2 = sort(unique(value), decreasing = TRUE)[2]) |>
  ungroup() |>
  mutate(value = cell_spec(sprintf("%.1f", value), format = "latex", bold = value == ming,
                           italic = value == ming2)) |>
  select(-ming, -ming2) |>
  arrange(err, h,name) |>
  pivot_wider(names_from = c("pval", "h")) |>
  select(-err, -type)

tmp <- unname(sapply(colnames(tab), function(x){
  tmp <- str_split(x, "_")[[1]]
  if(length(tmp)>1){
    paste0(c("$\\delta = ", "$h ="), c(100*(1-parse_number(tmp[1])), tmp[2]), c("\\%$", "$"))
  }else{
    c(" ", " ")
  }
}))

tab |>
  kbl(format = "latex", digits = 3, booktabs = TRUE, 
      linesep = "",
      align = "lcccccc",
      col.names = c("\\multicolumn{1}{l|}{\\textbf{Approach}}", tmp[1,-1]),
      escape = FALSE) |>
  add_header_above(table(tmp[2,])[unique(tmp[2,])], 
                   escape = FALSE, bold = TRUE, line_sep = 0) |>
  pack_rows(index = setNames(rep(length(comb), 2), c(paste0("Absolute loss - All ", 
                                                            nts[1], " time series"), 
                                                     paste0("Quadratic loss - All ", 
                                                            nts[1], " time series"))), 
            latex_gap_space = "0.3em", 
            indent = FALSE, escape = FALSE, latex_align = "c", bold = FALSE) |>
  
  row_spec(seq(0, nrow(tab), length(comb))[-c(1, nrow(tab)/length(comb) +1)], 
           extra_latex_after = "\\midrule") |>
  column_spec(c(1,4), border_right = T) |>
  pack_rows(index = rep(setNames(ncomb, c("Base (incoherent forecasts) and single model reconciliation",
                                          "Combination (incoherent forecasts)",
                                          #"Single model reconciliation",
                                          "Coherent combination")), 2), 
            latex_gap_space = "0em",
            indent = FALSE, escape = FALSE, latex_align = "l", bold = FALSE, italic = TRUE) |>
  save_kable(paste0("./tables/p_", name, "_mcs.tex"), self_contained = FALSE)

dfmcs <- data |>
  filter(name %in% comb_oa, h %in% c(0, 1)) |>
  pivot_longer(-any_of(c("h", "var", "name", "bts", "ite")), names_to = "err")

dfmcs <- bind_rows(dfmcs |>
                     mutate(type = as.numeric(var %in% var_bts)),
                   dfmcs |>
                     mutate(type = -1)) |>
  select(any_of(c("h", "name", "var", "value", "ite", "err", "type"))) |>
  group_by(h, err, type, var) |>
  nest() |>
  mutate(data = map(data, pivot_wider, names_from = name),
         data = map(data, function(x) select(x, any_of(unname(comb_oa)))),
         data = map(data, MCSfun)) |>
  unnest(cols = c(data)) |>
  ungroup() |>
  group_by(name, err, type, h) |>
  summarise("0.05" = 100*sum(pval>0.05)/length(pval),
            "0.10" = 100*sum(pval>0.10)/length(pval),
            "0.15" = 100*sum(pval>0.15)/length(pval),
            "0.20" = 100*sum(pval>0.20)/length(pval),
            "0.25" = 100*sum(pval>0.25)/length(pval), .groups = "drop") |>
  pivot_longer(-c("name", "err", "type", "h"), names_to = "pval") |>
  mutate(pval = as.numeric(pval))

for(errid in c("mse", "mae")){
  tab <- dfmcs |>
    mutate(
      h = recode(h, !!!setNames(names(hid), hid)),
      h = factor(h, names(hid), ordered = TRUE),
      name = recode(name, !!!setNames(names(comb_oa), comb_oa)),
      name = factor(name, names(comb_oa), ordered = TRUE)) |>
    filter(pval %in% c(0.05, 0.10, 0.20), 
           err == errid) |>
    group_by(pval, h, type, err) |>
    mutate(ming = sort(unique(value), decreasing = TRUE)[1],
           ming2 = sort(unique(value), decreasing = TRUE)[2]) |>
    ungroup() |>
    mutate(value = cell_spec(sprintf("%.1f", value), format = "latex", bold = value == ming,
                             italic = value == ming2)) |>
    select(-ming, -ming2) |>
    arrange(type, h,name) |>
    pivot_wider(names_from = c("pval", "h")) |>
    select(-err, -type)
  
  tmp <- unname(sapply(colnames(tab), function(x){
    tmp <- str_split(x, "_")[[1]]
    if(length(tmp)>1){
      paste0(c("$\\delta = ", "$h ="), c(100*(1-parse_number(tmp[1])), tmp[2]), c("\\%$", "$"))
    }else{
      c(" ", " ")
    }
  }))
  
  tab |>
    kbl(format = "latex", digits = 3, booktabs = TRUE, 
        linesep = "",
        align = "lcccccc",
        col.names = c("\\multicolumn{1}{l|}{\\textbf{Approach}}", tmp[1,-1]),
        escape = FALSE) |>
    add_header_above(table(tmp[2,])[unique(tmp[2,])], 
                     escape = FALSE, bold = TRUE, line_sep = 0) |>
    pack_rows(index = setNames(rep(length(comb_oa), 3), c(paste0("All ", nts[1], " time series"), 
                                                          paste0(nts[2], " upper time series"), 
                                                          paste0(nts[3], " bottom time series"))), 
              latex_gap_space = "0.3em", 
              indent = FALSE, escape = FALSE, latex_align = "c", bold = TRUE) |>
    
    row_spec(seq(0, nrow(tab), length(comb_oa))[-c(1, nrow(tab)/length(comb_oa) +1)], 
             extra_latex_after = "\\midrule") |>
    column_spec(c(1,4), border_right = T) |>
    pack_rows(index = rep(setNames(ncomb_oa, c("Base (incoherent forecasts)",
                                               "Single model reconciliation",
                                               "Combination (incoherent forecasts)",
                                               "Coherent combination")), 2), 
              latex_gap_space = "0em",
              indent = FALSE, escape = FALSE, latex_align = "l", bold = FALSE, italic = TRUE) |>
    save_kable(paste0("./tables/oa_", name, "_", errid, "_mcs.tex"), self_contained = FALSE)
}
