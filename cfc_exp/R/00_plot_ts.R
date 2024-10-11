suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))

name <- "energy"

if(!file.exists(paste0("RData/", name, ".RData"))){
  source("R/preprocess.R")
  preprocess(name)
}

load(file.path(".", "RData", paste0(name, ".RData")))
dir.create(file.path(".", "img", paste0("ts_", name)), recursive = TRUE, showWarnings = FALSE)

namets <- c("Total"="Total",
            "non-Renewable"="NRenew",
            "Renewable"="Renew",
            "Batteries"="Batt",
            "Hydro (inc. Pumps)"="HyproP",
            "Solar"="Solar",
            "Wind"="Wind",
            "Biomass"="Bio",
            "Coal"="Coal",
            "Gas"="Gas",
            "Battery (Discharging)"="BattD",
            "Battery (Charging)"="BattC",
            "Hydro"="Hydro",
            "Pumps"="Pumps",
            "Solar (Rooftop)"="SolarR",
            "Solar (Utility)"="SolarU",
            "Distillate"="Dist",
            "Black Coal"="BlCoal",
            "Brown Coal"="BrCoal",
            "Gas (Reciprocating)"="GasR",
            "Gas (OCGT)"="GasO",
            "Gas (CCGT)"="GasC",
            "Gas (Steam)"="GasS")

plot_ts_all <- as_tibble(data) |> 
  add_column(dt = info$inds, .before = 1) |>
  pivot_longer(-dt) |>
    mutate(name = factor(name, namets, ordered = TRUE),
           name = recode(name, !!!setNames(names(namets), namets))) |>
  ggplot(aes(x = dt, col = name, y = value)) + 
  geom_line() +  
    guides(col=guide_legend(ncol =1, title = NULL))+
    labs(y = "Electricity generation", x = "Date") + 
  theme_bw()

ggsave(filename = paste0("./img/ts_", name, "/plot_all_ts_energy.pdf"), 
       plot = plot_ts_all, width = 9, height = 7, dpi = 300)

# plot_ts_bts <- as_tibble(data) |>
#   add_column(dt = info$inds, .before = 1) |>
#   pivot_longer(-dt) |>
#   filter(name %in% colnames(agg_mat))|> 
#   mutate(name = factor(name, namets, ordered = TRUE),
#          name = recode(name, !!!setNames(names(namets), namets))) |>
#   ggplot(aes(x = dt, col = name, y = value)) + 
#   geom_line() +  
#   scale_colour_hue() + 
#   scale_y_sqrt() + 
#   guides(col=guide_legend(ncol =1, title = NULL))+
#   labs(y = "Electricity generation (bottom time series)", x = "Date") + 
#   theme_bw()
# 
# ggsave(filename = paste0("./img/ts_", name, "/plot_bts_ts_energy.pdf"), 
#        plot = plot_ts_bts, width = 9, height = 7, dpi = 300)
# 
# plot_ts_uts <- as_tibble(data) |>
#   add_column(dt = info$inds, .before = 1) |>
#   pivot_longer(-dt) |>
#   filter(name %in% rownames(agg_mat))|> 
#   mutate(name = factor(name, namets, ordered = TRUE),
#          name = recode(name, !!!setNames(names(namets), namets))) |>
#   ggplot(aes(x = dt, col = name, y = value)) + 
#   geom_line() +  
#   scale_colour_hue() + 
#   #scale_y_sqrt() + 
#   guides(col=guide_legend(ncol =1, title = NULL))+
#   labs(y = "Electricity generation (upper time series)", x = "Date") + 
#   theme_bw()
# ggsave(filename = paste0("./img/ts_", name, "/plot_uts_ts_energy.pdf"), 
#        plot = plot_ts_uts, width = 9, height = 7, dpi = 300)

namets_no <- c("Total"="Total",
            "non-Renew."="NRenew",
            "Renew."="Renew",
            "Batteries"="Batt",
            "Hydro (-Pumps)"="HyproP",
            "Solar"="Solar",
            "Wind"="Wind",
            "Biomass"="Bio",
            "Coal"="Coal",
            "Gas"="Gas",
            "Battery (Discharg.)"="BattD",
            "Battery (Charg.)"="BattC",
            "Hydro"="Hydro",
            "Pumps"="Pumps",
            "Solar (Rooftop)"="SolarR",
            "Solar (Utility)"="SolarU",
            "Distillate"="Dist",
            "Black Coal"="BlCoal",
            "Brown Coal"="BrCoal",
            "Gas (Recip.)"="GasR",
            "Gas (OCGT)"="GasO",
            "Gas (CCGT)"="GasC",
            "Gas (Steam)"="GasS")
colnames(agg_mat) <- recode(colnames(agg_mat), !!!setNames(names(namets_no), namets_no))
rownames(agg_mat) <- recode(rownames(agg_mat), !!!setNames(names(namets_no), namets_no))

plot_agg_mat <- as_tibble(as.matrix(agg_mat), rownames = "row") |>
  pivot_longer(-row, names_to = "col") |>
  mutate(row = factor(row, rev(names(namets_no)), ordered = TRUE),
         col = factor(col, names(namets_no), ordered = TRUE),
         label = ifelse(value == 0, "", ifelse(value>0, "+", "-"))) |>
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = as.character(value)), color = "black") +
  geom_text(aes(label = label), color = "white", size = 3) +
  scale_fill_manual(values = c("-1" = "red", "1" = "black", "0" = "white")) +
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  coord_fixed(expand = F)+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(h = 1, v = 0.5),
        axis.text.x = element_text(angle = 90, h = 1, v = 0.5))

ggsave(filename = paste0("./img/ts_", name, "/agg_mat_small.pdf"), 
       plot = plot_agg_mat + theme(legend.position = "top"), width = 4.5, height = 4)



