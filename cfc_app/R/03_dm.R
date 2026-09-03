suppressPackageStartupMessages(library(tidyverse))
dir.create(file.path(".", "img"), recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(TRUE)

if(length(args) == 0){
  name <- "energy"
}else{
  name <- args[1]
}

if(name %in% c("energy")){
  comb <- c("tbats"="tbats-none-base",
            "tbats[shr]"="tbats-cs-shr", 
            "ew"="tb+st+ar-sa-base", 
            "ow[var]"="tb+st+ar-var-base", 
            "ow[cov]"="tb+st+ar-cov0-base",
            "src"="tb+st+ar-src_sa-shr", 
            "scr[ew]"="tb+st+ar-scr_sa-shr", 
            "scr[var]"="tb+st+ar-scr_var-shr", 
            "scr[cov]"="tb+st+ar-scr_cov0-shr", 
            #"occ[wls]"="tb+st+ar-occ-wls",
            "occ[be]"="tb+st+ar-occ-shrbe")
  nts <- c(23, 8, 15)
  ncomb <- c(2, 3, 6)
  
  comb_oa <- c("stlf"="stlf-none-base", 
               "arima"="arima-none-base", 
               "tbats"="tbats-none-base",
               "stlf[shr]"="stlf-cs-shr", 
               "arima[shr]"="arima-cs-shr", 
               "tbats[shr]"="tbats-cs-shr", 
               "ew"="tb+st+ar-sa-base", 
               "ow[var]"="tb+st+ar-var-base", 
               "ow[cov]"="tb+st+ar-cov0-base",
               "src"="tb+st+ar-src_sa-shr", 
               "scr[ew]"="tb+st+ar-scr_sa-shr", 
               "scr[var]"="tb+st+ar-scr_var-shr", 
               "scr[cov]"="tb+st+ar-scr_cov0-shr", 
               "occ[bv]"="tb+st+ar-occ-shrbv", 
               "occ[shr]"="tb+st+ar-occ-shr", 
               "occ[wls]"="tb+st+ar-occ-wls",
               "occ[be]"="tb+st+ar-occ-shrbe")
  ncomb_oa <- c(3, 3, 3, 8)
}

data_se <- readRDS(file.path(".", "fc", name, paste0(name, "_se_dm_red.rds")))
data_abs <- readRDS(file.path(".", "fc", name, paste0(name, "_abs_dm_red.rds")))
var_bts <- readRDS(file.path(".", "fc", name, paste0(name, "_score.rds"))) |>
  select(var, bts) |>
  unique() |>
  filter(bts == 1) |>
  pull(var) |>
  as.character()

hmax <- max(data_se$h)

data <- bind_rows(data_se |>
                    add_column(err = "mse"),
                  data_abs |>
                    add_column(err = "mae"))

df1 <- bind_rows(data |>
                   mutate(type = as.numeric(var %in% var_bts)),
                 data |>
                   mutate(type = -1)) |>
  filter(h == 1) |>
  unique() |>
  group_by(h, e1, e2, alternative, type, err) |>
  summarise(value = sum(pvalue <= 0.05), 
            n = length(pvalue),
            value_perc = 100*value/n, .groups = "drop") |>
  mutate(facet = paste0(type, "-h==",h),
         h2 = paste0("h==",h))

df0 <- bind_rows(data |>
                   mutate(type = as.numeric(var %in% var_bts)),
                 data |>
                   mutate(type = -1)) |>
  unique() |>
  group_by(e1, e2, alternative, type, err) |>
  summarise(value = sum(pvalue <= 0.05), 
            n = length(pvalue),
            value_perc = 100*value/n, .groups = "drop") |>
  mutate(facet = paste0(type, "-h==1*ldots*", hmax),
         h2 = paste0("h==1*','*ldots*','*", hmax),
         h = 0)

df <- bind_rows(df1, df0) |>
  filter((e1 %in% comb & e2 %in%comb)) |>
  mutate(
    e1 = recode(e1, !!!setNames(names(comb), comb)),
    e2 = recode(e2, !!!setNames(names(comb), comb)))

lim_leg <- ceiling(max(filter(df, type == -1)$value_perc)/100)
plot_dm <- df |>
  filter(type == -1) |>
  mutate(err = recode(err, "mse" = "Quadratic~loss", "mae" = "Absolute~loss"),
         e1 = factor(e1, names(comb), ordered = TRUE),
         e2 = factor(e2, names(comb), ordered = TRUE)) |>
  filter(h %in% c(0,1),
         e1 != e2) |>
  arrange(type, h) |>
  mutate(facet = factor(facet, levels = unique(facet), ordered = TRUE)) |>
  ggplot(aes(x = e1, y = e2)) + 
  geom_tile(aes(fill = value_perc/100), color = "black") +
  geom_text(aes(label = round(value_perc, 0)), color = "black", size = 2) +
  geom_vline(xintercept = sum(ncomb[1:2])+0.5, linewidth = 0.75) + 
  geom_hline(yintercept = sum(ncomb[1:2])+0.5, linewidth = 0.75) + 
  geom_vline(xintercept = ncomb[1]+0.5, linewidth = 0.75) + 
  geom_hline(yintercept = ncomb[1]+0.5, linewidth = 0.75) + 
  scale_fill_gradient(low = "white", 
                      high = "#a7e7a7", limits = c(0, lim_leg),
                      breaks = seq(0, lim_leg, 0.1),
                      labels = scales::label_percent())+
  facet_grid(err~h2, labeller = label_parsed)+
  guides(fill = guide_colorbar(
    frame.colour = "black", ticks.colour = "black",frame.linewidth = 0.1,
    draw.llim = FALSE,
    draw.ulim = FALSE,
    display = "gradient",
    barwidth = unit(0.35, 'cm'),
    barheight = unit(10, 'cm'),
    title.position = 'top'))+
  labs(y = "M2", x = "M1", fill = NULL, caption = "M2 (y-axis) is more accurate than M1 (x-axis), p-value = 0.05") + 
  scale_y_discrete(labels = function(l) parse(text=l))+ 
  scale_x_discrete(labels = function(l) parse(text=l))+ 
  coord_fixed(expand = F)+
  theme_minimal() + 
  theme(panel.background = element_rect(linewidth = 0.1),
        legend.title=element_text(size=8, hjust = 0.5),
        legend.text = element_text(size=6),
        legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, v = 0.5, h = 1),
        legend.margin = margin(b = -12.5))

ggsave(filename = paste0("./img/p_", name, "_dm.pdf"), plot = plot_dm, width = 7, height = 6.5)


df <- bind_rows(df1, df0) |>
  filter((e1 %in% comb_oa & e2 %in%comb_oa)) |>
  mutate(
    e1 = recode(e1, !!!setNames(names(comb_oa), comb_oa)),
    e2 = recode(e2, !!!setNames(names(comb_oa), comb_oa)))

lim_leg <- ceiling(max(filter(df, type == -1)$value_perc)/100)
plot_oa <- df |>
  filter(type == -1) |>
  mutate(err = recode(err, "mse" = "Quadratic~loss", "mae" = "Absolute~loss"),
         e1 = factor(e1, names(comb_oa), ordered = TRUE),
         e2 = factor(e2, names(comb_oa), ordered = TRUE)) |>
  filter(h %in% c(0,1),
         e1 != e2) |>
  arrange(type, h) |>
  mutate(facet = factor(facet, levels = unique(facet), ordered = TRUE)) |>
  ggplot(aes(x = e1, y = e2)) + 
  geom_tile(aes(fill = value_perc/100), color = "black") +
  geom_text(aes(label = round(value_perc, 0)), color = "black", size = 2) +
  geom_vline(xintercept = sum(ncomb_oa[1:2])+0.5, linewidth = 0.75) + 
  geom_hline(yintercept = sum(ncomb_oa[1:2])+0.5, linewidth = 0.75) + 
  geom_vline(xintercept = ncomb_oa[1]+0.5, linewidth = 0.75) + 
  geom_hline(yintercept = ncomb_oa[1]+0.5, linewidth = 0.75) + 
  geom_vline(xintercept = sum(ncomb_oa[1:3])+0.5, linewidth = 0.75) + 
  geom_hline(yintercept = sum(ncomb_oa[1:3])+0.5, linewidth = 0.75) + 
  scale_fill_gradient(low = "white", 
                      high = "#a7e7a7", limits = c(0, lim_leg),
                      breaks = seq(0, lim_leg, 0.1),
                      labels = scales::label_percent())+
  facet_grid(err~h2, labeller = label_parsed)+
  guides(fill = guide_colorbar(
    frame.colour = "black", ticks.colour = "black",frame.linewidth = 0.1,
    draw.llim = FALSE,
    draw.ulim = FALSE,
    display = "gradient",
    barheight = unit(0.35, 'cm'),
    barwidth = unit(12, 'cm'),
    title.position = 'top'))+
  labs(y = "M2", x = "M1", fill = "M2 (y-axis) is more accurate than M1 (x-axis), p-value = 0.05") + 
  scale_y_discrete(labels = function(l) parse(text=l))+ 
  scale_x_discrete(labels = function(l) parse(text=l))+ 
  coord_fixed(expand = F)+
  theme_minimal() + 
  theme(panel.background = element_rect(linewidth = 0.1),
        legend.title=element_text(size=8, hjust = 0.5),
        legend.text = element_text(size=6),
        legend.position = "top",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, v = 0.5, h = 1),
        legend.margin = margin(b = -12.5))

ggsave(filename = paste0("./img/oa_", name, "_dm.pdf"), plot = plot_oa, width = 7, height = 7)


