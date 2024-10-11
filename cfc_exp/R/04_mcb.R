suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(forecast))
dir.create(file.path(".", "img"), recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(TRUE)

if(length(args) == 0){
  name <- "energy"
}else{
  name <- args[1]
}

if(name %in% c("energy")){
  comb <- c("stlf"="stlf-none-base", 
            "arima"="arima-none-base", 
            "tbats"="tbats-none-base",
            "ew"="tb+st+ar-sa-base", 
            "ow[var]"="tb+st+ar-var-base", 
            "ow[cov]"="tb+st+ar-cov0-base",
            "stlf[shr]"="stlf-cs-shr", 
            "arima[shr]"="arima-cs-shr", 
            "tbats[shr]"="tbats-cs-shr", 
            "src"="tb+st+ar-src_sa-shr", 
            "scr[ew]"="tb+st+ar-scr_sa-shr", 
            "scr[var]"="tb+st+ar-scr_var-shr", 
            "scr[cov]"="tb+st+ar-scr_cov0-shr", 
            "occ"="tb+st+ar-occ-shrbe")
  nts <- c(23, 8, 15)
}

data <- readRDS(file.path(".", "fc", name, paste0(name, "_mcb.rds")))
var_bts <- readRDS(file.path(".", "fc", name, paste0(name, "_score.rds"))) |>
  select(var, bts) |>
  unique() |>
  filter(bts == 1) |>
  pull(var) |>
  as.character()

nemenyi_fun <- function(data){
  nemenyi <- tsutils::nemenyi(data, plottype = "none")
  df_plot <- full_join(as_tibble(nemenyi$means, rownames = "name"), 
                       full_join(rename(as_tibble(nemenyi$means-nemenyi$cd/2, rownames = "name"), "l" = "value"),
                                 rename(as_tibble(nemenyi$means+nemenyi$cd/2, rownames = "name"), "u" = "value"), 
                                 by = "name"), by = "name") |>
    arrange(value) |>
    mutate(#name = gsub(" ", "", name),
      name_num = paste0(" - ", sprintf('"%.2f"', value))) |>
    add_column(fpval = nemenyi$fpval,
               fH = nemenyi$fH)
  df_plot$col <- df_plot$l <= df_plot$u[1]
  
  as_tibble(df_plot)
}

hmax <- max(data$h)

dfmcb <- data |>
  filter(name %in% comb, h %in% c(0, 1)) |>
  pivot_longer(-any_of(c("h", "var", "name", "bts", "ite")), names_to = "err")

dfmcb <- bind_rows(dfmcb |>
                      mutate(type = as.numeric(var %in% var_bts)),
                    dfmcb |>
                      mutate(type = -1)) |>
  select(any_of(c("h", "name", "var", "value", "ite", "err", "type"))) |>
  group_by(h, err, type) |>
  nest() |>
  mutate(data = map(data, pivot_wider, names_from = name),
         data = map(data, function(x) select(x, any_of(unname(comb)))),
         data = map(data, nemenyi_fun)) |>
  unnest(cols = c(data)) |>
  ungroup() |>
  arrange(value) |>
  mutate(pch_name_old = stringr::str_detect(name, "base"),
         pch_name = ifelse(pch_name_old, ifelse(stringr::str_detect(name, "stlf|arima|tbats"), "A", "B"), "C"),
         name = recode(factor(name), !!!setNames(names(comb), unname(comb))),
         name = paste0(name, name_num),
         err = toupper(err)) |>
  arrange(value) |>
  mutate(name = factor(name, unique(name), ordered = TRUE)) |>
  #mutate(h = factor(h, rev(input$h), ordered = TRUE)) |>
  arrange(desc(h), type) |> 
  mutate(type = factor(type, ordered = TRUE),
         type_name = recode(type, "-1" = paste0("All~", nts[1], "~ts"), 
                       "0" = paste0(nts[2], "~upper~ts"), 
                       "1" = paste0(nts[3], "~bottom~ts")),
         facet = paste0(err, "-", type_name, "-h==",
                        ifelse(h == 0, paste0("1*','*ldots*','*", hmax), h)),
         facet = factor(facet, unique(facet), ordered = TRUE)) 


for(errid in c("MSE", "MAE")){
  plot_mcb <- dfmcb |>
    filter(err == errid) |>
    ggplot() + 
    geom_rect(aes(xmin=l, xmax=u, fill = col), ymin=-Inf, ymax=Inf, alpha = 0.3, 
              data = function(x) summarise(group_by(x, facet), l = min(l), col = TRUE,
                                           u = min(u), .groups = "drop"))+
    geom_segment(aes(x = l, xend = u, yend = name, y = name), linewidth = 1) + 
    geom_point(aes(x = l, y = name), pch = "|", size = 3) + 
    geom_point(aes(x = u, y = name), pch = "|", size = 3) + 
    geom_point(aes(x = value, fill = col, y = name, pch = pch_name), size = 3) +
    geom_label(data = function(x) select(x, facet, fpval) |>
                 mutate(text = paste0("Friedman test p-value ", 
                                      ifelse(fpval<0.001, " < 0.001", round(fpval, 3)))),
               aes(x = Inf, y = -Inf, label = text), vjust = "inward", hjust = "inward", 
               size = 2.5,  label.size = NA) + 
    scale_shape_manual(values=c("A" = 24, "B" = 22, "C" = 21))+
    scale_fill_manual(values = c("TRUE" = "#a7e7a7","FALSE" = "white"))+
    facet_wrap(facet~., ncol = 2, scales = "free", dir = "v",
               labeller = label_parsed)+
    labs(y = NULL, x = NULL) + scale_y_discrete(labels = function(l) parse(text=l))+ 
    theme_minimal()+
    theme(legend.title = element_blank(),
          legend.position = "none",
          #text = element_text(size = utils$font),
          #strip.text = element_text(size = utils$font),
          legend.margin = margin())
  
  ggsave(filename = paste0("./img/oa_", name, "_", tolower(errid), "_mcb.pdf"), 
         plot = plot_mcb, width = 7, height = 9)
}

plot_mcb <- dfmcb |>
  filter(type == -1) |>
  ggplot() + 
  geom_rect(aes(xmin=l, xmax=u, fill = col), ymin=-Inf, ymax=Inf, alpha = 0.3, 
            data = function(x) summarise(group_by(x, facet), l = min(l), col = TRUE,
                                         u = min(u), .groups = "drop"))+
  geom_segment(aes(x = l, xend = u, yend = name, y = name), linewidth = 1) + 
  geom_point(aes(x = l, y = name), pch = "|", size = 3) + 
  geom_point(aes(x = u, y = name), pch = "|", size = 3) + 
  geom_point(aes(x = value, fill = col, y = name, pch = pch_name), size = 3) +
  geom_label(data = function(x) select(x, facet, fpval) |>
               mutate(text = paste0("Friedman test p-value ", 
                                    ifelse(fpval<0.001, " < 0.001", round(fpval, 3)))),
             aes(x = Inf, y = -Inf, label = text), vjust = "inward", hjust = "inward", 
             size = 2.5,  label.size = NA) + 
  scale_shape_manual(values=c("A" = 24, "B" = 22, "C" = 21))+
  scale_fill_manual(values = c("TRUE" = "#a7e7a7","FALSE" = "white"))+
  facet_wrap(facet~., ncol = 2, scales = "free", dir = "v",
             labeller = label_parsed)+
  labs(y = NULL, x = NULL) + scale_y_discrete(labels = function(l) parse(text=l))+ 
  theme_minimal()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        #text = element_text(size = utils$font),
        #strip.text = element_text(size = utils$font),
        legend.margin = margin())

ggsave(filename = paste0("./img/p_", name, "_mcb.pdf"), plot = plot_mcb, 
       width = 7, height = 5)
