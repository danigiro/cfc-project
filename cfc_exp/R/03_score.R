suppressPackageStartupMessages(library(forecast))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(TRUE)

if(length(args) == 0){
  name <- "energy"
}else{
  name <- args[1]
}

snaive <- TRUE

load(paste0("RData/", name, ".RData"))

app <- c("opt", "sa", "var", "cov0")
path_vec <- file.path(".", "fc", name, app)
path_vec <- list.dirs(path_vec, full.names = TRUE, recursive = FALSE)
list_tibble <- lapply(path_vec, function(path){
  files <- sort(list.files(path, full.names = TRUE))
  
  tmp <- lapply(files, function(id){
    load(id, temp_env <- new.env())
    as.list(temp_env)$fcdata
  })
  bind_rows(tmp)
})
df <- list_tibble[[1]]
for(i in 2:length(list_tibble)){
  tmp <- list_tibble[[i]]
  tmp <- tmp[, c("h", "var", "ite", "test", 
                 colnames(tmp)[-which(colnames(tmp) %in% colnames(df))])]
  df <- full_join(df, tmp)
}

if(snaive){
  files <- sort(list.files(file.path(".", "fc", name, "base", "snaive"), full.names = TRUE))
  
  list_naive <- lapply(files, function(id){
    load(id, temp_env <- new.env())
    base <- temp_env$base
    rownames(base) <- 1:NROW(base)
    reshape2::melt(base, 
                   value.name = "snaive-none-base") |>
      tibble() |>
      mutate(Var1 = abs(parse_number(as.character(Var1)))) |>
      add_column(ite = parse_number(basename(id)), .before = 1) |>
      left_join(reshape2::melt(temp_env$test, 
                               value.name = "test"), by = c("Var1", "Var2")) |>
      rename("h" = Var1, "var" = Var2)
    
  })
  df_naive <- bind_rows(list_naive)
  df <- full_join(df, df_naive)
}

dfh <- df |>
  pivot_longer(-c(ite, h, var, test)) |>
  group_by(h, var, name) |>
  summarise(mse = mean((test-value)^2),
            mae = mean(abs(test-value)), 
            .groups = "drop")

df0 <- df |>
  pivot_longer(-c(ite, h, var, test)) |>
  group_by(var, name) |>
  summarise(h = 0,
            mse = mean((test-value)^2),
            mae = mean(abs(test-value)), 
            .groups = "drop") |>
  rename(var = var)

score <- bind_rows(dfh, df0) |>
  mutate(bts = as.numeric(var %in% colnames(agg_mat)))
saveRDS(score, file =  file.path(".", "fc", name, paste0(name, "_score.rds")))

dfh_ite <- df |>
  pivot_longer(-c(ite, h, var, test)) |>
  group_by(ite, h, var, name) |>
  summarise(mse = mean((test-value)^2),
            mae = mean(abs(test-value)), 
            .groups = "drop")
df0_ite <- df |>
  pivot_longer(-c(ite, h, var, test)) |>
  group_by(ite, var, name) |>
  summarise(h = 0,
            mse = mean((test-value)^2),
            mae = mean(abs(test-value)), 
            .groups = "drop") |>
  rename(var = var)
score_ite <- bind_rows(dfh_ite, df0_ite) |>
  mutate(bts = as.numeric(var %in% colnames(agg_mat)))

saveRDS(score_ite, file =  file.path(".", "fc", name, paste0(name, "_mcb.rds")))

comb <- c("snaive-none-base", 
          "stlf-none-base", 
          "arima-none-base", 
          "tbats-none-base",
          "tb+st+ar-sa-base", 
          "tb+st+ar-var-base", 
          "tb+st+ar-cov-base", 
          "tb+st+ar-cov0-base",
          "stlf-cs-shr", 
          "arima-cs-shr", 
          "tbats-cs-shr",
          "tb+st+ar-src_sa-shr", 
          "tb+st+ar-scr_sa-shr", 
          "tb+st+ar-scr_var-shr",
          "tb+st+ar-scr_cov-shr",
          "tb+st+ar-scr_cov0-shr", 
          "tb+st+ar-occ-shr", 
          "tb+st+ar-occ-shrbe", 
          "tb+st+ar-occ-shrbv")

dmtest_se <- function(dfh){
  h <- unique(dfh$h2)[1]
  dfh <- dfh[, -which(colnames(dfh)=="h2")]
  out <- NULL
  for(i in 1:NCOL(dfh)){
    for(j in 1:NCOL(dfh)){
      e1 <- pull(dfh[,i])
      e2 <- pull(dfh[,j])
      if(all(e1 - e2 < 1e-6)){
        tmp <- NULL
        tmp$p.value <- 1
        tmp$alternative <- "greater"
      }else{
        tmp = forecast::dm.test(e1, e2, alternative="greater", h = h, varestimator = "bartlett")
      }
      
      out <- bind_rows(out, tibble(e1 = colnames(dfh)[i],
                                   e2 = colnames(dfh)[j],
                                   pvalue = tmp$p.value,
                                   alternative = tmp$alternative))
    }
  }
  out
}
dfdm_se <- df |>
  pivot_longer(-c(ite, h, var, test)) |>
  filter(name %in% comb) |>
  mutate(value = value-test,
         h2 = h) |>
  select(-test) |>
  group_by(h, var) |>
  nest() |>
  mutate(data = map(data, pivot_wider, names_from = name),
         data = map(data, function(x) x[, -c(1)]),
         data = map(data, dmtest_se)) |>
  unnest(data) |>
  group_by(h, e1, e2, alternative)
saveRDS(dfdm_se, file =  file.path(".", "fc", name, paste0(name, "_se_dm_red.rds")))

dmtest_abs <- function(dfh){
  h <- unique(dfh$h2)[1]
  dfh <- dfh[, -which(colnames(dfh)=="h2")]
  out <- NULL
  for(i in 1:NCOL(dfh)){
    for(j in 1:NCOL(dfh)){
      e1 <- pull(dfh[,i])
      e2 <- pull(dfh[,j])
      if(all(e1 - e2 < 1e-6)){
        tmp <- NULL
        tmp$p.value <- 1
        tmp$alternative <- "greater"
      }else{
        tmp = forecast::dm.test(e1, e2, alternative="greater", h = h, power = 1,
                                varestimator = "bartlett")
      }
      
      out <- bind_rows(out, tibble(e1 = colnames(dfh)[i],
                                   e2 = colnames(dfh)[j],
                                   pvalue = tmp$p.value,
                                   alternative = tmp$alternative))
    }
  }
  out
}
dfdm_abs <- df |>
  pivot_longer(-c(ite, h, var, test)) |>
  filter(name %in% comb) |>
  mutate(value = value-test,
         h2 = h) |>
  select(-test) |>
  group_by(h, var) |>
  nest() |>
  mutate(data = map(data, pivot_wider, names_from = name),
         data = map(data, function(x) x[, -c(1)]),
         data = map(data, dmtest_abs)) |>
  unnest(data) |>
  group_by(h, e1, e2, alternative)
saveRDS(dfdm_abs, file =  file.path(".", "fc", name, paste0(name, "_abs_dm_red.rds")))

