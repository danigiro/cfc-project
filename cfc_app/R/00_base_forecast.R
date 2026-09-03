rm(list = ls(all = TRUE))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(forecast))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Matrix))

source("R/preprocess.R")
source("R/base_model.R")
args <- commandArgs(TRUE)

if(length(args) == 0){
  # arima or ets
  name <- "energy"
  # log or lev
  model <- "snaive"
}else{
  name <- args[1]
  model <- args[2]
}

if(!file.exists(paste0("RData/", name, ".RData"))){
  preprocess(name)
}

load(paste0("RData/", name, ".RData"))
fixed_length <- info$fixed_length
forecast_horizon <- info$forecast_horizon

dir.create(file.path(".", "fc",name, "base", model), recursive = TRUE, 
           showWarnings = FALSE)

end_traing <- NROW(data)-fixed_length

pb <- progress_bar$new(format = paste0(name, " ", model, 
                                       " - Rep. :rep [:bar] :percent in :elapsed (ETA: :eta)"),
                       total = NCOL(data)*end_traing, clear = FALSE, width= 80, 
                       show_after = 0)
for(j in 0:(end_traing-1)){
  train <- data[1:(j+fixed_length), , drop = FALSE]
  h <- min(forecast_horizon, abs(NROW(data)-j-fixed_length))
  test <- data[j+fixed_length+c(1:h), , drop = FALSE]
  if(info$time_series == "ts"){
    train <- ts(train, frequency = info$freq, start = info$start_time)
  }else if(info$time_series == "zoo"){
    require(zoo)
    train2 = zoo(train, info$inds[1:(j+fixed_length)])
  }
  obj <- NULL
  for(i in 1:NCOL(data)){
    obj[[i]] <- base_model(model = model, x = train[,i], h = h)
    
    pb$tick(tokens = list(rep = formatC(j+1, width = nchar(end_traing+1), 
                                        format = "d", flag = "0")))
  }
  names(obj) <- colnames(data)
  base <- sapply(obj, extract_bm, type = "base")
  res <- sapply(obj, extract_bm, type = "res")
  fit <- lapply(obj, extract_bm, type = "fit")
  fitfc <- lapply(obj, extract_bm, type = "fc")
  
  if(h == 1){
    base <- rbind(base)
  }
  
  itername <- paste0("ite_", formatC(j+1, width = nchar(end_traing+1), 
                                     format = "d", flag = "0"), ".RData")
  save(base, res, fit, fitfc, test,
       file = file.path(".", "fc", name, "base", model, itername))
}
warnings()
