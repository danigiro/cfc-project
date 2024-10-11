preprocess <- function(name){
  data <- read.csv(file.path("data", name, "data.csv"))
  data <- as.matrix(data)
  tmp <- read.csv(file.path("data", name, "agg_mat.csv"))
  agg_mat <- Matrix::Matrix(as.matrix(tmp[, -1]), sparse = TRUE)
  rownames(agg_mat) <- tmp[,1]
  info <- NULL
  if(name == "energy"){
    # h = 1 in https://github.com/PuwasalaG/Probabilistic-Forecast-Reconciliation
    info$forecast_horizon <- 7
    # L+N in https://github.com/PuwasalaG/Probabilistic-Forecast-Reconciliation
    info$fixed_length <- 112+28 
    info$nn <- NULL
    info$time_series <- "zoo"
    info$inds <- readRDS(file.path("data", name, "inds.rds"))
    info$bounds = cbind(n = which(colnames(data)%in% colnames(agg_mat)), lb = 0, ub = Inf)
  }
  
  dir.create(file.path(".", "RData"), recursive = TRUE, showWarnings = FALSE)
  save(data, agg_mat, info, file = paste0("RData/", name, ".RData"))
  invisible(NULL)
}
