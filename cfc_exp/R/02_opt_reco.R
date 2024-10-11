rm(list = ls(all = TRUE))
suppressPackageStartupMessages(library(forecast))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(FoReco))
suppressPackageStartupMessages(library(FoCo2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))

args <- commandArgs(TRUE)

if(length(args) == 0){
  name <- "energy"
  models <- c("tbats", "stlf", "arima")
}else{
  name <- args[1]
  models <- args[-1]
}

load(paste0("RData/", name, ".RData"))
check_model <- list.dirs(file.path("fc", name, "base"), full.names = FALSE, recursive = FALSE)
models <- models[models %in% check_model]

cat(paste0("-----------------------------------------------------------------\n", 
           "Optimal reconciliation-combination \nBf: ", 
           paste0(models, collapse = " + "),
           "\n-----------------------------------------------------------------\n"))

folder_name <- paste0(sapply(strsplit(models, "_"), function(x){
  paste0(substr(x, 1, 2), collapse = "")
}), collapse = "+")

# reconciliation
dir.create(file.path("fc", name, "opt", folder_name), recursive = TRUE, 
           showWarnings = FALSE)
dir.create(file.path("fc", name, "time"), recursive = TRUE, 
           showWarnings = FALSE)

combs <- c("ols", "str", "wls", "shr", "shrbe", "shrbv")

files <- sapply(models, function(x){
  sort(list.files(file.path("fc", name, "base", x), 
                  full.names = TRUE))
})

time_opt <- matrix(NA, NROW(files), length(combs), dimnames = list(1:NROW(files),
                                                                   combs))
pb <- progress_bar$new(format = paste0(" [:bar] :percent in :elapsed (ETA: :eta)"),
                       total = NROW(files), clear = FALSE, width= 60, show_after = 0)
for(j in 1:NROW(files)){
  tmp <- lapply(files[j, ], function(x){
    e1 <- new.env()
    load(x, e1)
    out <- as.list(e1)
    rownames(out$base) <- rownames(out$test) <- NULL
    if(name %in% c("energy")){
      id <- out$base[, colnames(agg_mat)] < 0
      out$base[, colnames(agg_mat)][id] <- 0
    }
    return(list(base = out$base, res = out$res, test = out$test))
  })
  
  base <- lapply(tmp, function(x) x$base)
  res <- lapply(tmp, function(x) x$res)
  test <- tmp[[1]]$test
  rm(tmp)
  fcdata <- NULL
  for(comb in combs){
    start <- Sys.time()
    occ <- csocc(base = base, res = res, agg_mat = agg_mat,
                 nn = info$nn, comb = comb, bounds = info$bounds)
    end <- Sys.time()
    time_opt[j, comb] <- difftime(end, start, units = "secs")
    
    reco <- reshape2::melt(occ, 
                   value.name = paste0(folder_name, "-occ-", comb)) |>
      tibble() |>
      mutate(Var1 = abs(parse_number(as.character(Var1))))
    
    if(is.null(fcdata)){
      fcdata <- reco
    }else{
      fcdata <- left_join(fcdata, reco, by=c("Var1", "Var2"))
    }
  }
  
  fcdata <- fcdata |>
    add_column(ite = j, .before = 1) |>
    left_join(reshape2::melt(test, 
                             value.name = "test"), by = c("Var1", "Var2")) |>
    rename("h" = Var1, "var" = Var2)
  
  pb$tick()
  
  itername <- basename(files[j,1])
  save(fcdata, file = file.path("fc", name, "opt", folder_name, itername))
}
warnings()

save(time_opt, file = file.path("fc", name, "time", 
                                paste0("time_opt_", folder_name, ".RData")))
