rm(list = ls(all = TRUE))
suppressPackageStartupMessages(library(forecast))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(FoReco))
suppressPackageStartupMessages(library(FoCo2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))

args <- commandArgs(TRUE)

if (length(args) == 0) {
  name <- "energy1"
  type_arg <- "cov0"
  models <- c("tbats", "stlf", "arima")
} else {
  name <- args[1]
  type_arg <- args[2]
  models <- args[-c(1, 2)]
}

if (stringr::str_sub(type_arg, -1, -1) == "0") {
  type <- stringr::str_sub(type_arg, 1, -2)
  nn <- TRUE
} else {
  type <- type_arg
  nn <- FALSE
}

load(paste0("RData/", name, ".RData"))
check_model <- list.dirs(
  file.path("fc", name, "base"),
  full.names = FALSE,
  recursive = FALSE
)
models <- models[models %in% check_model]

cat(paste0(
  "-----------------------------------------------------------------\n",
  "Sequential reconciliation-combination (",
  type_arg,
  ")\nBf: ",
  paste0(models, collapse = " + "),
  "\n-----------------------------------------------------------------\n"
))

folder_name <- paste0(
  sapply(strsplit(models, "_"), function(x) {
    paste0(substr(x, 1, 2), collapse = "")
  }),
  collapse = "+"
)


# reconciliation
dir.create(
  file.path("fc", name, type_arg, folder_name),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path("fc", name, "time"),
  recursive = TRUE,
  showWarnings = FALSE
)

combs <- c("ols", "str", "wls", "shr")

files <- sapply(models, function(x) {
  sort(list.files(file.path("fc", name, "base", x), full.names = TRUE))
})

time_src <- matrix(
  NA,
  NROW(files),
  length(combs),
  dimnames = list(1:NROW(files), combs)
)
time_scr <- matrix(
  NA,
  NROW(files),
  length(combs),
  dimnames = list(1:NROW(files), combs)
)
pb <- progress_bar$new(
  format = paste0(" [:bar] :percent in :elapsed (ETA: :eta)"),
  total = NROW(files),
  clear = FALSE,
  width = 60,
  show_after = 0
)
for (j in 1:NROW(files)) {
  tmp <- lapply(files[j, ], function(x) {
    e1 <- new.env()
    load(x, e1)
    out <- as.list(e1)
    rownames(out$base) <- rownames(out$test) <- NULL
    if (name %in% c("energy")) {
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
  for (comb in combs) {
    start <- Sys.time()
    src <- FoCo2::cssrc(
      base = base,
      res = res,
      agg_mat = agg_mat,
      fc = type,
      nn = info$nn,
      nnw = nn,
      comb = comb,
      shrink = FALSE,
      bounds = info$bounds
    )
    end <- Sys.time()
    time_src[j, comb] <- difftime(end, start, units = "secs")

    start <- Sys.time()
    scr <- FoCo2::csscr(
      base = base,
      res = res,
      agg_mat = agg_mat,
      fc = type,
      nn = info$nn,
      nnw = nn,
      comb = comb,
      shrink = FALSE,
      bounds = info$bounds
    )
    end <- Sys.time()
    time_scr[j, comb] <- difftime(end, start, units = "secs")
    dimnames(scr) <- dimnames(src)

    reco <- as.list(attr(src, "FoReco"))$reco
    names(reco) <- paste0(names(res), "-cs-")
    reco[[paste0(folder_name, "-src_", type_arg, "-")]] <- src
    reco[[paste0(folder_name, "-scr_", type_arg, "-")]] <- scr

    reco <- lapply(1:length(reco), function(x) {
      reshape2::melt(
        as.matrix(reco[[x]]),
        value.name = paste0(names(reco)[x], comb)
      ) |>
        tibble() |>
        mutate(Var1 = abs(parse_number(as.character(Var1))))
    })
    reco <- Reduce(
      function(d1, d2) left_join(d1, d2, by = c("Var1", "Var2")),
      reco
    )
    if (is.null(fcdata)) {
      fcdata <- reco
    } else {
      fcdata <- left_join(fcdata, reco, by = c("Var1", "Var2"))
    }

    if (comb == combs[length(combs)]) {
      names(base) <- paste0(names(base), "-none-base")
      base[[paste0(folder_name, "-", type_arg, "-base")]] <- rbind(
        as.list(attr(scr, "FoReco"))$base
        #recoinfo(scr, verbose = FALSE)$base
      )
      dimnames(base[[paste0(
        folder_name,
        "-",
        type_arg,
        "-base"
      )]]) <- dimnames(base[[1]])
      base <- lapply(1:length(base), function(x) {
        reshape2::melt(as.matrix(base[[x]]), value.name = names(base)[x]) |>
          tibble() |>
          mutate(Var1 = abs(parse_number(as.character(Var1))))
      })
      base <- Reduce(
        function(d1, d2) left_join(d1, d2, by = c("Var1", "Var2")),
        base
      )
      fcdata <- left_join(fcdata, base, by = c("Var1", "Var2"))
    }
  }

  fcdata <- fcdata |>
    add_column(ite = j, .before = 1) |>
    left_join(
      reshape2::melt(test, value.name = "test"),
      by = c("Var1", "Var2")
    ) |>
    rename("h" = Var1, "var" = Var2)

  pb$tick()

  itername <- basename(files[j, 1])
  save(fcdata, file = file.path("fc", name, type_arg, folder_name, itername))
}
warnings()

save(
  time_src,
  time_scr,
  file = file.path(
    "fc",
    name,
    "time",
    paste0("time_", type_arg, "_", folder_name, ".RData")
  )
)
