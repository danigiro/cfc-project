suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(FoReco))
suppressPackageStartupMessages(library(FoCo2))
suppressPackageStartupMessages(library(progressr))
suppressPackageStartupMessages(library(doFuture))
registerDoFuture() ## tell foreach to use futures
plan(multisession, workers = 15)

args <- commandArgs(TRUE)

if (length(args) == 0) {
  type <- "sett1"
  p <- 4 # number of experts
  nr <- 50 # number of residuals
  nh <- 100 # number of test
  structure <- "hier"
  cor_hat <- "diag"
} else {
  type <- args[1]
  p <- as.numeric(args[2]) # number of experts
  nr <- as.numeric(args[3]) # number of residuals
  nh <- as.numeric(args[4]) # number of test
  if (length(args) <= 4) {
    structure <- "hier"
  } else {
    if (length(args) <= 5) {
      structure <- args[5]
      cor_hat <- "diag"
    } else {
      structure <- args[5]
      cor_hat <- args[6]
    }
  }
}
nsim <- 500

if (structure == "lc") {
  agg_mat <- matrix(
    c(0, 0, 1, 1, 0, -1, 1, 1, -1, -1, 1, 1),
    3,
    4,
    byrow = TRUE
  )

  colnames(agg_mat) <- c("AB", "B", "C", "D")
  rownames(agg_mat) <- c("Tot", "A", "AA")
} else {
  agg_mat <- matrix(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1), 3, 4, byrow = TRUE)
  colnames(agg_mat) <- c("AA", "AB", "BA", "BB")
  rownames(agg_mat) <- c("Tot", "A", "B")
}

dir.create(
  file.path(
    paste0("./results/unbalanced", "_", structure, "_", cor_hat),
    paste0(
      "sim_p",
      p,
      "_nr",
      nr,
      "_nh",
      nh,
      "_",
      type,
      "_",
      structure,
      "_",
      cor_hat
    )
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

df_full <- NULL
#cat("type:", type, "p:", p, "nr:", nr, "nh:", nh, "\n")
freq_mat <- matrix(c(0.59, 0.16, 0.41, 0.84), 2, 2)

unfreq_mat <- matrix(c(0.97, 0.31, 0.03, 0.69), 2, 2)
run_mc_sim <- function(P, num.iters = 50) {
  # number of possible states
  num.states <- nrow(P)

  # stores the states X_t through time
  states <- numeric(num.iters)

  # initialize variable for first state
  states[1] <- 2

  for (t in 2:num.iters) {
    # probability vector to simulate next state X_{t+1}
    p <- P[states[t - 1], ]

    ## draw from multinomial and determine state
    states[t] <- which(rmultinom(1, 1, p) == 1)
  }
  return(states - 1)
}

xs <- 1:nsim

handlers(list(
  handler_progress(
    format = "[:bar] :percent in :elapsed (ETA: :eta)",
    width = 60,
    complete = "+",
    clear = FALSE,
    enable = TRUE
  )
))

doRNG::registerDoRNG(367429)
with_progress(
  {
    progess_bar_iterator <- progressor(along = xs) ## create a 5-step progressor
    y <- foreach(
      id_sim = xs,
      .packages = c("tidyverse", "reshape2", "FoReco", "FoCo2")
    ) %dopar%
      {
        if (type == "sett1") {
          # Setting 1
          mu_y <- 0
          mu_i <- rep(0, p)
          beta_y <- rep(1, 2)
          beta_i <- matrix(1, p, 2)

          s2_y <- 1
          s2_f <- c(1, 1)
          s2_i <- rep(1, p)
          ar_f <- c(0, 0)
        } else if (type == "sett2") {
          # Setting 2
          mu_y <- 0
          mu_i <- rep(0, p)
          beta_y <- rep(1, 2)
          beta_i <- matrix(0.5, p, 2)

          s2_y <- 1
          s2_f <- c(1, 1)
          s2_i <- rep(1, p)
          ar_f <- c(0, 0)
        } else if (type == "sett3") {
          # Setting 3
          mu_y <- 0
          mu_i <- rep(0, p)
          beta_y <- rep(1, 2)
          beta_i <- matrix(0.5, p, 2)

          s2_y <- 1
          s2_f <- c(1, 1)
          s2_i <- rep(1, p)
          ar_f <- c(0.9, 0.9)
        } else if (type == "sett4") {
          # Setting 4
          mu_y <- 0
          mu_i <- rep(0, p)
          beta_y <- rep(1, 2)
          #beta_i <- matrix(rbeta(2*p, 5, 5), p, 2)
          beta_i <- matrix(rbeta(2 * p, 1, 1), p, 2)

          s2_y <- 1
          s2_f <- c(1, 1)
          s2_i <- rep(1, p)
          ar_f <- c(0, 0)
        } else if (type == "sett5") {
          # Setting 5
          mu_y <- 0
          mu_i <- rep(0, p)
          beta_y <- rep(1, 2)
          beta_i <- matrix(0.5, p, 2)

          s2_y <- 1
          s2_f <- rep(1, p)
          #s2_i <- 1/rgamma(n = p, shape = 2, rate = 2)
          s2_i <- 1 / rgamma(n = p, shape = 5, rate = 5)
          ar_f <- c(0, 0)
        } else if (type == "sett6") {
          # Setting 5
          mu_y <- 0
          mu_i <- rnorm(p)
          beta_y <- rep(1, 2)
          beta_i <- matrix(0.5, p, 2)

          s2_y <- 1
          s2_f <- rep(1, p)
          s2_i <- rep(1, p)
          ar_f <- c(0, 0)
        }

        R <- diag(1 / 2, 4)
        R[upper.tri(R)] <- runif(2 * 3, -1, 1)
        R <- R + t(R)
        Flist <- NULL
        for (j in 1:4) {
          Flist[[j]] <- sapply(1:2, function(i) {
            if (ar_f[i] != 0) {
              arima.sim(n = nr + nh, list(ar = ar_f[i]), sd = sqrt(s2_f[i]))
            } else {
              rnorm(nr + nh, 0, sqrt(s2_f[i]))
            }
          })
        }
        W <- diag(rep(sqrt(s2_y), 4)) %*% R %*% diag(rep(sqrt(s2_y), 4))
        err_y <- MASS::mvrnorm(
          n = nr + nh,
          mu = rep(0, 4),
          Sigma = Matrix::nearPD(W)$mat
        )

        bts <- sapply(1:4, function(j) {
          mu_y + Flist[[j]] %*% beta_y + err_y[, j]
        })
        y <- cbind(bts %*% t(agg_mat), bts)
        colnames(y) <- unlist(dimnames(agg_mat))

        fore <- NULL
        for (ip in 1:p) {
          fore_bts_noerr <- sapply(1:4, function(j) {
            mu_i[ip] + Flist[[j]] %*% beta_i[ip, ]
          })
          fore_noerr <- cbind(fore_bts_noerr %*% t(agg_mat), fore_bts_noerr)

          if (cor_hat != "diag") {
            Rp <- diag(1 / 2, 7)
            Rp[upper.tri(Rp)] <- runif(7 * 3, -1, 1)
            Rp <- Rp + t(Rp)
          } else {
            Rp <- diag(1, 7)
          }

          if (structure == "lc") {
            s2_err <- s2_i[ip] * diag(c(3, 2, 1, 1, 1, 1, 1))
          } else {
            s2_err <- diag(c(
              as.vector(rep(s2_i[ip], 4) %*% t(agg_mat)),
              rep(s2_i[ip], 4)
            ))
          }
          Sigma_err <- nearPD(sqrt(s2_err) %*% Rp %*% sqrt(s2_err))$mat

          fore[[ip]] <- fore_noerr +
            MASS::mvrnorm(n = nr + nh, mu = rep(0, 7), Sigma = Sigma_err)
          colnames(fore[[ip]]) <- unlist(dimnames(agg_mat))
        }

        mat_freq <- sapply(1:sum(dim(agg_mat)), function(x) {
          1:p %in% sample(1:p, size = round(p * 40 / 100))
        })
        colnames(mat_freq) <- 1:sum(dim(agg_mat))
        rownames(mat_freq) <- paste0("p", 1:p)

        mat_na <- lapply(1:p, function(x) {
          sapply(1:NCOL(mat_freq), function(z) {
            if (mat_freq[x, z]) {
              run_mc_sim(freq_mat, nr + nh)
            } else {
              run_mc_sim(unfreq_mat, nr + nh)
            }
          })
        })

        mat_na <- simplify2array(mat_na)
        for (i in 1:sum(dim(agg_mat))) {
          matres <- apply(mat_na[, i, ], 2, zoo::rollsum, k = nr)
          id <- which(
            rowSums(matres[-NROW(matres), ] * mat_na[-c(1:nr), i, ]) <= 0
          )
          if (length(id) > 0) {
            for (k in 1:length(id)) {
              mat_na[nr + id[k], i, sample(which(matres[id[k], ] != 0), 1)] <- 1
            }
          }
        }

        fore <- lapply(1:p, function(x) {
          tmp <- fore[[x]]
          tmp[mat_na[,, x] == 0] <- NA
          tmp
        })

        combs_occ <- c("ols", "wls", "shrbe", "shrbv", "shr")
        combs_scr <- c("ols", "wls", "shr", "none")
        df_test <- reshape2::melt(
          y[-c(1:nr), , drop = FALSE],
          value.name = "test"
        ) |>
          as_tibble() |>
          rename(i = Var1, var = Var2) |>
          mutate(sim = id_sim)
        df_base <- bind_rows(lapply(1:length(fore), function(x) {
          reshape2::melt(fore[[x]][-c(1:nr), , drop = FALSE]) |>
            as_tibble() |>
            rename(i = Var1, var = Var2) |>
            mutate(sim = id_sim, name = paste0("p", x, "-none"))
        })) |>
          pivot_wider()

        df_occ <- NULL
        df_src <- NULL
        df_scr <- NULL
        for (i in 1:nh) {
          res <- lapply(fore, function(x) {
            x[i:(i + nr - 1), , drop = FALSE] -
              y[i:(i + nr - 1), , drop = FALSE]
          })
          base <- lapply(fore, function(x) {
            x[nr + i, , drop = FALSE]
          })

          base_na <- lapply(1:p, function(x) {
            id <- apply(res[[x]], 2, function(y) {
              all(is.na(y))
            })
            out <- base[[x]]
            if (any(id)) {
              out[, id] <- NA
            }
            out
          })

          # OCC
          occ <- NULL
          for (comb in combs_occ) {
            if (comb %in% c("ols", "str")) {
              base_tmp <- base
            } else {
              base_tmp <- base_na
            }
            reco <- suppressWarnings(FoCo2::csocc(
              base = base_tmp,
              agg_mat = agg_mat,
              res = res,
              comb = comb
            ))
            if (is.null(occ)) {
              occ <- reshape2::melt(reco, value.name = paste0("occ-", comb)) |>
                as_tibble() |>
                add_column(i = i, sim = id_sim, .before = 1) |>
                rename(var = Var2) |>
                select(-Var1)
            } else {
              occ <- left_join(
                occ,
                reshape2::melt(reco, value.name = paste0("occ-", comb)) |>
                  as_tibble() |>
                  add_column(i = i, sim = id_sim, .before = 1) |>
                  rename(var = Var2) |>
                  select(-Var1),
                by = join_by(i, sim, var)
              )
            }
          }
          df_occ <- bind_rows(df_occ, occ)

          # SCR
          scr <- NULL
          for (comb in combs_scr) {
            for (app in c("sa", "var", "cov")) {
              reco <- suppressWarnings(FoCo2::csscr(
                base = base,
                res = res,
                agg_mat = agg_mat,
                fc = app,
                nnw = TRUE,
                comb = comb,
                shrink = FALSE
              ))

              if (is.null(scr)) {
                scr <- reshape2::melt(
                  reco,
                  value.name = paste0("scr_", app, "-", comb)
                ) |>
                  as_tibble() |>
                  add_column(i = i, sim = id_sim, .before = 1) |>
                  rename(var = Var2) |>
                  select(-Var1)
              } else {
                scr <- left_join(
                  scr,
                  reshape2::melt(
                    reco,
                    value.name = paste0(
                      ifelse(comb == "none", "", "scr_"),
                      app,
                      "-",
                      comb
                    )
                  ) |>
                    as_tibble() |>
                    add_column(i = i, sim = id_sim, .before = 1) |>
                    rename(var = Var2) |>
                    select(-Var1),
                  by = join_by(i, sim, var)
                )
              }
            }
          }
          df_scr <- bind_rows(df_scr, scr)
        }

        df_reco <- full_join(df_occ, df_scr, by = join_by(i, sim, var))

        # save(df_reco, df_test, df_base, fore, y, agg_mat, mu_y, mu_i, beta_y, beta_i, mat_freq, mat_na,
        #      s2_y, s2_f, s2_i, ar_f, R, type, p, nr, nh,
        #      file = file.path(paste0("./results/unbalanced", "_", structure, "_", cor_hat),
        #                       paste0("sim_p", p, "_nr", nr, "_nh", nh, "_",
        #                              type, "_", structure, "_", cor_hat),
        #                       paste0(formatC(id_sim, width = nchar(nsim), format = "d", flag = "0"),
        #                              "_sim_p", p, "_nr", nr, "_nh", nh, "_", type, "_", structure,
        #                              "_", cor_hat, ".RData")))

        tmp_full <- full_join(
          full_join(df_test, df_base, by = join_by(i, sim, var)),
          df_reco,
          by = join_by(i, var, sim)
        )

        progess_bar_iterator() ## signal a progression update

        return(tmp_full)
      }
  },
  enable = TRUE,
  delay_stdout = TRUE,
  delay_conditions = "condition",
  cleanup = FALSE
)
warnings()

df_full <- bind_rows(y)

df_full <- df_full |>
  add_column(setting = type, p = p, nr = nr, nh = nh, .before = 1)

saveRDS(
  df_full,
  file = file.path(
    paste0("./results/unbalanced", "_", structure, "_", cor_hat),
    paste0(
      "sim_p",
      p,
      "_nr",
      nr,
      "_nh",
      nh,
      "_",
      type,
      "_",
      structure,
      "_",
      cor_hat,
      ".rds"
    )
  )
)
