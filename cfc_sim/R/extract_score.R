suppressPackageStartupMessages(library(tidyverse))

dir.create("./score", recursive = TRUE, showWarnings = FALSE)
folder_vec <- list.dirs("./results", recursive = FALSE)
for(folder_id in folder_vec){
  data <- NULL
  for(pattern_id in c("sett1", "sett2", "sett3", "sett4", "sett5", "sett6")){
    files <- list.files(folder_id, pattern = paste0(pattern_id), full.names = TRUE)
    cat(length(files), " ")
    df_list <- lapply(files, function(x){
      readRDS(x) |>
        pivot_longer(-c("setting", "p", "nr", "nh", "var", "sim", "test", "i")) |>
        group_by(setting, p, nr, nh, var, sim, name) |>
        summarise(mse = mean((test-value)^2),
                  mae = mean(abs(test-value)),
                  .groups = "drop")
    })
    data <- bind_rows(data, df_list)
    print(pattern_id)
  }
  saveRDS(data, file = paste0("./score/", basename(folder_id), ".rds"))
  print(basename(folder_id))
}
