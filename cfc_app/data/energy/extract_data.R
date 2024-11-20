suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
daily <- read_csv("data/energy/daily.csv", show_col_types = FALSE)

# Code by https://github.com/PuwasalaG/Probabilistic-Forecast-Reconciliation/
alldata <- daily %>%
  head(-1)%>% #Remove last Observation
  select(date,contains(' -  GWh'))%>%
  rename_all(~gsub(' -  GWh','',.x))

# Biomass (should be always positive)
# (https://explore.openelectricity.org.au/energy/nem/?range=1y&interval=1d&view=discrete-time&group=Detailed)
#length(alldata$Biomass[alldata$Biomass<0])/length(alldata$Biomass) # 0.2650273
#summary(alldata$Biomass[alldata$Biomass<0]) # -0.06 -- -0.01
alldata$Biomass[alldata$Biomass<0] <- 0

# Distillate (should be always positive)
# (https://explore.openelectricity.org.au/energy/nem/?range=1y&interval=1d&view=discrete-time&group=Detailed)
#length(alldata$Distillate[alldata$Distillate<0])/length(alldata$Distillate) # 0.008196721
#summary(alldata$Distillate[alldata$Distillate<0]) # -0.01 -- -0.01
alldata$Distillate[alldata$Distillate<0] <- 0

alldata <- alldata%>%
  mutate(date=as.Date(date),
         Battery=rowSums(select(., contains("Battery"))),
         Gas = rowSums(select(., contains("Gas"))),
         Solar = rowSums(select(., contains("Solar"))),
         Coal = rowSums(select(., contains("Coal"))),
         `Hydro (inc. Pumps)` = Hydro + Pumps,
         Renewable=Biomass+Hydro+Solar+Wind + Pumps + Battery,
         `non-Renewable`=Coal+Distillate+Gas,
         Total=Renewable+`non-Renewable`)


# Matrix form
mat <- daily %>%
  head(-1)%>% #Remove last Observation
  select(date,contains(' -  GWh'))%>%
  rename_all(~gsub(' -  GWh','',.x))
mat$Biomass[mat$Biomass<0] <- 0
mat$Distillate[mat$Distillate<0] <- 0
mat$Pumps <- -mat$Pumps
mat$`Battery (Charging)` <- -mat$`Battery (Charging)`
inds <- mat$date
mat <- as.matrix(mat[,-1])
agg_mat = matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, # Total
                   1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1, # Renewable
                   0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, # non-Renewable
                   1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, # Solar
                   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, # Battery
                   0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, # Gas
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, # Coal
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0 # Hydro (inc. Pumps)
), nrow = 8, byrow = TRUE)
colnames(agg_mat) = c("SolarR", "SolarU", "Wind", "Hydro", "BattD", 
                      "GasR", "GasO", "GasC", "GasS", "Dist", "Bio", 
                      "BlCoal", "BrCoal", "Pumps", "BattC")
rownames(agg_mat) = c("Total", "Renew", "NRenew", "Solar", "Batt", "Gas", "Coal", "HyproP")
colnames(mat) <- colnames(agg_mat)

agg_mat <- agg_mat[c(1,3,5,2, 6, 7,4, 8), c(6,7,8,9,10,12,13, 5, 15, 1,2,3,11,4,14)]
mat <- mat[,  colnames(agg_mat)]
mat_all <- cbind(mat %*% t(agg_mat), mat)
# sum(abs(mat_all[, "Total"] - alldata$Total))
# sum(abs(mat_all[, "NRenew"] - alldata$`non-Renewable`))
# sum(abs(mat_all[, "Renew"] - alldata$`Renewable`))
# sum(abs(mat_all[, "Solar"] - alldata$Solar))
# sum(abs(mat_all[, "Batt"] - alldata$Battery))
# sum(abs(mat_all[, "Gas"] - alldata$Gas))
# sum(abs(mat_all[, "Coal"] - alldata$Coal))
# sum(abs(mat_all[, "HyproP"] - alldata$`Hydro (inc. Pumps)`))


dir.create(file.path(".", "data", "energy"), recursive = TRUE, showWarnings = FALSE)
write.csv(as.matrix(mat_all), "./data/energy/data.csv", row.names = FALSE)
write.csv(as.matrix(agg_mat), "./data/energy/agg_mat.csv")
saveRDS(inds, file = "./data/energy/inds.rds")
