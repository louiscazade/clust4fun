
source('Fonctions/FctAgrégées.R')
source('Fonctions/Simulations.R')

library(doParallel)
library(foreach)

# 0 - Paramètres Curta

# args <- commandArgs(TRUE)
# job <- as.numeric(args[1])
# rep <- as.numeric(args[2])
# type <- args[3]
# nb.mesures <- as.numeric(args[4])
# scen.start <- as.numeric(argq[5])
# scen.end <- as.numeric(args[6])

# set.seed(1000 + job * 10 + rep)

# 1 - Initialisation
n <- 250
tmax <- 20
nb.mesures <- c(10, 50, 100, 250)
type <- "sparse"
scen.start <- 1
scen.end <- 4

n.cores <- 4
cl <- makeCluster(n.cores)
registerDoParallel(cl)

# 2 - Clustering 
time <- system.time({
  
  # Boucle parallèle sur les scénarios
  resultats <- foreach(scen = scen.start:scen.end,
                       .packages = c(
                         "ggplot2", "patchwork", "dplyr", "tidyr", "lcmm",
                         "splines", "fdapace", "face", "cluster", "randomForest",
                         "mclust", "kml", "TSdist", "dtwclust", "fpc",
                         "ggplotify", "viridis", "ComplexHeatmap", "grid",
                         "mvtnorm", "nlraa"),
                       .export = c("generation.data", "comparaison.simulite")) %dopar% {
                         
     res.mes <- vector('list', length(nb.mesures))
     names(res.mes) <- paste0("nmes", nb.mesures)
     
     for (nmes in nb.mesures) {
       data.simu <- generation.data(n, nmes, tmax, type, scen)
       
       id <- rep(1:ncol(data.simu$t), each = nrow(data.simu$t))
       temps <- as.vector(data.simu$t)
       marqueur <- as.vector(data.simu$y)
       groupes <- data.simu$groupe
       
       res.mes[[paste0("nmes", nmes)]] <-
         comparaison.simulite(id, temps, marqueur, groupes)}
                         
     return(res.mes)}
  
})["elapsed"]

stopCluster(cl)
names(resultats) <- paste0("scen", scen.start:scen.end)
save(resultats, time, file = "resimu.RData")
