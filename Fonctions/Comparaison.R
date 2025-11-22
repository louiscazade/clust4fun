
source('Fonctions/FctAgrégées.R')
source('Fonctions/Simulations.R')

# 0 - Paramètres Curta
# args <- commandArgs(TRUE)
# job <- as.numeric(args[1])
# rep <- as.numeric(args[2])
# type <- args[3]
# n <- as.numeric(args[4])
# nmes <- as.numeric(args[5])
# scen <- as.numeric(args[6])
# method <- args[7]
# set.seed(1000 + job * 10 + rep)

# 1 - Initialisation
type <- "dense"
n <- 250
nmes <- 100
scen <- 2

# 2 - Clustering 
time <- system.time({
  
  # Génération des données 
  data.simu <- generation.data(n, nmes, 20, type, scen)
  id <- rep(1:ncol(data.simu$t), each = nrow(data.simu$t))
  temps <- as.vector(data.simu$t)
  marqueur <- as.vector(data.simu$y)
  groupes <- data.simu$groupe
  
  # Comparaison des méthodes
  resultats <- comparaison.dense_lite(id, temps, marqueur, groupes, k.max = 2, parallel = F)
  
})["elapsed"]

file.name <- paste(type, scen, n, nmes, sep = "_")
save(resultats, time, file = paste0(file.name, ".RData"))
