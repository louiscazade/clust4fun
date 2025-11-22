
# library(funHDDC)
# library(funFEM)
# library(LPCM)
# library(sasfunclust)

source('Fonctions/Simulations.R')
data.simu <- generation.data(250, 100, 20, "dense", 2)
id <- rep(1:ncol(data.simu$t), each = nrow(data.simu$t))
temps <- as.vector(data.simu$t)
marqueur <- as.vector(data.simu$y)
groupes <- data.simu$groupe
plot.data(data.simu, title = "Scenario 2", legend = F, x = " ", y = " ", title.pos = 0.5)

source('Fonctions/FctSimples.R')
acpf.fct(id, temps, marqueur, k = 5)
