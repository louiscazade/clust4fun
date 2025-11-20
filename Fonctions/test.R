
# library(funHDDC)
# library(funFEM)
# library(LPCM)
# library(sasfunclust)

source('Fonctions/Simulations.R')
data.simu <- generation.data(100, 10, 20, "sparse", 1)
data.simu <- dropout(data.simu, "incrMAR", 0.3, 1) 
id <- rep(1:ncol(data.simu$t), each = nrow(data.simu$t))
temps <- as.vector(data.simu$t)
marqueur <- as.vector(data.simu$y.do)
groupes <- data.simu$groupe

plot.data(data.simu, title = "Incr MAR", legend = F, x = " ", y = "Marqueur", title.pos = 0.5)
lcmm.fct(id, temps, marqueur, regression = "splines", lien = "splines")
dtw.fct(id, temps, marqueur, centroid = "pam")


source('Fonctions/FctSimples.R')
library(foreach)
library(doParallel)
library(dtw)

# --- Préparation des données ---
id <- as.numeric(id)
data <- data.frame(id, temps, marqueur)

# Fonction pour transformer données wide
data_wide <- pivot_wider(data, names_from = temps, values_from = marqueur) %>% 
  arrange(id) %>%
  dplyr::select(where(~ !all(is.na(.))))
data_wide <- data_wide[, as.character(sort(as.numeric(names(data_wide)[-1])))]

# --- 1. Version originale mat.dist ---
t0 <- Sys.time()
mat_orig <- mat.dist(temps, marqueur)
t1 <- Sys.time()
time_orig <- t1 - t0

# --- 2. dist() pour DTW ---
t0 <- Sys.time()
data_list <- split(data_wide, f = seq_len(nrow(data_wide)))
mat_dtw <- as.matrix(dist(data_list, method = 'DTW'))
t1 <- Sys.time()
time_dist <- t1 - t0

# --- 3. DTW parallèle ---
ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)

n <- nrow(data_wide)
pairs <- combn(n, 2)
t0 <- Sys.time()
dist_list <- foreach(p = 1:ncol(pairs), .combine='c') %dopar% {
  i <- pairs[1, p]
  j <- pairs[2, p]
  dtw::dtw(data_wide[i, ], data_wide[j, ])$distance
}
mat_dtw_par <- matrix(0, n, n)
mat_dtw_par[upper.tri(mat_dtw_par)] <- dist_list
mat_dtw_par <- mat_dtw_par + t(mat_dtw_par)
t1 <- Sys.time()
time_dtw_par <- t1 - t0
stopCluster(cl)

# --- 4. Fréchet parallèle ---
n <- length(unique(id))
data_f <- list(
  temps = split(data$temps, data$id),
  marqueur = split(data$marqueur, data$id)
)

cl <- makeCluster(ncores)
registerDoParallel(cl)

mat_frech_par <- matrix(0, n, n)
t0 <- Sys.time()
mat_frech_par[upper.tri(mat_frech_par)] <- foreach(p = 1:(n*(n-1)/2),
                                                   .combine = 'c',
                                                   .packages = 'longitudinalData') %dopar% {
                                                     comb <- combn(n, 2)[, p]
                                                     i <- comb[1]
                                                     j <- comb[2]
                                                     
                                                     distFrechet(data_f$temps[[i]],
                                                                 data_f$marqueur[[i]],
                                                                 data_f$temps[[j]],
                                                                 data_f$marqueur[[j]], 
                                                                 FrechetSumOrMax = 'sum', timeScale = 10e-6)}
mat_frech_par <- mat_frech_par + t(mat_frech_par)
t1 <- Sys.time()
time_frech_par <- t1 - t0
stopCluster(cl)

# --- Résumé des temps ---
times <- data.frame(
  method = c("mat.dist()", "dist(DTW)", "DTW parallèle", "Fréchet parallèle"),
  time_sec = c(time_orig, time_dist, time_dtw_par, time_frech_par)
)

print(times)


