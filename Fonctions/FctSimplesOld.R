
# Packages
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(lcmm) # LCMM
library(splines)
library(fdapace) # ACPF
library(face) # ACPF
library(cluster) # PAM
library(randomForest) # RF
library(mclust) # MM et Simulation
library(kml) # KML
library(TSdist)
library(dtwclust) # DTW
library(fpc)
library(ggplotify)
library(viridis)
library(ComplexHeatmap)
library(grid)

# Fonctions ACPF manuelle
interpol.fct <- function(y, t, tNew){
  # Function that interpolates 1d function from t grid to tNew grid
  
  yNew <- rep(NA, length(tNew))
  idx_tNew <- sapply(tNew, function(x) sum(x>=t))
  for (j in 1:length(tNew)){
    tj <- tNew[j]
    if (any(tj == t)){
      yNew[j] <- y[idx_tNew[j]]
    } else {
      idx_tauj <- idx_tNew[j]
      tauj <- t[idx_tauj]
      taujplus <- t[idx_tauj+1]
      yNew[j] <- (taujplus - tj)/(taujplus - tauj)*y[idx_tauj] + (tj - tauj)/(taujplus-tauj)*y[idx_tauj+1]
    }
  }
  return(yNew)
}
interpCovMat.fct <- function(CovMat, tMat, tNew){
  # Function that interpolates CovMat from tMat grid to tNew grid
  # add pre-treatment to remove tNew not in the range of tMat
  mat_tNew <- matrix(NA,nrow=length(tNew), ncol=length(tNew))
  idx_tNew <- sapply(tNew, function(x) sum(x>=tMat))
  
  for (j in 1:length(tNew)){
    
    tj <- tNew[j]
    idx_tauj <- idx_tNew[j]
    tauj <- tMat[idx_tauj]
    taujplus <- tMat[idx_tauj+1]
    # on parcourt la triangulaire inférieure
    for (k in 1:j){
      tk <- tNew[k]
      idx_tauk <- idx_tNew[k]
      tauk <- tMat[idx_tauk]
      taukplus <-  tMat[idx_tauk+1]
      # matrix interpolation
      Gjk <- CovMat[idx_tauj,idx_tauk] +
        (tj-tauj)/(taujplus-tauj)*(CovMat[idx_tauj+1, idx_tauk]-CovMat[idx_tauj,idx_tauk]) +
        (tj-tauj)*(tk-tauk)/((taujplus-tauj)*(taukplus-tauk))*
        (CovMat[idx_tauj+1, idx_tauk+1] - CovMat[idx_tauj+1, idx_tauk] - CovMat[idx_tauj,idx_tauk+1] + CovMat[idx_tauj,idx_tauk]) +
        (tk-tauk)/(taukplus-tauk)*(CovMat[idx_tauj,idx_tauk+1]-CovMat[idx_tauj,idx_tauk])
      if (j != k) {
        mat_tNew[j,k] <- mat_tNew[k,j] <- Gjk
      } else {
        mat_tNew[j,j] <- Gjk
      }
    }
  }
  return(mat_tNew)
}
scores.acpf <- function(FPCAobj, marqueur, temps, grilles.temps, type = 'fdapace'){
  
  if (type == 'fdapace'){
    
    scores <- matrix(NA, nrow = length(marqueur), ncol = FPCAobj$selectK)
    min_grilles.temps <- min(unlist(grilles.temps), na.rm = TRUE)
    max_grilles.temps <- max(unlist(grilles.temps), na.rm = TRUE)
    
    for (i in 1:length(marqueur)){
      
      Lti <- temps[[i]][temps[[i]]>=min_grilles.temps & temps[[i]]<max_grilles.temps]
      Lyi <- marqueur[[i]][temps[[i]]>=min_grilles.temps & temps[[i]]<max_grilles.temps]
      Lti <- Lti[!is.na(Lyi)]; Lyi <- Lyi[!is.na(Lyi)];
      
      if (length(Lti) > 0){
        
        interpol_mui <- interpol.fct(FPCAobj$mu, FPCAobj$workGrid, Lti)
        interpol_FPCi <- apply(FPCAobj$phi, 2, function(x) return(interpol.fct(x, FPCAobj$workGrid, Lti)))
        interpolCov <- interpCovMat.fct(FPCAobj$fittedCov, FPCAobj$workGrid, Lti)
        
        scores[i,] <- t(matrix(rep(FPCAobj$lambda, length(Lti)), nrow = length(Lti), byrow = TRUE) * interpol_FPCi) %*%
          solve(interpolCov + FPCAobj$sigma2*diag(length(Lti))) %*% (Lyi - interpol_mui)
      }
    }
  } else if (type == 'face'){
    
    scores <- matrix(NA, nrow = length(marqueur), ncol = FPCAobj$npc)
    min_grilles.temps <- min(unlist(grilles.temps), na.rm = TRUE)
    max_grilles.temps <- max(unlist(grilles.temps), na.rm = TRUE)
    
    for (i in 1:length(marqueur)){
      
      Lti <- temps[[i]][temps[[i]]>=min_grilles.temps & temps[[i]]<max_grilles.temps]
      Lyi <- marqueur[[i]][temps[[i]]>=min_grilles.temps & temps[[i]]<max_grilles.temps]
      Lti <- Lti[!is.na(Lyi)]; Lyi <- Lyi[!is.na(Lyi)];
      
      if (length(Lti) > 0){
        
        interpol_mui <- interpol.fct(FPCAobj$mu.new, FPCAobj$argvals.new, Lti)
        interpol_FPCi <- apply(FPCAobj$eigenfunctions, 2, function(x) return(interpol.fct(x, FPCAobj$argvals.new, Lti)))
        interpolCov <- interpCovMat.fct(FPCAobj$Chat.new, FPCAobj$argvals.new, Lti)
        
        scores[i,] <- t(matrix(rep(FPCAobj$eigenvalues, length(Lti)), nrow = length(Lti), byrow = TRUE) * interpol_FPCi) %*%
          solve(interpolCov + FPCAobj$sigma2*diag(length(Lti))) %*% (Lyi - interpol_mui)
      }
    }
  }
  
  return(scores)
}

# Méthodes
lcmm.fct <- function(id, temps, marqueur, k = 2, 
                     regression = 'linéaire', lien = 'linéaire', noeuds.reg = 3,
                     noeuds.lien = 3, B.init = "modele0", nmes.min = 3) {
  
  # 1 - Données
  
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur, noeuds.reg = rep(noeuds.reg, length(id)))
  data <- data %>%
    filter(!is.na(marqueur)) %>%
    group_by(id) %>%
    filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    ungroup() %>%
    as.data.frame()
  
  reg <- if (lien == 'linéaire') 'linear' else paste(noeuds.lien, '-equi-splines', sep = "")
  
  
  # 2 - Modèle.s 
  
  modele0 <- modele1 <- NULL
  
  # 2.1 - Linéaire
  if (regression == 'linéaire') {
    modele0 <- tryCatch(
      lcmm(marqueur ~ temps, 
           random = ~ temps, 
           subject = 'id', 
           ng = 1, idiag = TRUE, 
           data = data, link = reg),
      error = function(e) { message("Erreur modele0 (linéaire) : ", e$message); return(NULL) })
    
    if (!is.null(modele0)) {
      B_init <- if (!is.null(B.init) && (B.init == "modele0" || modele0$conv == 1)) modele0 else B.init
      modele1 <- tryCatch(
        lcmm(marqueur ~ temps, 
             mixture = ~ temps, 
             random = ~ temps,
             subject = 'id', ng = k, 
             idiag = TRUE, data = data, 
             link = reg, B = B_init),
        error = function(e) { message("Erreur modele1 (linéaire) : ", e$message); return(NULL) })}}
  
  # 2.2 - Quadratique
  if (regression == 'quadratique') {
    modele0 <- tryCatch(
      lcmm(marqueur ~ temps + I(temps^2), 
           random = ~ temps + I(temps^2), 
           subject = 'id',
           ng = 1, idiag = TRUE, 
           data = data, link = reg),
      error = function(e) { message("Erreur modele0 (quadratique) : ", e$message); return(NULL) })
    
    if (!is.null(modele0)) {
      B_init <- if (!is.null(B.init) && (B.init == "modele0" || modele0$conv == 1)) modele0 else B.init
      modele1 <- tryCatch(
        lcmm(marqueur ~ temps + I(temps^2), 
             mixture = ~ temps + I(temps^2),
             random = ~ temps + I(temps^2), 
             subject = 'id', ng = k, 
             idiag = TRUE, data = data,
             link = reg, B = B_init),
        error = function(e) { message("Erreur modele1 (quadratique) : ", e$message); return(NULL) })}}
  
  # 2.3 - Splines
  if (regression == 'splines') {
    modele0 <- tryCatch(
      lcmm(marqueur ~ ns(temps, df = unique(noeuds.reg)), 
           random = ~ ns(temps, df = unique(noeuds.reg)),
           subject = 'id', 
           ng = 1, idiag = TRUE, 
           data = data, link = reg),
      error = function(e) { message("Erreur modele0 (splines) : ", e$message); return(NULL) })
    
    if (!is.null(modele0)) {
      B_init <- if (!is.null(B.init) && (B.init == "modele0" || modele0$conv == 1)) modele0 else B.init
      modele1 <- tryCatch(
        lcmm(marqueur ~ ns(temps, df = unique(noeuds.reg)), 
             mixture = ~ ns(temps, df = unique(noeuds.reg)),
             random = ~ ns(temps, df = unique(noeuds.reg)), 
             subject = 'id', ng = k, idiag = TRUE, 
             data = data, link = reg, B = B_init),
        error = function(e) { message("Erreur modele1 (splines) : ", e$message); return(NULL) })}}
  
  
  # 3 - Résultats
  
  if (is.null(modele1)) {
    res.summary <- NULL; res.prob <- NULL
    res.clusters <- rep(NA, length(unique(data$id)))
    conv.status <- NA} 
  else if (is.null(modele1$conv) || modele1$conv != 1) {
    res.summary <- NULL; res.prob <- NULL
    res.clusters <- rep(NA, length(unique(data$id)))
    conv.status <- modele1$conv} 
  else {
    res.summary <- summary(modele1)
    res.prob <- postprob(modele1)
    res.clusters <- as.numeric(factor(modele1$pprob[, 2]))
    conv.status <- modele1$conv}
  
  
  # 4 - Graphique
  
  data.plot <- merge(data, data.frame(id = unique(data$id), class = res.clusters),
                     by = "id", all.x = TRUE)
  res.graph <- ggplot(data.plot, aes(x = temps, y = marqueur, group = id, color = factor(class))) +
    geom_line(alpha = 0.5) + geom_point() +
    scale_color_viridis(discrete = TRUE, option = 'A', begin = 0.4, end = 0.7, alpha = 0.9) +
    labs(title = paste("LCMM/", regression, " - Trajectoires selon cluster"), x = "Temps",
         y = "Marqueur", color = "Clusters") +
    theme_minimal()
  
  return(list(modele = res.summary, probabilite = res.prob,
              clusters = res.clusters, convergence = conv.status, graphique = res.graph))
}
hlme.fct <- function(id, temps, marqueur, k = 2, 
                     regression = 'linéaire', noeuds.reg = 3, 
                     B.init = "modele0", nmes.min = 3) {
  
  # 1 - Données
  
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur, noeuds.reg = rep(noeuds.reg, length(id)))
  data <- data %>%
    filter(!is.na(marqueur)) %>%
    group_by(id) %>%
    filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    ungroup() %>%
    as.data.frame()
  
  
  # 2 - Modèle.s
  
  modele0 <- modele1 <- NULL
  
  # 2.1 - Linéaire
  if (regression == 'linéaire') {
    modele0 <- tryCatch(
      hlme(fixed = marqueur ~ temps, 
           random = ~ temps, 
           subject = 'id', 
           ng = 1, data = data),
      error = function(e) { message("Erreur modele0 (linéaire) : ", e$message); return(NULL) })
    
    if (!is.null(modele0)) {
      B_init <- if (!is.null(B.init) && (B.init == "modele0" | modele0$conv == 1)) modele0 else B.init
      modele1 <- tryCatch(
        hlme(fixed = marqueur ~ temps, 
             mixture = ~ temps, 
             random = ~ temps, 
             subject = 'id', ng = k, 
             data = data, B = B_init),
        error = function(e) { message("Erreur modele1 (linéaire) : ", e$message); return(NULL) })}}
  
  
  # 2.2 - Quadratique
  if (regression == 'quadratique') {
    modele0 <- tryCatch(
      hlme(fixed = marqueur ~ temps + I(temps^2), 
           random = ~ temps + I(temps^2), 
           subject = 'id', 
           ng = 1, data = data),
      error = function(e) { message("Erreur modele0 (quadratique) : ", e$message); return(NULL) })
    
    if (!is.null(modele0)) {
      B_init <- if (!is.null(B.init) && (B.init == "modele0" | modele0$conv == 1)) modele0 else B.init
      modele1 <- tryCatch(
        hlme(fixed = marqueur ~ temps + I(temps^2), 
             mixture = ~ temps + I(temps^2),
             random = ~ temps + I(temps^2), 
             subject = 'id', ng = k, 
             data = data, B = B_init),
        error = function(e) { message("Erreur modele1 (quadratique) : ", e$message); return(NULL) })}}
  
  
  # 2.3 - Splines
  if (regression == 'splines') {
    modele0 <- tryCatch(
      hlme(fixed = marqueur ~ ns(temps, df = unique(noeuds.reg)), 
           random = ~ ns(temps, df = unique(noeuds.reg)),
           subject = 'id', 
           ng = 1, data = data),
      error = function(e) { message("Erreur modele0 (splines) : ", e$message); return(NULL) })
    
    if (!is.null(modele0)) {
      B_init <- if (!is.null(B.init) && (B.init == "modele0" | modele0$conv == 1)) modele0 else B.init
      modele1 <- tryCatch(
        hlme(fixed = marqueur ~ ns(temps, df = unique(noeuds.reg)), 
             mixture = ~ ns(temps, df = unique(noeuds.reg)),
             random = ~ ns(temps, df = unique(noeuds.reg)), 
             subject = 'id', ng = k, 
             data = data, B = B_init),
        error = function(e) { message("Erreur modele1 (splines) : ", e$message); return(NULL) })}}
  
  
  # 3 - Résultats
  
  if (is.null(modele1)) {
    res.summary <- NULL; res.prob <- NULL
    res.clusters <- rep(NA, length(unique(data$id)))
    conv.status <- NA} 
  else if (is.null(modele1$conv) || modele1$conv != 1) {
    res.summary <- NULL; res.prob <- NULL
    res.clusters <- rep(NA, length(unique(data$id)))
    conv.status <- modele1$conv} 
  else {
    res.summary <- summary(modele1)
    res.prob <- postprob(modele1)
    res.clusters <- as.numeric(factor(modele1$pprob[, 2]))
    conv.status <- modele1$conv}
  
  
  # 4 - Graphique
  
  data.plot <- merge(data, data.frame(id = unique(data$id), class = res.clusters),
                     by = "id", all.x = TRUE)
  res.graph <- ggplot(data.plot, aes(x = temps, y = marqueur, group = id, color = factor(class))) +
    geom_line(alpha = 0.5) + geom_point() +
    scale_color_viridis(discrete = TRUE, option = 'A', begin = 0.4, end = 0.7, alpha = 0.9) +
    labs(title = paste("HLME/", regression, " - Trajectoires selon cluster"), 
         x = "Temps", y = "Marqueur", color = "Clusters") +
    theme_minimal()
  
  return(list(modele = res.summary, probabilite = res.prob,
              clusters = res.clusters, convergence = conv.status, graphique = res.graph))
}
acpfm.fct <- function(id, temps, marqueur, k = 2, k.mult = TRUE, 
                      clustering = 'km', pam.type = F, pve = 0.99, 
                      grid = 51, binnedData = "AUTO", nmes.min = 3){
  
  # 1 - ACPF
  
  # Données
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur) 
  data <- data %>%
    group_by(id) %>%
    filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    ungroup() 
  data.fdapace <- list(
    temps = split(data$temps, data$id),
    marqueur = split(data$marqueur, data$id))
  
  # ACPF
  res.acpf <- FPCA(data.fdapace$marqueur, 
                   data.fdapace$temps, 
                   list(dataType = 'Sparse', 
                        imputeScores = F,
                        FVEthreshold = pve,
                        nRegGrid = grid,
                        useBinnedData = binnedData))
  grid.tps <- res.acpf$workGrid
  res.scores <- scores.acpf(res.acpf, data.fdapace$marqueur, data.fdapace$temps, grid.tps)
  res.lambda <- res.acpf$lambda
  
  
  # 2 - Clustering 
  
  if (clustering == 'km'){
    res.clust <- kmeans(res.scores, centers = k)
    res.class <- res.clust$cluster
  }
  
  if (clustering == 'wkm'){
    poids <- sqrt(res.lambda)
    scores.pond <- sweep(res.scores, 2, poids, '*')
    res.clust <- kmeans(scores.pond, centers = k)
    res.class <- res.clust$cluster
  }
  
  if (clustering == 'pam'){
    res.clust <- pam(res.scores, k = k, pamonce = pam.type)
    res.class <- res.clust$cluster
  }
  
  if (clustering == 'gmm'){
    if (k.mult == T){res.clust <- Mclust(res.scores, G = c(1:k))}
    else if (k.mult == F){res.clust <- Mclust(res.scores, G = k)}
    res.class <- as.numeric(factor(res.clust$classification))
  }
  
  if (clustering == 'rf'){
    RFP <- randomForest(res.scores, proximity = T)
    dist.matrix <- as.dist(1 - RFP$proximity)
    res.clust <- pam(dist.matrix, k = k)
    res.class <- res.clust$cluster
  }
  
  
  # 3 - Plots
  
  # Clustering des courbes
  data.clust <- data.frame(class = res.class,
                           id = unique(data$id))
  data.plot <- merge(data, data.clust, by = "id")
  res.graph <- ggplot(data.plot, aes(x = temps, y = marqueur, 
                                     group = id, color = factor(class))) +
    geom_line(alpha = 0.5) +
    geom_point() +
    scale_color_viridis(discrete = T, option = 'A', begin = 0.35, end = 0.75, alpha = 0.9) +
    labs(title = paste("ACPF/",clustering," - Trajectoires colorées selon le cluster", sep = ""),  
         x = "Temps", y = "Marqueur", color = "Clusters") +
    theme_minimal()
  
  plots.scores <- lapply(1:ncol(res.scores), function(i) {
    ggplot(data.frame(scores = res.scores[, i]), aes(x = scores)) + 
      geom_density(color = "white", fill = "darkblue", alpha = 0.8) + 
      labs(x = paste("Scores Dim", i), y = "Fréquence") +
      theme_minimal()
  })
  
  
  return(list(res.acpf = res.acpf,
              res.clust = res.clust, 
              clusters = res.class, 
              scores = list(xiEst = res.scores, 
                            distributions = plots.scores),
              graphique = res.graph))
}
kml.fct <- function(id, temps, marqueur, k = 2, 
                    distance = 'euclidienne', nmes.min = 3){
  
  # 1 - Données 
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur)
  data.wide <- data %>%
    pivot_wider(names_from = temps, values_from = marqueur) %>% 
    arrange(id) %>% 
    filter(rowSums(!is.na(across(-id))) >= nmes.min) %>%
    dplyr::select(id, where(~ !all(is.na(.))))
  data.wide <- data.wide[, c('id', as.character(sort(as.numeric(names(data.wide)[-1]))))]

  # 2 - KML
  long.data <- cld(as.data.frame(data.wide))
  assign("long.data", long.data, envir = .GlobalEnv) # ???
  if (distance == 'frechet'){
    kml(long.data, nbClusters = k, parAlgo = parALGO(distance = FrechetDistance(x,y)))} 
  else if (distance == 'dtw'){
    kml(long.data, nbClusters = k, parAlgo = parALGO(distance = DTWDistance(x,y)))}
  else {
    kml(long.data, nbClusters = k)}
  res.class <- as.numeric(getClusters(long.data, k))
  
  # 3 - Plot
  data.clust <- data.frame(class = res.class,
                           id = unique(data.wide$id))
  data.plot <- merge(data, data.clust, by = "id")
  res.graph <- ggplot(data.plot, aes(x = temps, y = marqueur, 
                                     group = id, color = factor(class))) +
                geom_line(alpha = 0.5) + 
                geom_point() +
                scale_color_viridis(discrete = T, option = 'A', begin = 0.4, end = 0.7, alpha = 0.9) +
                labs(title = paste("KML/",distance," - Trajectoires colorées selon le cluster", sep = ""),  
                     x = "Temps", y = "Marqueur", color = "Clusters") +
                theme_minimal()
  
  return(list(clusters = res.class, 
              graphique = res.graph))
}
dtw.fct <- function(id, temps, marqueur, k = 2, 
                    type = 'partitional', centroid = 'km', 
                    distance = 'dtw_basic', nmes.min = 3){
  
  # 1 - Données
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur)
  data <- data %>%
    filter(!is.na(marqueur)) %>%
    group_by(id) %>%                        
    filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    ungroup() 
  data.dtw <- list(
    temps = split(data$temps, data$id),
    marqueur = split(data$marqueur, data$id))
  
  # 2 - DTW
  res.clust <- tsclust(data.dtw$marqueur, k = k, 
                       type = type, 
                       centroid = centroid, 
                       distance = distance) 
  res.class <- res.clust@cluster
  
  # 3 - Plot
  data.clust <- data.frame(class = res.class,
                           id = unique(data$id))
  data.plot <- merge(data, data.clust)
  res.graph <- ggplot(data.plot, aes(x = temps, y = marqueur, 
                                     group = id, color = factor(class))) +
                geom_line(alpha = 0.5) +  
                geom_point() +
                scale_color_viridis(discrete = T, option = 'A', begin = 0.4, end = 0.7, alpha = 0.9) +
                labs(title = paste("DTW/",centroid,"/",distance," - Trajectoires colorées selon le cluster", sep = ""),
                     x = "Temps", y = "Marqueur", color = "Clusters") +
                theme_minimal() 

  
  return(list(res.clust = res.clust,
              clusters = res.class,
              graphique = res.graph))
}
matClust.fct <- function(id, temps, marqueur, k = 2, 
                        clustering = 'km', 
                        dist.method = "frechet", nmes.min = 3){
  
  dist.method <- match.arg(dist.method)
  
  # 1 - Données
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur) 
  data <- data %>%
    dplyr::filter(!is.na(marqueur)) %>%
    dplyr::group_by(id) %>%                        
    dplyr::filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    dplyr::ungroup()
  
  data.split <- list(
    temps = split(data$temps, data$id),
    marqueur = split(data$marqueur, data$id))
  
  # 2 - Distances
  if (dist.method == "frechet") {
    mat.dist <- mat.frechet(data.split$temps, data.split$marqueur, ts = 0.1) # Fréchet
  } else if (dist.method == "dtw") {
    mat.dist <- mat.frechet(data.split$temps, data.split$marqueur, ts = 1e-6) # DTW
  }
  
  # 3 - Clustering 
  if (clustering == 'km'){
    res.clust <- kmeans(mat.dist, centers = k)
    res.class <- res.clust$cluster
  } else if (clustering == 'pam'){
    res.clust <- pam(mat.dist, k = k)
    res.class <- res.clust$cluster
  } else {
    stop("Méthode de clustering non reconnue. Choisir 'km' ou 'pam'.")
  }
  
  # 4 - Graphiques
  id.valides <- as.numeric(names(data.split$temps))
  data.clust <- data.frame(class = res.class,
                           id = id.valides)
  data.plot <- merge(data, data.clust, by = "id")
  
  titre <- paste0( toupper(dist.method), "/", clustering,
                   " - Trajectoires colorées selon le cluster")
  
  res.graph <- ggplot(data.plot, aes(x = temps, y = marqueur, 
                                     group = id, color = factor(class))) +
    geom_line(alpha = 0.5) +  
    geom_point() +
    scale_color_viridis(discrete = TRUE, option = 'A', begin = 0.4, end = 0.7, alpha = 0.9) +  
    labs(title = titre, x = "Temps", y = "Marqueur", color = "Clusters") +
    theme_minimal()
  
  return(list(
    clusters = res.class, 
    graphique = res.graph,
    distance = mat.dist
  ))
}

# Distances Manuelles
mat.dist <- function(temps, marqueur, ts = 0.1){
  
  distances <- matrix(0, nrow = length(temps), 
                         ncol = length(temps)) 
  
  for (i in 1:(length(temps)-1)){
    for (j in (i+1):length(temps)){
      silent.res <- capture.output(distances[i,j] <- distances[j,i] <- distFrechet(temps[[i]], marqueur[[i]],
                                                                                   temps[[j]], marqueur[[j]], 
                                                                                   FrechetSumOrMax = 'sum', 
                                                                                   timeScale = ts))}}
  
  return(frechet)
}
dba.fct <- function(id, temps, marqueur, res.clust, data = data.frame(), 
                    data.list = data.frame(), data.wide = data.frame(), nmes.min = 3){
  
  # 1 - Données
  
  if (length(data) == 0){
    id <- as.numeric(id)
    data <- data.frame(id, temps, marqueur) 
    data <- data %>%
      filter(!is.na(marqueur)) %>%
      group_by(id) %>%
      filter(sum(!is.na(marqueur)) >= nmes.min) %>%
      ungroup()}
  
  if (length(data.list) == 0){
    data.list <- list(temps = split(data$temps, data$id),
                      marqueur = split(data$marqueur, data$id))}
  
  if (length(data.wide) == 0){
    data.wide <- data %>% 
      pivot_wider(names_from = temps, values_from = marqueur) %>% 
      dplyr::select(-id)
    data.wide <- data.wide[, as.character(sort(as.numeric(names(data.wide))))]
    colnames(data.wide) <- round(as.numeric(colnames(data.wide)), 7)}

  last.tps <- data %>% 
    group_by(id) %>% 
    filter(temps == max(temps)) %>% 
    ungroup() %>% 
    pull(temps)
  
  # 2 - Calcul des centroids 
  
  k <- length(unique(res.clust))
  cent.y <- vector('list', k)
  cent.t <- vector('list', k)
  cent.id <- vector('list', k)
  cols <- sort(as.numeric(names(data.wide)))
  centroids <- as.data.frame(matrix(NA, k, length(cols), dimnames = list(c(), cols)))
  
  for (i in 1:k){
    id.clust <- which(res.clust == i)
    if (length(id.clust) > 1){
      cent.y[[i]] <- dba(data.list$marqueur[id.clust],
                         centroid = data.list$marqueur[id.clust][[which(rowSums(!is.na(data.wide[id.clust,])) == median(rowSums(!is.na(data.wide[id.clust,]))))[1]]])
      avail.t <- sort(unique(unlist(data.list$temps[id.clust])))
      avail.t <- avail.t[avail.t <= median(last.tps[id.clust])]
      cent.t[[i]] <- quantile(avail.t, probs = seq(0, 1, length.out = length(cent.y[[i]])), names = F, type = 3)} 
    else {
      cent.y[[i]] <- as.numeric(data.wide[id.clust, which(!is.na(data.wide[id.clust,]))])
      cent.t[[i]] <- as.numeric(colnames(data.wide[,which(!is.na(data.wide[id.clust,]))]))
    }
    
    keep <- !duplicated(round(cent.t[[i]], 7))
    cent.t[[i]] <- cent.t[[i]][keep]
    cent.y[[i]] <- cent.y[[i]][keep]
    cent.id[[i]] <- rep(i, length(cent.y[[i]]))
    tps.clust <- as.character(round(cent.t[[i]], 7))
    centroids[i, tps.clust] <- cent.y[[i]]
  }
  
  # 3 - Plots
  
  data.plot <- merge(data, 
                     data.frame(id = unique(data$id), class = res.clust), 
                     by = 'id')
  
  data.cent <- data.frame(centroid = unlist(cent.id),
                          t = unlist(cent.t),
                          y = unlist(cent.y))

  res.plot <- ggplot(data.plot, aes(x = temps, y = marqueur, color = factor(class))) +
                geom_line(aes(group = id), alpha = 0.5) +  
                geom_line(data = data.cent, aes(x = t, y = y, group = factor(centroid)), color = 'black', linewidth = 0.8) + 
                labs(x = "Temps", y = "Marqueur", color = "Clusters") +
                theme_minimal()
  
  return(list(centroids = centroids, 
              graphique = res.plot))
}

# Évaluation
sil.fct <- function(id, temps, marqueur, res.clust, 
                    matrice.dist = matrix(), distance = 'DTW', nmes.min = 3){
  
  if (length(matrice.dist) == 1){
    
    if (distance == 'fréchet'){
      
      # 1 - Données
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur) 
      data <- data %>%
        filter(!is.na(marqueur)) %>%
        group_by(id) %>%                        
        filter(sum(!is.na(marqueur)) >= nmes.min) %>%
        ungroup() 
      data.list <- list(
        temps = split(data$temps, data$id),
        marqueur = split(data$marqueur, data$id))
      
      # 2 - Fréchet
      matrice.dist <- mat.frechet(data.list$temps, data.list$marqueur) 
      
    } else if (distance == 'DTW') {
      
      # 1 - Données 
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur)
      data.wide <- data %>%
        pivot_wider(names_from = temps, values_from = marqueur) %>% 
        arrange(id) %>% 
        filter(rowSums(!is.na(across(-id))) >= nmes.min) %>%
        dplyr::select(where(~ !all(is.na(.))))
      data.wide <- data.wide[, as.character(sort(as.numeric(names(data.wide)[-1])))]
      
      # 2 - Matrice de distance (DTW, DTW2)
      data.list <- split(data.wide, f = seq_len(nrow(data.wide)))
      matrice.dist <- dist(data.wide, method = distance)
      matrice.dist <- as.matrix(matrice.dist)
      
    } else if (distance == 'DTWB'){
      
      # 1 - Données
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur) 
      data <- data %>%
        group_by(id) %>%                        
        filter(sum(!is.na(marqueur)) >= nmes.min) %>%
        ungroup() 
      data.list <- list(
        temps = split(data$temps, data$id),
        marqueur = split(data$marqueur, data$id))
      
      # 2 - Matrice de distance 
      matrice.dist <- mat.dist(data.list$temps, data.list$marqueur, ts = 10^(-6)) 
      
    }
  }
  # 3 - Indice de Silhouette
  if (any(is.na(res.clust))){return(NA)}
  ind.silhouette <- silhouette(res.clust, matrice.dist)
  res.silhouette <- tryCatch(
    mean(ind.silhouette[,3]),
    error = function(e) {
      message("Erreur dans le calcul de l'indice de Silhouette, les individus appartiennent tous au même groupe. ")
      return(NA)})
  
  return(res.silhouette)
}
dbi.fct <- function(id, temps, marqueur, res.clust, 
                    distance = 'DTW', nmes.min = 3){
  
  # 1 - Données
  
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur) 
  data <- data %>%
    filter(!is.na(marqueur)) %>%
    group_by(id) %>%
    filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    ungroup()
  
  data.wide <- data %>% 
    pivot_wider(names_from = temps, values_from = marqueur) %>% 
    dplyr::select(-id)
  data.wide <- data.wide[, as.character(sort(as.numeric(names(data.wide))))]
  colnames(data.wide) <- round(as.numeric(colnames(data.wide)), 7)
  
  data.list <- list(temps = split(data$temps, data$id),
                    marqueur = split(data$marqueur, data$id))
  
  
  # 2 - Indice de Davies-Bouldin
  
  if (any(is.na(res.clust)) | length(unique(res.clust)) == 1)
    {return(list(score = NA, centroids = NA))}
  
  # Centroids
  k <- length(unique(res.clust))
  res.dba <- dba.fct(res.clust = res.clust,
                     data = data,
                     data.list = data.list, 
                     data.wide = data.wide)
  centroids <- res.dba$centroids

  # Distances inter.
  if (distance == 'DTW'){
    cent.list <- split(centroids, f = 1:k)
    dist.inter <- as.matrix(dist(cent.list, method = distance))
  } 
  else if (distance == 'fréchet'){
    cent.list <- list(temps = vector('list', k),
                      marqueur = vector('list', k))
    for (i in 1:k){
      cent.list$temps[[i]] <- as.numeric(colnames(centroids[i,which(!is.na(centroids[i,]))]))
      cent.list$marqueur[[i]] <- as.numeric(centroids[i,which(!is.na(centroids[i,]))])
    }
    dist.inter <- mat.frechet(cent.list$temps, cent.list$marqueur)
  } 
  else if (distance == 'DTWB'){
    cent.list <- list(temps = vector('list', k),
                      marqueur = vector('list', k))
    for (i in 1:k){
      cent.list$temps[[i]] <- as.numeric(colnames(centroids[i,which(!is.na(centroids[i,]))]))
      cent.list$marqueur[[i]] <- as.numeric(centroids[i,which(!is.na(centroids[i,]))])
    }
    dist.inter <- mat.frechet(cent.list$temps, cent.list$marqueur, ts = 10^(-6))}
  
  # Distance intra.
  var.intra <- numeric(k)
  for (i in 1:k){
    id.clust <- which(res.clust == unique(res.clust)[i])
    dist.intra <- numeric(length(id.clust))
    for (j in 1:length(id.clust)){
      if (distance == "DTW"){
        dist.intra[j] <- dtw(data.wide[id.clust[j],], centroids[i,])$distance
      }
      else if (distance == 'fréchet'){
        silent.res <- capture.output(dist.intra[j] <- distFrechet(data.list$temps[[id.clust[j]]], 
                                                                  data.list$marqueur[[id.clust[j]]],
                                                                  cent.list$temps[[i]], 
                                                                  cent.list$marqueur[[i]], 
                                                                  FrechetSumOrMax = 'sum'))
      }
      else if (distance == 'DTWB'){
        silent.res <- capture.output(dist.intra[j] <- distFrechet(data.list$temps[[id.clust[j]]], 
                                                                  data.list$marqueur[[id.clust[j]]],
                                                                  cent.list$temps[[i]], 
                                                                  cent.list$marqueur[[i]], 
                                                                  FrechetSumOrMax = 'sum',
                                                                  timeScale = 10^(-6)))}}
    var.intra[i] <- mean(dist.intra, na.rm = T)
  }
  
  dbi <- sapply(1:k, function(i) {
    max(sapply(setdiff(1:k, i), function(j) (var.intra[i] + var.intra[j]) / (dist.inter[i, j])), na.rm = T)})
  
  return(list(score = mean(dbi, na.rm = T),
              centroids = res.dba$graphique))
}
chi.fct <- function(id, temps, marqueur, res.clust, 
                    matrice.dist = matrix(), distance = 'DTW', nmes.min = 3){
  
  if (length(matrice.dist) == 1){
    
    if (distance == 'fréchet'){
      
      # 1 - Données
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur) 
      data <- data %>%
        filter(!is.na(marqueur)) %>%
        group_by(id) %>%                        
        filter(sum(!is.na(marqueur)) >= nmes.min) %>%
        ungroup() 
      data.frechet <- list(
        temps = split(data$temps, data$id),
        marqueur = split(data$marqueur, data$id))
      
      # 2 - Fréchet
      matrice.dist <- mat.frechet(data.frechet$temps, data.frechet$marqueur) 
      
    } else if (distance == 'DTW') {
      
      # 1 - Données 
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur)
      data.wide <- data %>%
        pivot_wider(names_from = temps, values_from = marqueur) %>% 
        arrange(id) %>% 
        filter(rowSums(!is.na(across(-id))) >= nmes.min) %>%
        dplyr::select(where(~ !all(is.na(.))))
      data.wide <- data.wide[, as.character(sort(as.numeric(names(data.wide)[-1])))]
      
      # 2 - Matrice de distance (DTW, DTW2)
      data.list <- split(data.wide, f = seq_len(nrow(data.wide)))
      matrice.dist <- dist(data.wide, method = distance)
      matrice.dist <- as.matrix(matrice.dist)
      
    } else if (distance == 'DTWB'){
      
      # 1 - Données
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur) 
      data <- data %>%
        group_by(id) %>%                        
        filter(sum(!is.na(marqueur)) >= nmes.min) %>%
        ungroup() 
      data.frechet <- list(
        temps = split(data$temps, data$id),
        marqueur = split(data$marqueur, data$id))
      
      # 2 - Matrice de distance 
      matrice.dist <- mat.frechet(data.frechet$temps, data.frechet$marqueur, ts = 10^(-6)) 
      
    }
  }
  
  # 3 - Indice de Calinski-Harabasz
  if (any(is.na(res.clust)) | length(unique(res.clust)) == 1){return(NA)}
  CH <- cluster.stats(matrice.dist, res.clust)$ch
  
  return(CH)
}
indices.fct <- function(id, temps, marqueur, methodes, 
                        matrice.dist = matrix(), distance = 'DTW', nmes.min = 3, 
                        indices = c('Silhouette','Davies-Bouldin','Calinski-Harabasz')){
  
  if (length(matrice.dist) == 1){
    
    if (distance == 'fréchet'){
      
      # 1 - Données
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur) 
      data <- data %>%
        filter(!is.na(marqueur)) %>%
        group_by(id) %>%                        
        filter(sum(!is.na(marqueur)) >= nmes.min) %>%
        ungroup() 
      data.frechet <- list(
        temps = split(data$temps, data$id),
        marqueur = split(data$marqueur, data$id))
      
      # 2 - Matrice de distance (Fréchet)
      matrice.dist <- mat.frechet(data.frechet$temps, data.frechet$marqueur) 
      
    } else if (distance == 'DTW') {
      
      # 1 - Données 
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur)
      data.wide <- data %>%
        pivot_wider(names_from = temps, values_from = marqueur) %>% 
        arrange(id) %>% 
        filter(rowSums(!is.na(across(-id))) >= nmes.min) %>%
        dplyr::select(where(~ !all(is.na(.))))
      data.wide <- data.wide[, as.character(sort(as.numeric(names(data.wide)[-1])))]
      
      # 2 - Matrice de distance (DTW, DTW2)
      data.list <- split(data.wide, f = seq_len(nrow(data.wide)))
      matrice.dist <- dist(data.wide, method = distance)
      matrice.dist <- as.matrix(matrice.dist)
      
    } else if (distance == 'DTWB'){
      
      # 1 - Données
      id <- as.numeric(id)
      data <- data.frame(id, temps, marqueur) 
      data <- data %>%
        group_by(id) %>%                        
        filter(sum(!is.na(marqueur)) >= nmes.min) %>%
        ungroup() 
      data.frechet <- list(
        temps = split(data$temps, data$id),
        marqueur = split(data$marqueur, data$id))
      
      # 2 - Matrice de distance 
      matrice.dist <- mat.frechet(data.frechet$temps, data.frechet$marqueur, ts = 10^(-6))
      
    }
  }
  
  # 3 - Évaluation des clusters
  m <- length(methodes)
  scores <- list(si = numeric(m),
                 dbi = numeric(m),
                 chi  = numeric(m))
  for (i in 1:m){
    if ('Silhouette' %in% indices){
      scores$si[i] <- sil.fct(res.clust = methodes[[i]], 
                              matrice.dist = matrice.dist, 
                              distance = distance)}
    if ('Davies-Bouldin' %in% indices){
      scores$dbi[i] <- dbi.fct(id, temps, marqueur, methodes[[i]], distance)$score}
    if ('Calinski-Harabasz' %in% indices){
      scores$chi[i]  <- chi.fct(res.clust = methodes[[i]], 
                                matrice.dist = matrice.dist, 
                                distance = distance)}}

  noms <- c()
  for (i in 1:length(methodes)){noms <- c(noms, names(methodes)[i])}
  plots.scores <- lapply(1:length(scores), function(i) {
    data <- data.frame(methode = noms, score = scores[[i]])
    if (i == 1){
      y.limits <- c(-1,1)
      indice.nom <- 'Silhouette'}
    else if (i == 2){
      indice.nom <- 'Davies-Bouldin'
      if (distance == 'frechet'){
        y.limits <- c(0,30)}
      else {
        y.limits <- c(0,20)}}
    else if (i == 3){
      indice.nom <- 'Calinski-Harabasz'
      if (distance == 'fréchet')
      {y.limits <- c(0,1000)}
      else
      {y.limits <- c(0,100)}}
    ggplot(data, aes(x = methode, y = score)) +
      geom_bar(stat = "identity", fill = "steelblue", color = "steelblue", width = 0.6) +  
      geom_text(aes(label = round(score,2)), vjust = -1) +
      theme_minimal() +  
      labs(title = paste("Indice de", indice.nom), x = '', y = '') + 
      coord_cartesian(ylim = y.limits) +
      theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 8))})
  
  return(list(scores = scores, 
              graphique = plots.scores))
}
eval.fct <- function(cluster.ref, clusters.test){
  
  if (is.list(clusters.test)){
    
    L <- length(clusters.test)
    resultats <- numeric(L)
    
    for (i in 1:L){
      if (any(is.na(clusters.test[[i]]))){resultats[i] <- NA}
      else {resultats[i] <- adjustedRandIndex(cluster.ref, clusters.test[[i]])}
    }
    
    noms <- c()
    for (i in 1:length(clusters.test)){noms <- c(noms, names(clusters.test)[i])}
    data <- data.frame(methodes = noms, scores = resultats)
    graphique <- ggplot(data, aes(x = methodes, y = scores)) +
      geom_bar(stat = "identity", fill = "steelblue", color = "steelblue", width = 0.6) +  
      geom_text(aes(label = round(scores, 2)), vjust = -1) +
      theme_minimal() +  
      labs(title = paste("Indice ARI"), x = '', y = '') + 
      coord_cartesian(ylim = c(0,1.2)) +
      theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 8))
  } 
  
  else {
    if (any(is.na(clusters.test))){resultats <- NA}
    else {resultats <- adjustedRandIndex(cluster.ref, clusters.test)}
    graphique <- NULL
  }
  
  return(list(scores = resultats, 
              graphique = graphique))
}

# Comparaison
sim.fct <- function(class1, class2){
  
  if (any(is.na(class1)) | any(is.na(class2))){return(NA)}
  
  # 1 - Matrice de similarité pour chaque méthode
  
  mat.sim <- list()
  
  k <- 1
  for (class in list(class1, class2)){
    
    mat <- matrix(0, ncol = length(class), nrow = length(class))
    
    for (i in 1:(length(class)-1)){
      for (j in (i+1):length(class)){ 
        if (class[[i]] == class[[j]]){mat[i,j] <- 1}
      }
    }
    
    mat.sim[[k]] <- mat
    k <- k+1
  }
  
  # 2 - Calcul du score de similarité 
  
  sim.score <- 1 - sum((mat.sim[[1]]-mat.sim[[2]])^2) / (length(class)*(length(class)-1)/2)
  
  return(sim.score)
}
comp.fct <- function(methodes = list(), type = 'sim', color = 1){
  
  # 1 - Création de la matrice des scores (avec la diagonale)
  scoresSIM <- diag(1, ncol = length(methodes), nrow = length(methodes))
  scoresSIM[upper.tri(scoresSIM)] <- NA
  scoresARI <- diag(1, ncol = length(methodes), nrow = length(methodes))
  scoresARI[upper.tri(scoresARI)] <- NA
  for (i in 2:length(methodes)){
    for (j in 1:(i-1)){
      scoresSIM[i,j] <- sim.fct(methodes[[i]], methodes[[j]])
      scoresARI[i,j] <- eval.fct(methodes[[i]], methodes[[j]])$scores
    }
  }
  
  # 2 - Noms des méthodes
  noms <- names(methodes)
  rownames(scoresSIM) <- colnames(scoresSIM) <- noms
  rownames(scoresARI) <- colnames(scoresARI) <- noms
  
  # 3 - Palette de couleurs
  if (color == 1){colors <- viridis(50, option = "G", begin = 1, end = 0.5)}
  else if (color == 2){colors <- viridis(50, option = "G", begin = 0.5, end = 1)}
  else if (color == 3){colors <- viridis(50, begin = 1, end = 0.5, option = "A")}
  
  # 4 - Création de la heatmap
  graphiqueSIM <- Heatmap(scoresSIM,
                          name = " ",
                          column_title = "Matrice de similarité entre les méthodes (SIM)",
                          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          col = colors,
                          cluster_rows = F,
                          cluster_columns = F,
                          row_names_side = "left", 
                          row_names_gp = gpar(fontsize = 7),
                          column_names_side = "bottom", 
                          column_names_gp = gpar(fontsize = 7),
                          column_names_rot = 45,
                          na_col = 'white',
                          cell_fun = function(j, i, x, y, width, height, fill) { 
                            if (!is.na(scoresSIM[i, j])) { 
                              grid.text(round(scoresSIM[i, j], 2), x, y, gp = gpar(fontsize = 6, col = "black")) 
                            }
                          })
  graphiqueARI <- Heatmap(scoresARI,
                          name = " ",
                          column_title = "Matrice de similarité entre les méthodes (ARI)",
                          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          col = colors,
                          cluster_rows = F,
                          cluster_columns = F,
                          row_names_side = "left", 
                          row_names_gp = gpar(fontsize = 7),
                          column_names_side = "bottom", 
                          column_names_gp = gpar(fontsize = 7),
                          column_names_rot = 45,
                          na_col = 'white',
                          cell_fun = function(j, i, x, y, width, height, fill) { 
                            if (!is.na(scoresARI[i, j])) { 
                              grid.text(round(scoresARI[i, j], 2), x, y, gp = gpar(fontsize = 6, col = "black")) 
                            }
                          })
  
  return(list(scores = list(SIM = scoresSIM,
                            ARI = scoresARI), 
              graphique = list(SIM = graphiqueSIM,
                               ARI = graphiqueARI)))
}


