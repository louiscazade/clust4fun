
# Chargement des packages 
packages <- c( "ggplot2", "dplyr", "tidyr", "lcmm", "splines", 
               "fdapace", "face", "cluster", "randomForest", 
               "mclust", "kml", "dtwclust", "dtw", "fpc", "TSdist",
               "future", "future.apply", "parallel", 
               "viridis", "ComplexHeatmap", "grid")

invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)}
  library(pkg, character.only = TRUE)}))

# Fonctions ACPF manuelle
interpol.fct     <- function(y, t, tNew){
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
scores.acpf      <- function(FPCAobj, marqueur, temps, grilles.temps, type = 'fdapace'){
  
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

# Fonctions utilitaires
filter.mes <- function(df, id, valeurs, nmes.min) {
  df %>% 
    group_by({{id}}) %>% 
    filter(sum(!is.na({{valeurs}})) >= nmes.min) %>% 
    ungroup()}
plot.clust <- function(data, clusters, titre) {
  data.plot <- merge(data, data.frame(id = unique(data$id), class = clusters), by = "id", all.x = TRUE)
  ggplot(data.plot, aes(x = temps, y = marqueur, group = id, color = factor(class))) +
    geom_line(alpha = 0.5) + geom_point() +
    scale_color_viridis(discrete = TRUE, option = 'A', begin = 0.4, end = 0.7, alpha = 0.9) +
    labs(title = titre, x = "Temps", y = "Marqueur", color = "Clusters") +
    theme_minimal()}

# Méthodes de clustering
lcmm.fct <- function(id, temps, marqueur, k = 2, 
                     regression = 'linéaire', noeuds.reg = 3,
                     lien = 'linéaire', noeuds.lien = 3, 
                     B.init = "modele0", nmes.min = 3) {
  
  # 1 - Données
  id <- as.integer(factor(id))
  data <- data.frame(id, temps, marqueur)
  data <- filter.mes(data, id, marqueur, nmes.min)
  
  # 2 - Choix de la régression
  formules <- switch(regression,
                     'linéaire' = list(fixed = marqueur ~ temps, 
                                       random = ~ temps, 
                                       mixture = ~ temps),
                     'quadratique' = list(fixed = marqueur ~ temps + I(temps^2), 
                                          random = ~ temps + I(temps^2), 
                                          mixture = ~ temps + I(temps^2)),
                     'splines' = list(fixed = as.formula(paste0("marqueur ~ ns(temps, df = ", noeuds.reg, ")")),
                                      random = as.formula(paste0("~ ns(temps, df = ", noeuds.reg, ")")),
                                      mixture = as.formula(paste0("~ ns(temps, df = ", noeuds.reg, ")"))),
                     stop("Régression non reconnue"))
  
  # 3 - Choix du lien 
  link <- if (lien == 'linéaire') 'linear' else paste0(noeuds.lien, "-equi-splines")
  
  # 4 - Modèles
  modele0 <- tryCatch(
    lcmm(formules$fixed, random = formules$random, subject = 'id', ng = 1, data = data, link = link),
    error = function(e) { message("Erreur modele0 : ", e$message); return(NULL) })
  
  B_init <- if (!is.null(B.init) && (B.init == "modele0" || modele0$conv == 1)) modele0 else B.init
  
  modele1 <- tryCatch(
    lcmm(formules$fixed, mixture = formules$mixture, random = formules$random, subject = 'id', ng = k, data = data, link = link, B = B_init),
    error = function(e) { message("Erreur modele1 : ", e$message); return (NULL) })
  
  # 5 - Résultats
  if (is.null(modele1) || is.null(modele1$conv) || modele1$conv != 1) {
    res.summary <- NULL
    res.prob <- NULL
    res.clusters <- rep(NA, length(unique(data$id)))
    conv.status <- if(!is.null(modele1)) modele1$conv else NA
  } else {
    res.summary <- summary(modele1)
    res.prob <- postprob(modele1)
    res.clusters <- as.numeric(factor(modele1$pprob[,2]))
    conv.status <- modele1$conv}
  
  # 6 - Graphique 
  res.graph <- plot.clust(data, res.clusters, paste("LCMM/", regression, " - Trajectoires selon cluster"))
  
  return(list(modele = res.summary, probabilite = res.prob,
              clusters = res.clusters, convergence = conv.status, graphique = res.graph))
}
hlme.fct <- function(id, temps, marqueur, k = 2, 
                     regression = 'linéaire', noeuds.reg = 3, 
                     B.init = "modele0", nmes.min = 3) {
  
  # 1 - Données
  id <- as.integer(factor(id))
  data <- data.frame(id, temps, marqueur)
  data <- filter.mes(data.frame(id, temps, marqueur), id, marqueur, nmes.min)
  
  # 2 - Choix de la régression
  formules <- switch(regression,
                     'linéaire' = list(fixed = marqueur ~ temps, 
                                       random = ~ temps, 
                                       mixture = ~ temps),
                     'quadratique' = list(fixed = marqueur ~ temps + I(temps^2), 
                                          random = ~ temps + I(temps^2), 
                                          mixture = ~ temps + I(temps^2)),
                     'splines' = list(fixed = as.formula(paste0("marqueur ~ ns(temps, df = ", noeuds.reg, ")")),
                                      random = as.formula(paste0("~ ns(temps, df = ", noeuds.reg, ")")),
                                      mixture = as.formula(paste0("~ ns(temps, df = ", noeuds.reg, ")"))),
                     stop("Régression non reconnue"))
  
  # 3 - Modèles
  modele0 <- tryCatch(
    hlme(fixed = formules$fixed, random = formules$random, subject = 'id', ng = 1, data = data),
    error = function(e) { message("Erreur modele0 : ", e$message); return(NULL) })
  
  B_init <- if (!is.null(B.init) && (B.init == "modele0" || modele0$conv == 1)) modele0 else B.init
  
  modele1 <- tryCatch(
    hlme(fixed = formules$fixed, mixture = formules$mixture, random = formules$random, subject = 'id', ng = k, data = data, B = B_init),
    error = function(e) { message("Erreur modele1 : ", e$message); return(NULL) })
  
  # 4 - Résultats
  if (is.null(modele1) || is.null(modele1$conv) || modele1$conv != 1) {
    res.summary <- NULL
    res.prob <- NULL
    res.clusters <- rep(NA, length(unique(data$id)))
    conv.status <- if(!is.null(modele1)) modele1$conv else NA
  } else {
    res.summary <- summary(modele1)
    res.prob <- postprob(modele1)
    res.clusters <- as.numeric(factor(modele1$pprob[,2]))
    conv.status <- modele1$conv}
  
  # 5 - Graphique
  res.graph <- plot.clust(data, res.clusters, paste("HLME/", regression, " - Trajectoires selon cluster"))
  
  return(list(modele = res.summary, probabilite = res.prob,
              clusters = res.clusters, convergence = conv.status, graphique = res.graph))
}
acpf.fct <- function(id, temps, marqueur, k = 2, k.mult = TRUE, 
                     clustering = 'km', pam.type = FALSE, pve = 0.99, 
                     grid = 51, binnedData = "AUTO", nmes.min = 3) {
  
  # 1 - Données
  id <- as.integer(factor(id))
  data <- filter.mes(data.frame(id, temps, marqueur), id, marqueur, nmes.min)
  data.fdapace <- list(
    temps = split(data$temps, data$id),
    marqueur = split(data$marqueur, data$id))
  
  # 2 - FPCA
  res.acpf <- FPCA(data.fdapace$marqueur, data.fdapace$temps, 
                   list(dataType = 'Sparse', imputeScores = FALSE, 
                        FVEthreshold = pve, nRegGrid = grid, 
                        useBinnedData = binnedData))
  grid.tps <- res.acpf$workGrid
  res.scores <- scores.acpf(res.acpf, data.fdapace$marqueur, data.fdapace$temps, grid.tps)
  res.lambda <- res.acpf$lambda
  
  # 3 - Clustering
  if (clustering == 'km') {
    res.clust <- kmeans(res.scores, centers = k)
    res.class <- res.clust$cluster
  } else if (clustering == 'wkm') {
    scores.pond <- sweep(res.scores, 2, sqrt(res.lambda), '*')
    res.clust <- kmeans(scores.pond, centers = k)
    res.class <- res.clust$cluster
  } else if (clustering == 'pam') {
    res.clust <- pam(res.scores, k = k, pamonce = pam.type)
    res.class <- res.clust$cluster
  } else if (clustering == 'gmm') {
    G <- if (k.mult) 1:k else k
    res.clust <- Mclust(res.scores, G = G)
    res.class <- as.numeric(factor(res.clust$classification))
  } else if (clustering == 'rf') {
    RFP <- randomForest(res.scores, proximity = TRUE)
    dist.matrix <- as.dist(1 - RFP$proximity)
    res.clust <- pam(dist.matrix, k = k)
    res.class <- res.clust$cluster
  } else stop("Méthode de clustering non reconnue")
  
  # 4 - Graphiques
  res.graph <- plot.clust(data, res.class, paste("ACPF/", clustering, " - Trajectoires selon cluster"))
  
  plots.scores <- lapply(1:ncol(res.scores), function(i) {
    ggplot(data.frame(scores = res.scores[, i]), aes(x = scores)) +
      geom_density(color = "white", fill = "darkblue", alpha = 0.8) +
      labs(x = paste("Scores Dim", i), y = "Fréquence") +
      theme_minimal()})
  
  return(list(res.acpf = res.acpf,
              res.clust = res.clust,
              clusters = res.class,
              scores = list(xiEst = res.scores, distributions = plots.scores),
              graphique = res.graph))
}
kml.fct <- function(id, temps, marqueur, k = 2, 
                    distance = 'euclidienne', nmes.min = 3) {
  
  # 1 - Données
  id <- as.integer(factor(id))
  data <- filter.mes(data.frame(id, temps, marqueur), id, marqueur, nmes.min)
  data.wide <- data %>% 
    pivot_wider(names_from = temps, values_from = marqueur) %>%
    arrange(id) %>% 
    dplyr::select(id, where(~ !all(is.na(.))))
  
  long.data <- cld(as.data.frame(data.wide))
  # assign("long.data", long.data, envir = .GlobalEnv)
  
  # 2 - KML
  options <- NULL
  if (distance == 'fréchet') {
    options <- parALGO(distanceName = "fréchet", distance = function(x, y){FrechetDistance(x, y)})} 
  else if (distance == 'dtw') {
    options <- parALGO(distanceName = "dtw", distance = function(x, y){DTWDistance(x, y)})}
  
  tryCatch({
    if (!is.null(options)) {
      kml(long.data, nbClusters = k, parAlgo = options)} 
    else {
      kml(long.data, nbClusters = k)}}, 
    error = function(e) {message("Erreur KML : ", e$message)
    return(NA)})
  
  # 4 - Résultats
  res.class <- tryCatch({
    as.numeric(getClusters(long.data, k))}, 
    error = function(e) {message("Impossible de récupérer les clusters : ", e$message)
    return(rep(NA, nrow(data.wide)))})
  
  # 5 - Graphique
  res.graph <- plot.clust(data, res.class, paste("KML/", distance, " - Trajectoires selon cluster"))
  
  return(list(clusters = res.class, graphique = res.graph))
}
matClust.fct <- function(id, temps, marqueur, k = 2, 
                         clustering = 'km', distance = "fréchet", 
                         mat.parallel = T, nmes.min = 3) {
  
  # 1 - Données
  id <- as.integer(factor(id))
  data <- filter.mes(data.frame(id, temps, marqueur), id, marqueur, nmes.min)
  
  # 2 - Distance
  mat.dist <- switch(distance,
                     "fréchet" = mat.dist(data, ts = 0.1, mat.parallel),
                     "dtw"     = mat.dist(data, ts = 1e-6, mat.parallel),
                     stop("Méthode de distance non reconnue"))
  
  # 3 - Clustering
  if (clustering == 'km') res.clust <- kmeans(mat.dist, centers = k)
  else if (clustering == 'pam') res.clust <- pam(mat.dist, k = k)
  else stop("Méthode de clustering non reconnue")
  
  res.class <- res.clust$cluster
  
  # 4 - Graphique
  res.graph <- plot.clust(data, res.class, paste0(toupper(distance), "/", clustering, " - Trajectoires selon cluster"))
  
  return(list(clusters = res.class, graphique = res.graph, distance = mat.dist))
}

# Distances manuelles 
mat.dist <- function(data, ts = 0.1, parallel = TRUE, n.cores = NULL) {
  
  # 1 - Données 
  data.split <- list(
    temps = split(data$temps, data$id),
    marqueur = split(data$marqueur, data$id))
  
  n <- length(data.split$temps)
  distances <- matrix(0, n, n)
  comb <- combn(n, 2)
  
  # 2 - Calcul des distances
  if (parallel) {
    
    if (is.null(n.cores)) n.cores <- parallel::detectCores() - 1
    cl <- makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    
    upper.vals <- foreach(p = 1:ncol(comb), .combine = 'c', .packages = "longitudinalData") %dopar% {
      i <- comb[1, p]; j <- comb[2, p]
      suppressMessages(
        distFrechet(
          data.split$temps[[i]], data.split$marqueur[[i]],
          data.split$temps[[j]], data.split$marqueur[[j]],
          FrechetSumOrMax = 'sum', timeScale = ts))}
    
    stopCluster(cl)
    
  } else {
    
    upper.vals <- numeric(ncol(comb))
    for (p in 1:ncol(comb)) {
      i <- comb[1, p]; j <- comb[2, p]
      upper.vals[p] <- suppressMessages(
        distFrechet(
          data.split$temps[[i]], data.split$marqueur[[i]],
          data.split$temps[[j]], data.split$marqueur[[j]],
          FrechetSumOrMax = 'sum', timeScale = ts))
    }
  }
  
  # 3 - Construction de la matrice
  distances[upper.tri(distances)] <- upper.vals
  distances <- distances + t(distances)
  
  return(distances)
}
dtw.dist <- function(data) {

  # 1 - Données
  data.wide <- data %>%
    pivot_wider(names_from = temps, values_from = marqueur) %>%
    arrange(id) %>%
    select(where(~ !all(is.na(.))))
  
  temps_cols <- names(data.wide)
  data.wide <- data.wide[, as.character(sort(as.numeric(temps_cols)))]
  data.list <- split(data.wide, seq_len(nrow(data.wide)))

  # 2 - Matrice des distances
  dist.mat <- dist(data.list, method = distance)
  dist.mat <- as.matrix(dist.mat)
  
  return(dist.mat)
}
dba.fct  <- function(id, temps, marqueur, res.clust, data = NULL, data.list = NULL, data.wide = NULL, nmes.min = 3) {
  
  # 1 - Données
  if (is.null(data)) {
    data <- filter.mes(data.frame(id = id, temps = temps, marqueur = marqueur), id, marqueur, nmes.min)}
  
  if (is.null(data.list)) {
    data.list <- list(
      temps = split(data$temps, data$id),
      marqueur = split(data$marqueur, data$id))}
  
  if (is.null(data.wide)) {
    data.wide <- data %>% pivot_wider(names_from = temps, values_from = marqueur) %>% dplyr::select(-id)
    data.wide <- data.wide[, as.character(sort(as.numeric(names(data.wide))))]
    colnames(data.wide) <- round(as.numeric(colnames(data.wide)), 7)}
  
  last.tps <- data %>% group_by(id) %>% filter(temps == max(temps)) %>% ungroup() %>% pull(temps)
  
  # 2 - Centroids
  k <- length(unique(res.clust))
  cent.y <- vector('list', k)
  cent.t <- vector('list', k)
  cent.id <- vector('list', k)
  cols <- sort(as.numeric(names(data.wide)))
  centroids <- as.data.frame(matrix(NA, k, length(cols), dimnames = list(NULL, cols)))
  
  for (i in 1:k) {
    id.clust <- which(res.clust == i)
    if (length(id.clust) > 1) {
      med_row <- which(rowSums(!is.na(data.wide[id.clust,])) == median(rowSums(!is.na(data.wide[id.clust,]))))[1]
      cent.y[[i]] <- dba(data.list$marqueur[id.clust], centroid = data.list$marqueur[id.clust][[med_row]])
      avail.t <- sort(unique(unlist(data.list$temps[id.clust])))
      avail.t <- avail.t[avail.t <= median(last.tps[id.clust])]
      cent.t[[i]] <- quantile(avail.t, probs = seq(0, 1, length.out = length(cent.y[[i]])), names = FALSE, type = 3)
    } else {
      cent.y[[i]] <- as.numeric(data.wide[id.clust, which(!is.na(data.wide[id.clust,]))])
      cent.t[[i]] <- as.numeric(colnames(data.wide[, which(!is.na(data.wide[id.clust,]))]))
    }
    
    keep <- !duplicated(round(cent.t[[i]], 7))
    cent.t[[i]] <- cent.t[[i]][keep]
    cent.y[[i]] <- cent.y[[i]][keep]
    cent.id[[i]] <- rep(i, length(cent.y[[i]]))
    tps.clust <- as.character(round(cent.t[[i]], 7))
    centroids[i, tps.clust] <- cent.y[[i]]
  }
  
  # 3 - Graphique
  data.plot <- merge(data, data.frame(id = unique(data$id), class = res.clust), by = 'id')
  data.cent <- data.frame(centroid = unlist(cent.id), t = unlist(cent.t), y = unlist(cent.y))
  
  res.plot <- ggplot(data.plot, aes(x = temps, y = marqueur, color = factor(class))) +
    geom_line(aes(group = id), alpha = 0.5) +
    geom_line(data = data.cent, aes(x = t, y = y, group = factor(centroid)), color = 'black', linewidth = 0.8) +
    labs(x = "Temps", y = "Marqueur", color = "Clusters") +
    theme_minimal()
  
  return(list(centroids = centroids, graphique = res.plot))
}

# Évaluation
sil.fct <- function(id, temps, marqueur, res.clust, 
                    matrice.dist = matrix(), distance = 'DTW', nmes.min = 3){
  
  # 1 - Données 
  if(length(matrice.dist) == 1){
    
    id <- as.numeric(id)
    data <- data.frame(id, temps, marqueur) 
    data <- filter.mes(data, id, marqueur, nmes.min)
    data.list <- list(temps = split(data$temps, data$id),
                      marqueur = split(data$marqueur, data$id))
    
    matrice.dist <- switch(distance,
                           'fréchet' = mat.dist(data.list$temps, data.list$marqueur),
                           'DTWB'     = mat.dist(data.list$temps, data.list$marqueur, ts = 1e-6),
                           'DTW'      = {
                             data.wide <- pivot_wider(data, names_from = temps, values_from = marqueur) %>%
                               arrange(id) %>%
                               dplyr::select(where(~ !all(is.na(.))))
                             data.wide <- data.wide[, as.character(sort(as.numeric(names(data.wide))))]
                             as.matrix(dist(data.wide, method = distance))})}
  
  # 2 - Indices de Silhouette
  if(any(is.na(res.clust))) return(NA)
  ind <- silhouette(res.clust, matrice.dist)
  tryCatch(mean(ind[,3]), error = function(e) NA)
}
dbi.fct <- function(id, temps, marqueur, res.clust, 
                    distance = 'DTW', nmes.min = 3){
  
  # 1 - Données 
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur)
  data <- filter.mes(data, id, marqueur, nmes.min)
  data.list <- list(temps = split(data$temps, data$id),
                    marqueur = split(data$marqueur, data$id))
  data.wide <- pivot_wider(data, names_from = temps, values_from = marqueur) %>%
    dplyr::select(-id)
  data.wide <- data.wide[, as.character(sort(as.numeric(names(data.wide))))]
  colnames(data.wide) <- round(as.numeric(colnames(data.wide)),7)
  
  # 2 - Centroïdes
  if(any(is.na(res.clust)) | length(unique(res.clust))==1) return(list(score = NA, centroids = NA))
  k <- length(unique(res.clust))
  res.dba <- dba.fct(id = id, temps = temps, marqueur = marqueur, res.clust = res.clust, 
                     data = data, data.list = data.list, data.wide = data.wide)
  centroids <- res.dba$centroids
  
  # 3 - Distances inter.
  if(distance == 'DTW'){
    dist.inter <- as.matrix(dist(split(centroids, f=1:k), method=distance))} 
  else {
    cent.list <- list(temps=vector('list',k), marqueur=vector('list',k))
    for(i in 1:k){
      cent.list$temps[[i]] <- as.numeric(colnames(centroids[i, which(!is.na(centroids[i,]))]))
      cent.list$marqueur[[i]] <- as.numeric(centroids[i, which(!is.na(centroids[i,]))])}
    if (distance == 'fréchet') dist.inter <- mat.dist(cent.list$temps, cent.list$marqueur)
    else if (distance == "DTWB") dist.inter <- mat.dist(cent.list$temps, cent.list$marqueur, ts=1e-6)}
  
  # 4 - Distances intra.
  var.intra <- numeric(k)
  for (i in 1:k){
    id.clust <- which(res.clust==unique(res.clust)[i])
    dist.intra <- numeric(length(id.clust))
    for (j in seq_along(id.clust)){
      if (distance == "DTW") dist.intra[j] <- dtw(data.wide[id.clust[j],], centroids[i,])$distance
      else dist.intra[j] <- suppressMessages(
        distFrechet(data.list$temps[[id.clust[j]]], data.list$marqueur[[id.clust[j]]],
                    cent.list$temps[[i]], cent.list$marqueur[[i]], 
                    FrechetSumOrMax='sum', timeScale = if (distance == 'DTWB') 1e-6 else 0.1))}
    var.intra[i] <- mean(dist.intra, na.rm=TRUE)}
  
  # 5 - Indice de DB
  dbi <- sapply(1:k, function(i){
    max(sapply(setdiff(1:k,i), function(j) (var.intra[i]+var.intra[j])/dist.inter[i,j]), na.rm=TRUE)
  })
  
  list(score = mean(dbi, na.rm=TRUE), centroids = res.dba$graphique)
}
chi.fct <- function(id, temps, marqueur, res.clust, 
                    matrice.dist = matrix(), distance = 'DTW', nmes.min = 3){
  
  # 1 - Données
  if(length(matrice.dist) == 1){
    id <- as.numeric(id)
    data <- data.frame(id, temps, marqueur)
    data <- filter.mes(data, id, marqueur, nmes.min)
    data.list <- list(temps = split(data$temps, data$id),
                      marqueur = split(data$marqueur, data$id))
    matrice.dist <- switch(distance,
                           'fréchet' = mat.dist(data.list$temps, data.list$marqueur),
                           'DTWB'    = mat.dist(data.list$temps, data.list$marqueur, ts=1e-6),
                           'DTW'     = {
                             data.wide <- pivot_wider(data, names_from=temps, values_from=marqueur)
                             as.matrix(dist(data.wide))})}
  
  # 2 - Indice de CH
  if(any(is.na(res.clust)) | length(unique(res.clust))==1) return(NA)
  cluster.stats(matrice.dist, res.clust)$ch
}
indices.fct <- function(id, temps, marqueur, methodes, 
                        matrice.dist = matrix(), distance = 'DTW', nmes.min = 3, 
                        indices = c('Silhouette','Davies-Bouldin','Calinski-Harabasz')){
  
  # 1 - Données
  if(length(matrice.dist)==1){
    id <- as.numeric(id)
    data <- data.frame(id, temps, marqueur)
    data <- filter.mes(data, id, marqueur, nmes.min)
    data.list <- list(temps = split(data$temps, data$id),
                      marqueur = split(data$marqueur, data$id))
    matrice.dist <- switch(distance,
                           'fréchet' = mat.dist(data.list$temps, data.list$marqueur),
                           'DTWB'    = mat.dist(data.list$temps, data.list$marqueur, ts=1e-6),
                           'DTW'     = as.matrix(dist(pivot_wider(data, names_from=temps, values_from=marqueur))))}
  
  # 2 - Comparaison
  m <- length(methodes)
  scores <- list(si=numeric(m), dbi=numeric(m), chi=numeric(m))
  noms <- names(methodes)
  for(i in seq_along(methodes)){
    if('Silhouette' %in% indices) scores$si[i] <- sil.fct(res.clust = methodes[[i]], matrice.dist = matrice.dist, distance=distance)
    if('Davies-Bouldin' %in% indices) scores$dbi[i] <- dbi.fct(id, temps, marqueur, methodes[[i]], distance)$score
    if('Calinski-Harabasz' %in% indices) scores$chi[i] <- chi.fct(id, temps, marqueur, methodes[[i]], matrice.dist, distance)}
  
  plots.scores <- lapply(1:length(scores), function(i){
    data.plot <- data.frame(methode=noms, score=scores[[i]])
    y.lim <- switch(i, c(-1,1), if(distance=='fréchet') c(0,30) else c(0,20), if(distance=='fréchet') c(0,1000) else c(0,100))
    indice.nom <- c('Silhouette','Davies-Bouldin','Calinski-Harabasz')[i]
    ggplot(data.plot, aes(x=methode, y=score)) +
      geom_bar(stat="identity", fill="steelblue", color="steelblue", width=0.6) +
      geom_text(aes(label=round(score,2)), vjust=-1) +
      theme_minimal() +
      labs(title=paste("Indice de", indice.nom), x='', y='') +
      coord_cartesian(ylim=y.lim) +
      theme(axis.text.x = element_text(angle=75,hjust=1,size=8))
  })
  
  list(scores=scores, graphique=plots.scores)
}
eval.fct    <- function(cluster.ref, clusters.test){
  
  if(is.list(clusters.test)){
    
    resultats <- sapply(clusters.test, 
                        function(x) if(any(is.na(x))) NA else adjustedRandIndex(cluster.ref, x))
    noms <- names(clusters.test)
    data.plot <- data.frame(methodes = noms, scores = resultats)
    
    graphique <- ggplot(data.plot, aes(x = methodes, y = scores)) +
      geom_bar(stat="identity", fill="steelblue", color="steelblue", width=0.6) +
      geom_text(aes(label=round(scores,2)), vjust=-1) +
      theme_minimal() +
      labs(title="Indice ARI", x='', y='') +
      coord_cartesian(ylim=c(0,1.2)) +
      theme(axis.text.x=element_text(angle=75,hjust=1,size=8))} 
  
  else {
    resultats <- if(any(is.na(clusters.test))) NA else adjustedRandIndex(cluster.ref, clusters.test)
    graphique <- NULL}
  
  list(scores=resultats, graphique=graphique)
}

# Comparaison
sim.fct  <- function(class1, class2){
  
  if (any(is.na(class1)) | any(is.na(class2))) return(NA)
  
  mat.sim <- lapply(list(class1, class2), function(cl) {
    mat <- outer(cl, cl, FUN = "==") * 1
    mat[lower.tri(mat, diag = TRUE)] <- 0
    mat
  })
  
  sim.score <- 1 - sum((mat.sim[[1]] - mat.sim[[2]])^2) / 
                    (length(class1) * (length(class1) - 1) / 2)
  
  return(sim.score)
}
comp.fct <- function(methodes = list(), type = 'sim', color = 1){
  
  # 1 - Matrice des scores
  n <- length(methodes)
  scoresSIM <- diag(1, n)
  scoresSIM[upper.tri(scoresSIM)] <- NA
  scoresARI <- diag(1, n)
  scoresARI[upper.tri(scoresARI)] <- NA
  
  for (i in 2:n){
    for (j in 1:(i-1)){
      scoresSIM[i,j] <- sim.fct(methodes[[i]], methodes[[j]])
      scoresARI[i,j] <- eval.fct(methodes[[i]], methodes[[j]])$scores}}
  
  # 2 - Graphiques
  noms <- names(methodes)
  rownames(scoresSIM) <- colnames(scoresSIM) <- noms
  rownames(scoresARI) <- colnames(scoresARI) <- noms
  
  colors <- switch(as.character(color),
                   "1" = viridis(50, option = "G", begin = 1, end = 0.5),
                   "2" = viridis(50, option = "G", begin = 0.5, end = 1),
                   "3" = viridis(50, begin = 1, end = 0.5, option = "A"))
  
  graphiqueSIM <- Heatmap(scoresSIM,
                          name = " ",
                          column_title = "Matrice de similarité entre les méthodes (SIM)",
                          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          col = colors,
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          row_names_side = "left", row_names_gp = gpar(fontsize = 7),
                          column_names_side = "bottom", column_names_gp = gpar(fontsize = 7),
                          column_names_rot = 45, na_col = 'white',
                          cell_fun = function(j, i, x, y, width, height, fill) { 
                            if (!is.na(scoresSIM[i, j])) 
                              grid.text(round(scoresSIM[i, j], 2), x, y, gp = gpar(fontsize = 6, col = "black"))
                          })
  
  graphiqueARI <- Heatmap(scoresARI,
                          name = " ",
                          column_title = "Matrice de similarité entre les méthodes (ARI)",
                          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                          col = colors,
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          row_names_side = "left", row_names_gp = gpar(fontsize = 7),
                          column_names_side = "bottom", column_names_gp = gpar(fontsize = 7),
                          column_names_rot = 45, na_col = 'white',
                          cell_fun = function(j, i, x, y, width, height, fill) { 
                            if (!is.na(scoresARI[i, j])) 
                              grid.text(round(scoresARI[i, j], 2), x, y, gp = gpar(fontsize = 6, col = "black"))
                          })
  
  list(scores = list(SIM = scoresSIM, ARI = scoresARI),
       graphique = list(SIM = graphiqueSIM, ARI = graphiqueARI))
}

