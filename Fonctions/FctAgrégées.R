
source('FctSimples.R')

comparaison.denseOLD <- function(id, temps, marqueur, 
                              cluster.ref, k.min = 2, k.max = 4, 
                              distance = 'DTW', nmes.min = 3){
  
  k.test <- k.max - k.min + 1
  
  methodes <- list(acpf.km9      = vector('list', k.test),
                   acpf.gmm9     = vector('list', k.test),
                   lcmm.quad     = vector('list', k.test),
                   lcmm.s3       = vector('list', k.test),
                   kml.frech     = vector('list', k.test),
                   kml.dtw       = vector('list', k.test),
                   dtw.pam       = vector('list', k.test), 
                   dtw.dba       = vector('list', k.test), 
                   mat.dtw       = vector('list', k.test),
                   mat.frech     = vector('list', k.test))
  
  
  # 1 - Clustering pour toutes les valeurs de k 
  
  i <- 1
  for (k.test in k.min:k.max){
    
    # ACPF / PVE : 0.99
    methodes$acpf.km9[[i]]  <- acpfm.fct(id, temps, marqueur, k.test, clustering = 'km', nmes.min = nmes.min)$clusters
    methodes$acpf.gmm9[[i]] <- acpfm.fct(id, temps, marqueur, k.test, clustering = 'gmm', nmes.min = nmes.min)$clusters
    
    # LCMM
    methodes$lcmm.quad[[i]] <- lcmm.fct(id, temps, marqueur, k.test, regression = 'quadratique', nmes.min = nmes.min)$clusters
    methodes$lcmm.s3[[i]]   <- lcmm.fct(id, temps, marqueur, k.test, regression = 'splines', nmes.min = nmes.min)$clusters
    
    # KML
    methodes$kml.frech[[i]] <- kml.fct(id, temps, marqueur, k.test, distance = 'frechet', nmes.min = nmes.min)$clusters
    methodes$kml.dtw[[i]]   <- kml.fct(id, temps, marqueur, k.test, distance = 'dtw', nmes.min = nmes.min)$clusters
    
    # DTW
    methodes$dtw.pam[[i]]   <- dtw.fct(id, temps, marqueur, k.test, centroid = 'pam', nmes.min = nmes.min)$clusters
    methodes$dtw.dba[[i]]   <- dtw.fct(id, temps, marqueur, k.test, centroid = 'dba', nmes.min = nmes.min)$clusters
    
    # Matrice de distance
    methodes$mat.dtw[[i]]   <- matClust.fct(id, temps, marqueur, k.test, distance = 'dtw', nmes.min = nmes.min)$clusters
    methodes$mat.frech[[i]] <- matClust.fct(id, temps, marqueur, k.test, distance = 'frechet', nmes.min = nmes.min)$clusters
    
    i <- i+1
  }
  
  # 2.1 - Calcul de la matrice de distance (1x pour toutes)
  
  if (distance == 'fréchet'){
    
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
    matrice.dist <- mat.frechet(data.frechet$temps, data.frechet$marqueur) 
    
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
    
    # 2 - Matrice de distance 
    data.list <- split(data.wide, f = seq_len(nrow(data.wide)))
    matrice.dist <- dist(data.list, method = 'DTW')
    matrice.dist <- as.matrix(matrice.dist)
  }
  
  # 2.2 - Calcul du k optimal pour chaque méthode
  
  k.test <- k.max - k.min + 1
  
  index <- list(silhouette = NULL,
                dbi.cent   = NULL,
                chi        = NULL)
  
  for (crit in names(index)) {
    index[[crit]] <- setNames(vector("list", length(methodes)), names(methodes))
    for (meth in names(methodes)) {
      index[[crit]][[meth]] <- numeric(k.test)
    }
  }
  
  for (i in 1:length(methodes)){
    for (id.k in 1:k.test){
      index$silhouette[[i]][id.k] <- sil.fct(matrice.dist = matrice.dist, res.clust = methodes[[i]][[id.k]])
      index$dbi.cent[[i]][id.k]   <- dbiCent.fct(id, temps, marqueur, res.clust = methodes[[i]][[id.k]], distance = distance)
      index$chi[[i]][id.k]        <- chi.fct(matrice.dist = matrice.dist, res.clust = methodes[[i]][[id.k]])
    }
    index$silhouette[[i]] <- which.max(index$silhouette[[i]])+k.min-1
    index$dbi.cent[[i]]   <- which.min(index$dbi.cent[[i]])+k.min-1
    index$chi[[i]]        <- which.max(index$chi[[i]])+k.min-1
  }
  
  # 3 - Clustering avec le véritable k
  
  k.true <- length(unique(cluster.ref))
  temps.exe <- list()
  
  # ACPF
  temps.exe$acpf.km9  <- system.time(
    methodes$acpf.km9 <- acpfm.fct(id, temps, marqueur, k.true, clustering = 'km', nmes.min = nmes.min))["elapsed"]
  temps.exe$acpf.gmm9 <- system.time(
    methodes$acpf.gmm9 <- acpfm.fct(id, temps, marqueur, k.true, clustering = 'gmm', nmes.min = nmes.min))["elapsed"]
  
  # LCMM
  temps.exe$lcmm.quad <- system.time(
    methodes$lcmm.quad <- lcmm.fct(id, temps, marqueur, k.true, regression = 'quadratique', nmes.min = nmes.min))["elapsed"]
  temps.exe$lcmm.s3 <- system.time(
    methodes$lcmm.s3 <- lcmm.fct(id, temps, marqueur, k.true, regression = 'splines', nmes.min = nmes.min))["elapsed"]
  
  # KML
  temps.exe$kml.frech <- system.time(
    methodes$kml.frech <- kml.fct(id, temps, marqueur, k.true, distance = 'frechet', nmes.min = nmes.min))["elapsed"]
  temps.exe$kml.dtw <- system.time(
    methodes$kml.dtw <- kml.fct(id, temps, marqueur, k.true, distance = 'dtw', nmes.min = nmes.min))["elapsed"]
  
  # DTW
  temps.exe$dtw.pam <- system.time(
    methodes$dtw.pam <- dtw.fct(id, temps, marqueur, k.true, centroid = 'pam', nmes.min = nmes.min))["elapsed"]
  temps.exe$dtw.dba <- system.time(
    methodes$dtw.dba <- dtw.fct(id, temps, marqueur, k.true, centroid = 'dba', nmes.min = nmes.min))["elapsed"]
  
  # Matrices de distance
  temps.exe$mat.dtw <- system.time(
    methodes$mat.dtw <- matClust.fct(id, temps, marqueur, k.true, distance = 'dtw', nmes.min = nmes.min))["elapsed"]
  temps.exe$mat.frech <- system.time(
    methodes$mat.frech <- matClust.fct(id, temps, marqueur, k.true, distance = 'frechet', nmes.min = nmes.min))["elapsed"]
  
  clusters <- list(acpf.km9    = methodes$acpf.km9$clusters,
                   acpf.gmm9   = methodes$acpf.gmm9$clusters,
                   lcmm.quad   = methodes$lcmm.quad$clusters,
                   lcmm.s3     = methodes$lcmm.s3$clusters,
                   kml.frech   = methodes$kml.frech$clusters,
                   kml.dtw     = methodes$kml.dtw$clusters,
                   dtw.pam     = methodes$dtw.pam$clusters,
                   dtw.dba     = methodes$dtw.dba$clusters,
                   mat.dtw     = methodes$mat.dtw$clusters,
                   mat.frech   = methodes$mat.frech$clusters)
  
  # 4 - Évaluation du clustering
  eval.clusters <- eval.fct(cluster.ref, clusters)
  eval.indexes <- indices.fct(id, temps, marqueur,
                              matrice.dist = matrice.dist,
                              methodes = clusters,
                              distance = distance)
  
  # 5 - Comparaison entre méthodes
  comp.methodes <- comp.fct(clusters)
  
  return(list(evaluation = list(clusters = eval.clusters$scores,
                                indices = eval.indexes$scores),
              comparaison = comp.methodes$scores,
              k.true = k.true,
              k.opt = index,
              temps = temps.exe))
}

comparaison.dense <- function(id, temps, marqueur, 
                              cluster.ref, k.min = 2, k.max = 4, 
                              distance = 'fréchet', nmes.min = 3,
                              parallel = TRUE, n.cores = NULL){
  
  # 0 - Initialisation de la parallélisation
  if (parallel){
    if (is.null(n.cores)) n.cores <- parallel::detectCores() - 1
    plan(multisession, workers = n.cores)
  } else plan(sequential)
  
  methodes.info <- list(
    acpf.km9  = list(fun = acpfm.fct, args = list(clustering = "km")),
    acpf.gmm9 = list(fun = acpfm.fct, args = list(clustering = "gmm")),
    lcmm.quad = list(fun = lcmm.fct, args = list(regression = "quadratique")),
    lcmm.s3   = list(fun = lcmm.fct, args = list(regression = "splines")),
    kml.eucl  = list(fun = kml.fct, args = list(distance = "euclidienne")),
    kml.frech = list(fun = kml.fct, args = list(distance = "fréchet")),
    kml.dtw   = list(fun = kml.fct, args = list(distance = "dtw")),
    mat.dtw   = list(fun = matClust.fct, args = list(distance = "dtw", mat.parallel = F)),
    mat.frech = list(fun = matClust.fct, args = list(distance = "fréchet", mat.parallel = F)))
  
  # 1 - Clustering pour toutes les valeurs de k 
  k.seq <- k.min:k.max
  methodes <- lapply(names(methodes.info), function(x) vector("list", length(k.seq)))
  names(methodes) <- names(methodes.info)
  
  combinaisons <- expand.grid(methode = names(methodes.info), k = k.seq, stringsAsFactors = FALSE)
  
  run.meth <- function(meth, k){
    fun <- meth$fun
    args <- c(list(id = id, temps = temps, marqueur = marqueur, k = k, nmes.min = nmes.min), meth$args)
    do.call(fun, args)$clusters}
  
  res.list <- future_lapply(seq_len(nrow(combinaisons)), function(i){
    m <- combinaisons$methode[i]
    k.val <- combinaisons$k[i]
    idx <- which(k.seq == k.val)
    list(methode = m, k.idx = idx, clusters = run.meth(methodes.info[[m]], k.val))})
  
  for(res in res.list) methodes[[res$methode]][[res$k.idx]] <- res$clusters
  
  # 2 - Calcul du k optimal
  data <- data.frame(id = as.numeric(id), temps, marqueur) %>%
    group_by(id) %>%
    filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    ungroup()
  
  if (distance %in% c('fréchet','dtwb')) {
    if (distance == 'fréchet'){
      matrice.dist <- mat.dist(data, parallel = parallel, n.cores = n.cores)} 
    else if (distance == 'dtwb'){
      matrice.dist <- mat.dist(data, ts = 1e-6, parallel = parallel, n.cores = n.cores)}} 
  else if (distance == 'dtw') {
    matrice.dist <- dtw.dist(data, nmes.min = nmes.min)}
  
  # 2.2 Calcul des indices
  k.test <- length(k.seq)
  index <- list(silhouette = list(), dbi.cent = list(), chi = list())
  for(crit in names(index)) index[[crit]] <- lapply(methodes, function(x) numeric(k.test))
  
  comb.indices <- expand.grid(methode = names(methodes), k.idx = seq_along(k.seq), stringsAsFactors = FALSE)
  
  res.index <- future_lapply(seq_len(nrow(comb.indices)), function(i){
    m <- comb.indices$methode[i]
    idx <- comb.indices$k.idx[i]
    clust <- methodes[[m]][[idx]]
    list(methode = m,
         k.idx = idx,
         silhouette = sil.fct(matrice.dist = matrice.dist, res.clust = clust),
         dbi.cent = dbiCent.fct(id, temps, marqueur, res.clust = clust, distance = distance),
         chi = chi.fct(matrice.dist = matrice.dist, res.clust = clust))})
  
  for(res in res.index){
    index$silhouette[[res$methode]][res$k.idx] <- res$silhouette
    index$dbi.cent[[res$methode]][res$k.idx]   <- res$dbi.cent
    index$chi[[res$methode]][res$k.idx]        <- res$chi}
  
  for(m in names(methodes)){
    index$silhouette[[m]] <- which.max(index$silhouette[[m]]) + k.min - 1
    index$dbi.cent[[m]]   <- which.min(index$dbi.cent[[m]]) + k.min - 1
    index$chi[[m]]        <- which.max(index$chi[[m]]) + k.min - 1}
  
  # 3 - Clustering avec le k optimal
  k.true <- length(unique(cluster.ref))
  run.true <- function(meth){
    fun <- meth$fun
    args <- c(list(id = id, temps = temps, marqueur = marqueur, k = k.true, nmes.min = nmes.min), meth$args)
    t <- system.time(res <- do.call(fun, args))["elapsed"]
    list(clusters = res$clusters, time = t)}
  
  res.true <- future_lapply(methodes.info, run.true)
  names(res.true) <- names(methodes.info)
  methodes <- lapply(res.true, `[[`, "clusters")
  temps.exe <- lapply(res.true, `[[`, "time")
  
  # 4 - Évaluation du clustering
  clusters <- methodes
  eval.clusters <- eval.fct(cluster.ref, clusters)
  eval.indexes <- indices.fct(id, temps, marqueur, matrice.dist = matrice.dist,
                              methodes = clusters, distance = distance)
  
  # 5 - Comparaison du clustering
  comp.methodes <- comp.fct(clusters)
  
  return(list(evaluation = list(clusters = eval.clusters$scores,
                                indices = eval.indexes$scores),
              comparaison = comp.methodes$scores,
              k.true = k.true,
              k.opt = index,
              temps = temps.exe))
}




