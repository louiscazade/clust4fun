
source('Fonctions/FctSimples.R')

comparaison.denseOLD <- function(id, temps, marqueur, 
                                 cluster.ref, k.min = 2, k.max = 5, 
                                 distance = 'dtw', nmes.min = 3){
  
  k.test <- k.max - k.min + 1
  
  methodes <- list(acpf.km9      = vector('list', k.test),
                   acpf.gmm9     = vector('list', k.test),
                   hlme.quad     = vector('list', k.test),
                   hlme.s3       = vector('list', k.test),
                   kml.eucl      = vector('list', k.test),
                   kml.frech     = vector('list', k.test),
                   kml.dtw       = vector('list', k.test), 
                   mat.dtw       = vector('list', k.test),
                   mat.frech     = vector('list', k.test))
  
  # 1 - Clustering pour toutes les valeurs de k 
  
  i <- 1
  for (k in k.min:k.max){
    
    # ACPF / PVE : 0.99
    methodes$acpf.km9[[i]]  <- acpf.fct(id, temps, marqueur, k, clustering = 'km', nmes.min = nmes.min)$clusters
    methodes$acpf.gmm9[[i]] <- acpf.fct(id, temps, marqueur, k, clustering = 'gmm', nmes.min = nmes.min)$clusters
    
    # LCMM
    methodes$hlme.quad[[i]] <- hlme.fct(id, temps, marqueur, k, regression = 'quadratique', nmes.min = nmes.min)$clusters
    methodes$hlme.s3[[i]]   <- hlme.fct(id, temps, marqueur, k, regression = 'splines', nmes.min = nmes.min)$clusters
    
    # KML
    methodes$kml.eucl[[i]]  <- kml.fct(id, temps, marqueur, k, distance = 'euclidienne', nmes.min = nmes.min)$clusters
    methodes$kml.frech[[i]] <- kml.fct(id, temps, marqueur, k, distance = 'fréchet', nmes.min = nmes.min)$clusters
    methodes$kml.dtw[[i]]   <- kml.fct(id, temps, marqueur, k, distance = 'dtw', nmes.min = nmes.min)$clusters
    
    # Matrice de distance
    methodes$mat.dtw[[i]]   <- matClust.fct(id, temps, marqueur, k, distance = 'dtw', nmes.min = nmes.min)$clusters
    methodes$mat.frech[[i]] <- matClust.fct(id, temps, marqueur, k, distance = 'fréchet', nmes.min = nmes.min)$clusters
    
    i <- i+1
  }
  
  # 2.1 - Calcul de la matrice de distance (1x pour toutes)
  
  id <- as.numeric(id)
  data <- data.frame(id, temps, marqueur) 
  data <- filter.mes(data, id, marqueur, nmes.min)
  
  if (distance == 'fréchet'){matrice.dist <- mat.dist(data, parallel = F)} 
  else if (distance == 'dtw') {matrice.dist <- dtw.dist(data)}
  
  # 2.2 - Calcul du k optimal pour chaque méthode
  
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
    methodes$acpf.km9 <- acpf.fct(id, temps, marqueur, k.true, clustering = 'km', nmes.min = nmes.min))["elapsed"]
  temps.exe$acpf.gmm9 <- system.time(
    methodes$acpf.gmm9 <- acpf.fct(id, temps, marqueur, k.true, clustering = 'gmm', nmes.min = nmes.min))["elapsed"]
  
  # LCMM
  temps.exe$hlme.quad <- system.time(
    methodes$hlme.quad <- hlme.fct(id, temps, marqueur, k.true, regression = 'quadratique', nmes.min = nmes.min))["elapsed"]
  temps.exe$hlme.s3 <- system.time(
    methodes$hlme.s3 <- hlme.fct(id, temps, marqueur, k.true, regression = 'splines', nmes.min = nmes.min))["elapsed"]
  
  # KML
  temps.exe$kml.eucl <- system.time(
    methodes$kml.eucl <- kml.fct(id, temps, marqueur, k.true, distance = 'euclidienne', nmes.min = nmes.min))["elapsed"]
  temps.exe$kml.frech <- system.time(
    methodes$kml.frech <- kml.fct(id, temps, marqueur, k.true, distance = 'fréchet', nmes.min = nmes.min))["elapsed"]
  temps.exe$kml.dtw <- system.time(
    methodes$kml.dtw <- kml.fct(id, temps, marqueur, k.true, distance = 'dtw', nmes.min = nmes.min))["elapsed"]
  
  # Matrices de distance
  temps.exe$mat.dtw <- system.time(
    methodes$mat.dtw <- matClust.fct(id, temps, marqueur, k.true, distance = 'dtw', nmes.min = nmes.min))["elapsed"]
  temps.exe$mat.frech <- system.time(
    methodes$mat.frech <- matClust.fct(id, temps, marqueur, k.true, distance = 'fréchet', nmes.min = nmes.min))["elapsed"]
  
  clusters <- list(acpf.km9    = methodes$acpf.km9$clusters,
                   acpf.gmm9   = methodes$acpf.gmm9$clusters,
                   hlme.quad   = methodes$hlme.quad$clusters,
                   hlme.s3     = methodes$hlme.s3$clusters,
                   kml.eucl    = methodes$kml.eucl$clusters,
                   kml.frech   = methodes$kml.frech$clusters,
                   kml.dtw     = methodes$kml.dtw$clusters,
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
                              cluster.ref, k.min = 2, k.max = 5, 
                              distance = 'fréchet', nmes.min = 3,
                              parallel = TRUE, n.cores = NULL) {
  # 0 - Initialisation
  if (parallel) {
    if (is.null(n.cores)) n.cores <- parallel::detectCores() - 2} 
  else {n.cores <- 1}
  
  methodes.info <- list(
    acpf.km9  = list(fun = acpf.fct, args = list(clustering = "km")),
    acpf.gmm9 = list(fun = acpf.fct, args = list(clustering = "gmm")),
    hlme.quad = list(fun = hlme.fct, args = list(regression = "quadratique")),
    hlme.s3   = list(fun = hlme.fct, args = list(regression = "splines")),
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
  
  run.meth <- function(meth_name, k_val) {
    meth <- methodes.info[[meth_name]]
    args <- c(list(id = id, temps = temps, marqueur = marqueur, k = k_val, nmes.min = nmes.min), meth$args)
    do.call(meth$fun, args)$clusters}
  
  res.list <- mclapply(seq_len(nrow(combinaisons)), function(i){
    m <- combinaisons$methode[i]
    k_val <- combinaisons$k[i]
    idx <- which(k.seq == k_val)
    clusters <- run.meth(m, k_val)
    list(methode = m, k.idx = idx, clusters = clusters)}, mc.cores = n.cores)
  
  for (res in res.list) methodes[[res$methode]][[res$k.idx]] <- res$clusters
  
  # 2 - Calcul du k optimal
  data <- data.frame(id = as.numeric(id), temps, marqueur) %>%
    group_by(id) %>%
    filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    ungroup()
  
  if (distance == 'fréchet'){
    matrice.dist <- mat.dist(data, parallel = parallel, n.cores = n.cores)}
  else if (distance == 'dtw'){
    matrice.dist <- dtw.dist(data)}
  
  # 2.2 - Calcul des indices
  k.test <- length(k.seq)
  index <- list(silhouette = list(), dbi.cent = list(), chi = list())
  for (crit in names(index)) index[[crit]] <- lapply(methodes, function(x) numeric(k.test))
  
  comb.indices <- expand.grid(methode = names(methodes), k.idx = seq_along(k.seq), stringsAsFactors = FALSE)
  
  res.index <- mclapply(seq_len(nrow(comb.indices)), function(i){
    m <- comb.indices$methode[i]
    idx <- comb.indices$k.idx[i]
    clust <- methodes[[m]][[idx]]
    list(methode = m, k.idx = idx,
         silhouette = sil.fct(matrice.dist = matrice.dist, res.clust = clust),
         dbi.cent   = dbiCent.fct(id, temps, marqueur, res.clust = clust, distance = distance),
         chi        = chi.fct(matrice.dist = matrice.dist, res.clust = clust))}, 
    mc.cores = n.cores)
  
  for (res in res.index) {
    index$silhouette[[res$methode]][res$k.idx] <- res$silhouette
    index$dbi.cent[[res$methode]][res$k.idx]   <- res$dbi.cent
    index$chi[[res$methode]][res$k.idx]        <- res$chi}
  
  for (m in names(methodes)) {
    index$silhouette[[m]] <- which.max(index$silhouette[[m]]) + k.min - 1
    index$dbi.cent[[m]]   <- which.min(index$dbi.cent[[m]]) + k.min - 1
    index$chi[[m]]        <- which.max(index$chi[[m]]) + k.min - 1}
  
  # 3 - Clustering avec le k optimal
  k.true <- length(unique(cluster.ref))
  
  run.true <- function(meth_name) {
    meth <- methodes.info[[meth_name]]
    args <- c(list(id = id, temps = temps, marqueur = marqueur, k = k.true, nmes.min = nmes.min), meth$args)
    t <- system.time(res <- do.call(meth$fun, args))["elapsed"]
    list(clusters = res$clusters, time = t)}
  
  res.true <- mclapply(names(methodes.info), run.true, mc.cores = n.cores)
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
  
  return(list(
    evaluation = list(clusters = eval.clusters$scores, indices = eval.indexes$scores),
    comparaison = comp.methodes$scores,
    k.true = k.true,
    k.opt = index,
    temps = temps.exe
  ))
}

comparaison.dense_lite <- function(id, temps, marqueur, 
                              cluster.ref, k.min = 2, k.max = 5, 
                              distance = 'fréchet', nmes.min = 3,
                              parallel = TRUE, n.cores = NULL) {
  # 0 - Initialisation
  if (parallel) {
    if (is.null(n.cores)) n.cores <- parallel::detectCores() - 2} 
  else {n.cores <- 1}
  
  methodes.info <- list(
    acpf.gmm9 = list(fun = acpf.fct, args = list(clustering = "gmm")),
    hlme.s3   = list(fun = hlme.fct, args = list(regression = "splines")),
    kml.eucl  = list(fun = kml.fct, args = list(distance = "euclidienne")),
    mat.frech = list(fun = matClust.fct, args = list(distance = "fréchet", mat.parallel = F)))
  
  # 1 - Clustering pour toutes les valeurs de k
  k.seq <- k.min:k.max
  methodes <- lapply(names(methodes.info), function(x) vector("list", length(k.seq)))
  names(methodes) <- names(methodes.info)
  
  combinaisons <- expand.grid(methode = names(methodes.info), k = k.seq, stringsAsFactors = FALSE)
  
  run.meth <- function(meth_name, k_val) {
    meth <- methodes.info[[meth_name]]
    args <- c(list(id = id, temps = temps, marqueur = marqueur, k = k_val, nmes.min = nmes.min), meth$args)
    do.call(meth$fun, args)$clusters}
  
  res.list <- mclapply(seq_len(nrow(combinaisons)), function(i){
    m <- combinaisons$methode[i]
    k_val <- combinaisons$k[i]
    idx <- which(k.seq == k_val)
    clusters <- run.meth(m, k_val)
    list(methode = m, k.idx = idx, clusters = clusters)}, mc.cores = n.cores)
  
  for (res in res.list) methodes[[res$methode]][[res$k.idx]] <- res$clusters
  
  # 2 - Calcul du k optimal
  data <- data.frame(id = as.numeric(id), temps, marqueur) %>%
    group_by(id) %>%
    filter(sum(!is.na(marqueur)) >= nmes.min) %>%
    ungroup()
  
  if (distance == 'fréchet'){
    matrice.dist <- mat.dist(data, parallel = parallel, n.cores = n.cores)}
  else if (distance == 'dtw'){
    matrice.dist <- dtw.dist(data)}
  
  # 2.2 - Calcul des indices
  k.test <- length(k.seq)
  index <- list(silhouette = list(), dbi.cent = list(), chi = list())
  for (crit in names(index)) index[[crit]] <- lapply(methodes, function(x) numeric(k.test))
  
  comb.indices <- expand.grid(methode = names(methodes), k.idx = seq_along(k.seq), stringsAsFactors = FALSE)
  
  res.index <- mclapply(seq_len(nrow(comb.indices)), function(i){
    m <- comb.indices$methode[i]
    idx <- comb.indices$k.idx[i]
    clust <- methodes[[m]][[idx]]
    list(methode = m, k.idx = idx,
         silhouette = sil.fct(matrice.dist = matrice.dist, res.clust = clust),
         dbi.cent   = dbiCent.fct(id, temps, marqueur, res.clust = clust, distance = distance),
         chi        = chi.fct(matrice.dist = matrice.dist, res.clust = clust))}, 
    mc.cores = n.cores)
  
  for (res in res.index) {
    index$silhouette[[res$methode]][res$k.idx] <- res$silhouette
    index$dbi.cent[[res$methode]][res$k.idx]   <- res$dbi.cent
    index$chi[[res$methode]][res$k.idx]        <- res$chi}
  
  for (m in names(methodes)) {
    index$silhouette[[m]] <- which.max(index$silhouette[[m]]) + k.min - 1
    index$dbi.cent[[m]]   <- which.min(index$dbi.cent[[m]]) + k.min - 1
    index$chi[[m]]        <- which.max(index$chi[[m]]) + k.min - 1}
  
  # 3 - Clustering avec le k optimal
  k.true <- length(unique(cluster.ref))
  
  run.true <- function(meth_name) {
    meth <- methodes.info[[meth_name]]
    args <- c(list(id = id, temps = temps, marqueur = marqueur, k = k.true, nmes.min = nmes.min), meth$args)
    t <- system.time(res <- do.call(meth$fun, args))["elapsed"]
    list(clusters = res$clusters, time = t)}
  
  res.true <- mclapply(names(methodes.info), run.true, mc.cores = n.cores)
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
  
  return(list(
    evaluation = list(clusters = eval.clusters$scores, indices = eval.indexes$scores),
    comparaison = comp.methodes$scores,
    k.true = k.true,
    k.opt = index,
    temps = temps.exe
  ))
}








