
# Packages
library(mvtnorm)
library(nlraa)

generation.data <- function(n, nmes, tmax, type = c("sparse", "dense"), scenario = 1, 
                            ngroups = 2, prob.group = NULL,
                            rho = 0.7, sd_bruit = 0.1){
  
  type <- match.arg(type)
  
  # 1 - Initialisation
  y <- matrix(NA, nrow = nmes, ncol = n)
  t <- matrix(0, nrow = nmes, ncol = n)
  tdef <- seq(0, tmax, length.out = nmes)
  
  # 2 - Génération des temps
  if (type == "sparse") {
    for (i in seq_len(n)) {
      t[, i] <- tdef + c(0, runif(nmes - 1, -tmax/(2*(nmes - 1)), tmax/(2*(nmes - 1))))}}
  else if (type == "dense"){
    for (i in seq_len(n)) {
      t[, i] <- tdef} # Pas de bruit autour des mesures
  } 
  
  # 3 - Attribution des groupes 
  if (type == "sparse"){ngroups <- ifelse(scenario %in% c(1,2,3), 2, 3)}
  if (is.null(prob.group)){prob.group <- rep(1/ngroups, ngroups)}
  groups <- sample(0:(ngroups-1), n, replace = TRUE, prob = prob.group)
  
  # 4 - Bruit AR (dense)
  bruit_AR <- function(nmes, rho, sd){
    e <- rnorm(nmes, 0, sd)
    for (i in 2:nmes) e[i] <- rho * e[i-1] + sqrt(1 - rho^2) * e[i]
    return(e)
  }
  
  # 5 - Génération des données selon scénario et type
  for (i in seq_len(n)) {
    g <- groups[i]
    
    if (type == "dense") { 
     
      # Quadratique
      if (scenario == 1) { 
        if (g == 0) {
          a <- -abs(0.03 + rnorm(1, 0, 0.02))
          b <- -abs(2 + rnorm(1, 0, 0.3))
          c0 <- 500 + rnorm(1, 0, 10)
        } else if (g == 1) {
          a <- -abs(0.01 + rnorm(1, 0, 0.02))
          b <- -abs(3 + rnorm(1, 0, 0.3))
          c0 <- 475 + rnorm(1, 0, 10)
        } else {
          a <- -abs(0.05 + rnorm(1, 0, 0.02))
          b <- -abs(1.5 + rnorm(1, 0, 0.3))
          c0 <- 525 + rnorm(1, 0, 10)}
        y[, i] <- a * t[, i]^2 + b * t[, i] + c0} 
      
      # Sinusoïdes 
      else if (scenario == 2) { 
        if (g == 0) { 
          A <- 80; B <- 2*pi/tmax; C <- 0} 
        else if (g == 1) { 
          A <- 60; B <- 2*pi/tmax; C <- pi/4} 
        else { 
          A <- 100; B <- 2*pi/tmax; C <- pi/2}
        y[, i] <- 500 + A * sin(B * t[, i] + C)} 
      
      # Sigmoïdes
      else if (scenario == 3) {  
        if (g == 0) { 
          L <- 600; k <- 0.8; x0 <- tmax/2 }
        else if (g == 1) { 
          L <- 550; k <- 1.3; x0 <- tmax/2.5 }
        else { 
          L <- 650; k <- 1.0; x0 <- tmax/1.75 }
        y[, i] <- L / (1 + exp(-k * (t[, i] - x0)))} 
      
      # Mix quadratique + sinusoïdal
      else if (scenario == 4) { 
        if (g == 0) {
          y[, i] <- 400 + 0.5*(t[, i]^2) - 5*t[, i] + 40*sin(2*pi*t[, i]/tmax)} 
        else if (g == 1) {
          y[, i] <- 420 + 0.4*(t[, i]^2) - 4*t[, i] + 60*sin(2*pi*t[, i]/tmax + 1)} 
        else {
          y[, i] <- 450 + 0.6*(t[, i]^2) - 6*t[, i] + 20*sin(2*pi*t[, i]/tmax - 1)}} 
      
      else {stop("Scenario inconnu pour données denses")}
      
      y[, i] <- y[, i] + bruit_AR(nmes, rho, sd_bruit * mean(abs(y[, i])))
      
    } else {
      
      # Quadratique VS Quadratique
      if (scenario == 1) { 
        if (g == 0) {
          tcut <- tmax * 0.9
          a <- -abs(0.15 + rnorm(1, 0, 0.05))
          b <- -abs(1.5 + rnorm(1, 0, 1.5))
          c0 <- 500 + rnorm(1, 0, 15)
          y[, i] <- ifelse(t[, i] < tcut,
                           a * t[, i]^2 + b * t[, i] + c0,
                           a * t[, i]^2 + b * t[, i] + c0 - 2*(t[, i] - tcut)^2) + rnorm(nmes, 0, 15)
        } else {
          tcut <- tmax * 0.5
          a <- -abs(0.01 + rnorm(1, 0, 0.15))
          b <- -abs(2.5 + rnorm(1, 0, 1))
          c0 <- 500 + rnorm(1, 0, 15)
          y[, i] <- ifelse(t[, i] < tcut,
                           a * t[, i]^2 + b * t[, i] + c0,
                           a * t[, i]^2 + b * t[, i] + c0 - (t[, i] - tcut)^2) + rnorm(nmes, 0, 15)}} 
      
      # Sinusoïde VS Quadratique
      else if (scenario == 2) { 
        if (g == 0) {
          bi <- rmvnorm(1, mean = rep(0,5), sigma = diag(rep(1,5)))
          y[, i] <- 600 - (SSlogis5(t[,i], 10 + bi[1], 16 + bi[2], 13 + bi[3], 10+bi[4], 1 + bi[4]/3))^2 + rnorm(nmes, 0, 10)} 
        else {
          tcut <- tmax * 0.4
          a <- -0.2 + rnorm(1, 0, 0.125)  
          b <- -2 + rnorm(1, 0, 0.5)      
          c <- 500 + rnorm(1, 0, 10)       
          y[, i] <- ifelse(t[,i] < tcut,
                           a * t[,i]^2 + b * t[,i] + c,
                           a * t[,i]^2 + b * t[,i] + c - (t[,i] - tcut)^2) + rnorm(nmes, 0, 10)}}
        
      # Quadratique VS Sinusoïde
      else if (scenario == 3) { 
        if (g == 0) {
          tcut <- tmax * 0.75
          a <- -abs(0.05 + rnorm(1, 0, 0.2))
          b <- -abs(2 + rnorm(1, 0, 2))    
          c <- 500 + rnorm(1, 0, 15)       
          y[, i] <- ifelse(t[,i] < tcut,
                           a * t[,i]^2 + b * t[,i] + c,
                           a * t[,i]^2 + b * t[,i] + c - 3*(t[,i] - tcut)^2) + rnorm(nmes, 0, 20)} 
        else {
          bi <- rmvnorm(1, mean = rep(0,5), sigma = diag(rep(1,5)))
          y[, i] <- 550 - (SSlogis5(t[,i], 7 + bi[1], 18 + bi[2], 10 + bi[3], 5+bi[4], 1 + bi[4]/3))^2 + rnorm(nmes, 0, 15)}}
        
      # Quadratique VS Quadratique VS Sinusoïde (3 groupes)
      else if (scenario == 4) {
        if (g == 0) {
          tcut <- tmax * 0.9
          a <- -abs(0.025 + rnorm(1, 0, 0.15))
          b <- -abs(1.5 + rnorm(1, 0, 1.5))
          c0 <- 500 + rnorm(1, 0, 15)
          y[, i] <- ifelse(t[, i] < tcut,
                           a * t[, i]^2 + b * t[, i] + c0,
                           a * t[, i]^2 + b * t[, i] + c0 - 3*(t[, i] - tcut)^2) + rnorm(nmes, 0, 12.5)} 
        else if (g == 1) {
          tcut <- tmax * 0.5
          a <- -abs(0.15 + rnorm(1, 0, 0.1))
          b <- -abs(2.5 + rnorm(1, 0, 2))
          c0 <- 500 + rnorm(1, 0, 15)
          y[, i] <- ifelse(t[, i] < tcut,
                           a * t[, i]^2 + b * t[, i] + c0,
                           a * t[, i]^2 + b * t[, i] + c0 - (t[, i] - tcut)^2) + rnorm(nmes, 0, 15)} 
        else if (g == 2) {
          bi <- rmvnorm(1, mean = rep(0,5), sigma = diag(rep(0.75,5)))
          y[, i] <- 550 - (SSlogis5(t[,i], 7 + bi[1], 18 + bi[2], 10 + bi[3], 5+bi[4], 1 + bi[4]/3))^2 + rnorm(nmes, 0, 15)}}
        
      else {stop("Scenario inconnu pour data sparse")}
    }
  } 
  
  # 6 - Normalisation globale 
  y <- (y - min(y)) / (max(y) - min(y))
  
  return(list(t = t, y = y, groupes = groups, 
              n = n, nmes = nmes, scenario = scenario, n.grps = ngroups))
}
dropout <- function(data, miss.pattern, txdo, scenario){
  
  # Initialisation
  n <- data$n
  nmes <- data$nmes
  
  # MAR - Fixed
  if (miss.pattern == "fixed"){
    
    dtpast <- rbind(rep(500, n), data$y[-nmes,])
    dropout <- matrix(NA, ncol = n, nrow = nmes)

    if (scenario == 1)
      dropout <- ifelse(dtpast < 0.70 - 0.2*(txdo==0.3) + 0.1*(nmes == 7) - 0.1*(nmes == 7)*(txdo==0.3) + 0.15*(nmes == 5), 0, 1) 
    
    if (scenario == 2) 
      dropout <- ifelse(dtpast < 0.55 - 0.165*(txdo==0.3) + 0.05*(nmes == 7) - 0.05*(nmes == 7)*(txdo==0.3) + 0.15*(nmes == 5), 0, 1) 
    
    if (scenario == 3)
      dropout <- ifelse(dtpast < 0.6 - 0.3*(txdo==0.3) + 0.05*(nmes == 7) - 0.05*(nmes == 7)*(txdo==0.3) + 0.15*(nmes == 5), 0, 1)
    
    if (scenario == 4)
      dropout <- ifelse(dtpast < 0.64 - 0.3*(txdo==0.3) + 0.05*(nmes == 7) - 0.05*(nmes == 7)*(txdo==0.3) + 0.15*(nmes == 5), 0, 1)
    
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    data$y.do <- ifelse(dropout == 0, NA, data$y) 
    
    txevent <- sum(!apply(dropout, 2, prod))/n 
    meanNA <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x)))))
    meanNAdo <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x))))[which(apply(data$y.do, 2, function(x) return(sum(is.na(x))))>0)])}
  
  
  # MAR - Threshold
  if (miss.pattern == "threshMAR"){
    
    if (scenario == 1)
      prob_Yobs_linpred <- 4.5*(data$y < 0.60) - 3.25 + 1*(nmes == 5) + 0.75*(nmes == 7) - 2*(txdo == 0.3) - 0.75*(nmes == 5)*(txdo == 0.3) - 0.5*(nmes==7)*(txdo == 0.3)
    
    if (scenario == 2) 
      prob_Yobs_linpred <- 4.5*(data$y < 0.55) - 3.25 + 1*(nmes == 5) + 0.75*(nmes == 7) - 1.85*(txdo == 0.3) - 0.75*(nmes == 5)*(txdo == 0.3) - 0.5*(nmes==7)*(txdo == 0.3) 
    
    if (scenario == 3) 
      prob_Yobs_linpred <- 4*(data$y < 0.60) - 3.5 + 1*(nmes == 5) + 0.75*(nmes == 7) - 2.1*(txdo == 0.3) - 0.75*(nmes == 5)*(txdo == 0.3) - 0.5*(nmes==7)*(txdo == 0.3)
    
    if (scenario == 4) 
      prob_Yobs_linpred <- 4*(data$y < 0.64) - 3.5 + 1*(nmes == 5) + 0.75*(nmes == 7) - 2.1*(txdo == 0.3) - 0.75*(nmes == 5)*(txdo == 0.3) - 0.5*(nmes==7)*(txdo == 0.3)
    
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    prob_Yobs <- rbind(matrix(0, ncol = n, nrow = 1), prob_Yobs[-nmes,])
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    data$y.do <- ifelse(dropout == 0, NA, data$y) 
    
    data$txevent <- sum(!apply(dropout, 2, prod))/n
    data$meanNA <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x)))))
    data$meanNAdo <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x))))[which(apply(data$y.do, 2, function(x) return(sum(is.na(x))))>0)])}
  
  
  # MAR - Increasing
  if (miss.pattern == "incrMAR"){
    
    if (scenario == 1) 
      prob_Yobs_linpred <- -2.7 * data$y + 0.5*(nmes==7) + 1.5*(nmes == 5) - 1.25*(txdo == 0.3) - 1*(nmes==5)*(txdo==0.3)
    
    if (scenario == 2) 
      prob_Yobs_linpred <- -2.7 * data$y + 0.5*(nmes==7) + 1.5*(nmes == 5) - 1.25*(txdo == 0.3) - 0.5*(nmes==5)*(txdo==0.3) 
    
    if (scenario == 3) 
      prob_Yobs_linpred <- -2.8 * data$y + 0.5*(nmes==7) + 1.5*(nmes == 5) - 1.35*(txdo == 0.3) - 0.5*(nmes==5)*(txdo==0.3) 
    
    if (scenario == 4) 
      prob_Yobs_linpred <- -2.75 * data$y + 0.5*(nmes==7) + 1.5*(nmes == 5) - 1.35*(txdo == 0.3) - 0.5*(nmes==5)*(txdo==0.3) 
    
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    prob_Yobs <- rbind(matrix(0, ncol = n, nrow = 1), prob_Yobs[-nmes,])
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    data$y.do <- ifelse(dropout == 0, NA, data$y) 
    
    data$txevent <- sum(!apply(dropout, 2, prod))/n
    data$meanNA <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x)))))
    data$meanNAdo <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x))))[which(apply(data$y.do, 2, function(x) return(sum(is.na(x))))>0)])}
  
  
  # MNAR - Threshold
  if (miss.pattern == "threshMNAR"){
     
    if (scenario == 1)
      prob_Yobs_linpred <- 4*(data$y < 0.50) - 3 + 1*(nmes == 5) + 0.75*(nmes == 7) - 1.75*(txdo == 0.3) - 0.75*(nmes == 5)*(txdo == 0.3) - 0.5*(nmes==7)*(txdo == 0.3)
    
    if (scenario == 2) 
      prob_Yobs_linpred <- 4*(data$y < 0.50) - 3.25 + 1*(nmes == 5) + 0.75*(nmes == 7) - 1.75*(txdo == 0.3) - 0.75*(nmes == 5)*(txdo == 0.3) - 0.5*(nmes==7)*(txdo == 0.3) 
    
    if (scenario == 3) 
      prob_Yobs_linpred <- 4*(data$y < 0.50) - 4 + 1*(nmes == 5) + 0.75*(nmes == 7) - 1.6*(txdo == 0.3) - 0.75*(nmes == 5)*(txdo == 0.3) - 0.5*(nmes==7)*(txdo == 0.3)
    
    if (scenario == 4) 
      prob_Yobs_linpred <- 4*(data$y < 0.55) - 4.1 + 1*(nmes == 5) + 0.75*(nmes == 7) - 1.75*(txdo == 0.3) - 0.75*(nmes == 5)*(txdo == 0.3) - 0.5*(nmes==7)*(txdo == 0.3)
    
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    data$y.do <- ifelse(dropout == 0, NA, data$y) 
    
    data$txevent <- sum(!apply(dropout, 2, prod))/n
    data$meanNA <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x)))))
    data$meanNAdo <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x))))[which(apply(data$y.do, 2, function(x) return(sum(is.na(x))))>0)])}
  
  
  # MNAR - Increasing
  if (miss.pattern == "incrMNAR"){
    
    if (scenario == 1)
      prob_Yobs_linpred <- -2.75 * data$y + 0.5*(nmes==7) + 1.5*(nmes == 5) - 1.3*(txdo == 0.3) - 1*(nmes==5)*(txdo==0.3)
    
    if (scenario == 2)
      prob_Yobs_linpred <- -2.9 * data$y + 0.5*(nmes==7) + 1.5*(nmes == 5) - 1.3*(txdo == 0.3) - 0.5*(nmes==5)*(txdo==0.3)
    
    if (scenario == 3) 
      prob_Yobs_linpred <- -3.25 * data$y + 0.5*(nmes==7) + 1.5*(nmes == 5) - 1.35*(txdo == 0.3) - 0.5*(nmes==5)*(txdo==0.3) 
    
    if (scenario == 4) 
      prob_Yobs_linpred <- -3.1 * data$y + 0.5*(nmes==7) + 1.5*(nmes == 5) - 1.35*(txdo == 0.3) - 0.5*(nmes==5)*(txdo==0.3) 
    
    prob_Yobs <- exp(prob_Yobs_linpred)/(1+exp(prob_Yobs_linpred))
    dropout <- !matrix(rbinom(n = n*nmes, size = 1, prob = prob_Yobs), ncol = n, nrow = nmes)
    dropout <- apply(dropout, 2, cumprod)
    dropout[1,] <- rep(1, n)
    data$y.do <- ifelse(dropout == 0, NA, data$y) 
    
    data$txevent <- sum(!apply(dropout, 2, prod))/n
    data$meanNA <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x)))))
    data$meanNAdo <- mean(apply(data$y.do, 2, function(x) return(sum(is.na(x))))[which(apply(data$y.do, 2, function(x) return(sum(is.na(x))))>0)])}
  
  return(data)
}
plot.data <- function(data, title = "", legend = FALSE, legend.name = "", x = "", y = "", title.pos = 0) {
  
  name.x <- if (x == "") "Temps" else x
  name.y <- if (y == "") "Marqueur" else y
  legend.pos <- if (legend) "right" else "none"
  legend.name <- if (legend.name == "") "Type d'évolution" else legend.name 
  
  if ("y.do" %in% names(data)) {
    
    # Données incomplètes (NA)
    data.long <- data.frame(
      id = rep(1:ncol(data$t), each = nrow(data$t)),
      temps = as.vector(data$t),
      y.do = as.vector(data$y.do),
      groupe = factor(rep(data$groupes, each = nrow(data$t)))
    )
    
    cols <- if (length(unique(data$groupes)) == 2) {
      c("0" = "#434343", "1" = "#85C1E9")
    } else {
      c("0" = "#434343", "1" = "#85C1E9", "2" = "#FA9ECD")
    }
    
    # Derniers points observés
    data.red <- data.long %>%
      group_by(id) %>%
      filter(any(is.na(y.do))) %>%
      filter(!is.na(y.do)) %>%
      filter(row_number() == max(row_number())) %>%
      mutate(shape_point = "Dernier point observé") %>%
      ungroup()
    
    if (title == "") title <- "Données simulées"
    
    graphique <- ggplot() +
      geom_line(data = data.long, aes(x = temps, y = y.do, group = id, color = groupe), alpha = 0.75) +
      geom_point(data = data.long, aes(x = temps, y = y.do, group = id, color = groupe)) +
      geom_point(
        data = data.red,
        aes(x = temps, y = y.do, shape = shape_point),
        color = "red",
        size = 2,
        show.legend = FALSE
      ) +
      scale_color_manual(name = legend.name, values = cols) +
      scale_shape_manual(name = "", values = c("Dernier point observé" = 4)) +
      theme_minimal() +
      theme(legend.position = legend.pos, plot.title = element_text(hjust = title.pos)) +
      labs(title = title, x = name.x, y = name.y)}
  
  else { 
    
    # Données complètes
    data.long <- data.frame(
      id = rep(1:ncol(data$t), each = nrow(data$t)),
      temps = as.vector(data$t),
      y = as.vector(data$y),
      groupe = factor(rep(data$groupes, each = nrow(data$t)))
    )
    
    cols <- if (length(unique(data$groupes)) == 2) {
      c("0" = "#434343", "1" = "#85C1E9")
    } else {
      c("0" = "#434343", "1" = "#85C1E9", "2" = "#AC4EBF")
    }
    
    if (title == "") title <- "Mesures longitudinales par individu"
    
    graphique <- ggplot(data.long, aes(x = as.numeric(temps), y = y, group = id, color = groupe)) +
      geom_line(alpha = 0.2) + 
      scale_color_manual(name = legend.name, values = cols) +
      theme_minimal() +
      theme(legend.position = legend.pos, plot.title = element_text(hjust = title.pos)) +
      labs(title = title, x = name.x, y = name.y)} 
  
  if (!legend) {
    graphique <- graphique + guides(color = "none", shape = "none")
  }
  
  return(graphique)
}

# # Test de funFEM
# data.simu <- generation.data(250, 200, 20, type = "dense", scenario = 2, ngroups = 3)
# temps <- 1:data.simu$nmes
# y <- data.simu$y
# plot.data(data.simu)
# 
# library(fda)
# nbasis <- 20
# basis <- create.bspline.basis(rangeval = range(temps), nbasis = nbasis)
# 
# fdParobj <- fdPar(basis, Lfdobj = 2, lambda = 1e-4)
# smoothed <- smooth.basis(argvals = temps, y = y, fdParobj = fdParobj)
# fdobj <- smoothed$fd
# plot(fdobj)
# 
# res <- funFEM(fdobj, K = c(2,2)) 
# summary(res)
# plot(res)
