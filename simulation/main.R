library(bridgedist)
library(grf)
library(rgenoud)


simData <- function(n = 1000) {
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Z <- rbinom(n, size = 1, prob = 0.5)
  T <- rbinom(n, size = 1, prob = 0.5)
  
  U0 <- rbridge(n)
  U1 <- rbridge(n)
  
  A0 <- rbinom(n, size = 1, prob = plogis(2 - 7 * Z + 0.2 * U0 + 2 * X1))
  A1 <- rbinom(n, size = 1, prob = plogis(-1.5 + 5 * Z - 0.15 * U1 + 1.5 * X2))
  
  Y0 <- rnorm(n, mean = 200 + 10 * ((-0.5 + 1.5 * X1 + 2 * X2) * A0 + 0.5 * U0 + 2 * Z + 1.5 * X1 + 2 * X2), sd = 1)
  Y1 <- rnorm(n, mean = 240 + 10 * ((-0.5 + 1.5 * X1 + 2 * X2) * A1 + 0.5 * U1 + 2 * Z + 2 * X1 + 1.5 * X2), sd = 1)
  
  A <- T * A1 + (1 - T) * A0
  Y <- T * Y1 + (1 - T) * Y0
  
  out <- list(X1 = X1, X2 = X2,
              Z = Z, T = T,
              A = A, Y = Y)
  return(out)
}


learn_policy <- function(data,
                         model = "LR",
                         V = NULL) {
  X1 <- data$X1
  X2 <- data$X2
  Z <- data$Z
  T <- data$T
  A <- data$A
  Y <- data$Y
  
  policy <- matrix(nrow = 7, ncol = 3)
  rownames(policy) <- c("IV.t0", "IV.t1", "IPW1", "IPW2", "Wald", "MR1", "MR2")
  colnames(policy) <- paste0("beta", 0:2)
  
  PCD <- numeric(7)
  names(PCD) <- c("IV.t0", "IV.t1", "IPW1", "IPW2", "Wald", "MR1", "MR2")
  
  X1.eval <- rnorm(1e6)
  X2.eval <- rnorm(1e6)
  compute.PCD <- function(d) {
    Ad <- as.numeric(c(cbind(1, X1.eval, X2.eval) %*% d) > 0)
    Aopt <- as.numeric(-0.5 + 1.5 * X1.eval + 2 * X2.eval > 0)
    PCD <- mean(Ad == Aopt)
    return(PCD)
  }
  
  #####################################################################
  ## standard IV methods without considering the longitudinal structure
  ## Cui, Yifan, and Eric Tchetgen Tchetgen. "A semiparametric instrumental variable approach to optimal treatment regimes under endogeneity." Journal of the American Statistical Association 116.533 (2021): 162-173.
  
  # IV using only data at time 0
  X1.t0 <- X1[T == 0]
  X2.t0 <- X2[T == 0]
  Z.t0 <- Z[T == 0]
  A.t0 <- A[T == 0]
  Y.t0 <- Y[T == 0]
  
  if (model == "LR") {
    glmZ <- glm(Z.t0 ~ X1.t0 + X2.t0, family = binomial(link = "logit"))
    pi.Z <- Z.t0 * glmZ$fitted.values + (1 - Z.t0) * (1 - glmZ$fitted.values)
    
    A.t0.z1 <- A.t0[Z.t0 == 1]
    X1.t0.z1 <- X1.t0[Z.t0 == 1]
    X2.t0.z1 <- X2.t0[Z.t0 == 1]
    glmA.z1 <- glm(A.t0.z1 ~ X1.t0.z1 + X2.t0.z1, family = binomial(link = "logit"))
    
    A.t0.z0 <- A.t0[Z.t0 == 0]
    X1.t0.z0 <- X1.t0[Z.t0 == 0]
    X2.t0.z0 <- X2.t0[Z.t0 == 0]
    glmA.z0 <- glm(A.t0.z0 ~ X1.t0.z0 + X2.t0.z0, family = binomial(link = "logit"))
    
    delta.X <- predict(glmA.z1, data.frame(X1.t0.z1 = X1.t0, X2.t0.z1 = X2.t0), type = "response") - predict(glmA.z0, data.frame(X1.t0.z0 = X1.t0, X2.t0.z0 = X2.t0), type = "response")
  } else if (model == "RF") {
    grfZ <- probability_forest(X = cbind(X1.t0, X2.t0), Y = as.factor(Z.t0))
    pi.Z <- Z.t0 * grfZ$predictions[, 2] + (1 - Z.t0) * grfZ$predictions[, 1]
    
    grfA <- probability_forest(X = cbind(X1.t0, X2.t0, Z.t0), Y = as.factor(A.t0))
    delta.X <- predict(grfA, cbind(X1.t0, X2.t0, 1))$predictions[, "1"] - predict(grfA, cbind(X1.t0, X2.t0, 0))$predictions[, "1"]
  } else {
    stop("model must be set to LR or RF")
  }
  
  Vhat <- function(d) {
    Ad <- as.numeric(c(cbind(1, X1.t0, X2.t0) %*% d) > 0)
    mean(as.numeric(Ad == A.t0) * (2 * Z.t0 - 1) * (2 * A.t0 - 1) * Y.t0 / pi.Z / delta.X)
  }
  Vopt <- genoud(Vhat, default.domains = 1000, nvars = 3,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  
  policy["IV.t0", ] <- Vopt$par
  PCD["IV.t0"] <- compute.PCD(Vopt$par)
  
  
  # IV using only data at time 1
  X1.t1 <- X1[T == 1]
  X2.t1 <- X2[T == 1]
  Z.t1 <- Z[T == 1]
  A.t1 <- A[T == 1]
  Y.t1 <- Y[T == 1]
  
  if (model == "LR") {
    glmZ <- glm(Z.t1 ~ X1.t1 + X2.t1, family = binomial(link = "logit"))
    pi.Z <- Z.t1 * glmZ$fitted.values + (1 - Z.t1) * (1 - glmZ$fitted.values)
    
    A.t1.z1 <- A.t1[Z.t1 == 1]
    X1.t1.z1 <- X1.t1[Z.t1 == 1]
    X2.t1.z1 <- X2.t1[Z.t1 == 1]
    glmA.z1 <- glm(A.t1.z1 ~ X1.t1.z1 + X2.t1.z1, family = binomial(link = "logit"))
    
    A.t1.z0 <- A.t1[Z.t1 == 0]
    X1.t1.z0 <- X1.t1[Z.t1 == 0]
    X2.t1.z0 <- X2.t1[Z.t1 == 0]
    glmA.z0 <- glm(A.t1.z0 ~ X1.t1.z0 + X2.t1.z0, family = binomial(link = "logit"))
    
    delta.X <- predict(glmA.z1, data.frame(X1.t1.z1 = X1.t1, X2.t1.z1 = X2.t1), type = "response") - predict(glmA.z0, data.frame(X1.t1.z0 = X1.t1, X2.t1.z0 = X2.t1), type = "response")
  } else if (model == "RF") {
    grfZ <- probability_forest(X = cbind(X1.t1, X2.t1), Y = as.factor(Z.t1))
    pi.Z <- Z.t1 * grfZ$predictions[, 2] + (1 - Z.t1) * grfZ$predictions[, 1]
    
    grfA <- probability_forest(X = cbind(X1.t1, X2.t1, Z.t1), Y = as.factor(A.t1))
    delta.X <- predict(grfA, cbind(X1.t1, X2.t1, 1))$predictions[, "1"] - predict(grfA, cbind(X1.t1, X2.t1, 0))$predictions[, "1"]
  } else {
    stop("model must be set to LR or RF")
  }
  
  Vhat <- function(d) {
    Ad <- as.numeric(c(cbind(1, X1.t1, X2.t1) %*% d) > 0)
    mean(as.numeric(Ad == A.t1) * (2 * Z.t1 - 1) * (2 * A.t1 - 1) * Y.t1 / pi.Z / delta.X)
  }
  Vopt <- genoud(Vhat, default.domains = 1000, nvars = 3,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  
  policy["IV.t1", ] <- Vopt$par
  PCD["IV.t1"] <- compute.PCD(Vopt$par)
  
  
  #########################################
  ## Instrumented Difference-in-Differences
  
  # learn nuisance parameters
  if (model == "LR") {
    Y11 <- Y[T == 1 & Z == 1]
    X1.11 <- X1[T == 1 & Z == 1]
    X2.11 <- X2[T == 1 & Z == 1]
    glmY11 <- glm(Y11 ~ X1.11 + X2.11)
    muY11 <- predict(glmY11, data.frame(X1.11 = X1, X2.11 = X2))
    
    Y01 <- Y[T == 0 & Z == 1]
    X1.01 <- X1[T == 0 & Z == 1]
    X2.01 <- X2[T == 0 & Z == 1]
    glmY01 <- glm(Y01 ~ X1.01 + X2.01)
    muY01 <- predict(glmY01, data.frame(X1.01 = X1, X2.01 = X2))
    
    Y10 <- Y[T == 1 & Z == 0]
    X1.10 <- X1[T == 1 & Z == 0]
    X2.10 <- X2[T == 1 & Z == 0]
    glmY10 <- glm(Y10 ~ X1.10 + X2.10)
    muY10 <- predict(glmY10, data.frame(X1.10 = X1, X2.10 = X2))
    
    Y00 <- Y[T == 0 & Z == 0]
    X1.00 <- X1[T == 0 & Z == 0]
    X2.00 <- X2[T == 0 & Z == 0]
    glmY00 <- glm(Y00 ~ X1.00 + X2.00)
    muY00 <- predict(glmY00, data.frame(X1.00 = X1, X2.00 = X2))
    
    A11 <- A[T == 1 & Z == 1]
    glmA11 <- glm(A11 ~ X1.11 + X2.11, family = binomial(link = "logit"))
    pi.A11 <- predict(glmA11, data.frame(X1.11 = X1, X2.11 = X2), type = "response")
    
    A01 <- A[T == 0 & Z == 1]
    glmA01 <- glm(A01 ~ X1.01 + X2.01, family = binomial(link = "logit"))
    pi.A01 <- predict(glmA01, data.frame(X1.01 = X1, X2.01 = X2), type = "response")
    
    A10 <- A[T == 1 & Z == 0]
    glmA10 <- glm(A10 ~ X1.10 + X2.10, family = binomial(link = "logit"))
    pi.A10 <- predict(glmA10, data.frame(X1.10 = X1, X2.10 = X2), type = "response")
    
    A00 <- A[T == 0 & Z == 0]
    glmA00 <- glm(A00 ~ X1.00 + X2.00, family = binomial(link = "logit"))
    pi.A00 <- predict(glmA00, data.frame(X1.00 = X1, X2.00 = X2), type = "response")
    
    glmZ <- glm(Z ~ X1 + X2, family = binomial(link = "logit"))
    pi.TZ <- (Z * glmZ$fitted.values + (1 - Z) * (1 - glmZ$fitted.values)) / 2
  } else if (model == "RF") {
    grfY <- regression_forest(X = cbind(X1, X2, T, Z), Y = Y)
    muY11 <- predict(grfY, cbind(X1, X2, 1, 1))$predictions
    muY10 <- predict(grfY, cbind(X1, X2, 1, 0))$predictions
    muY01 <- predict(grfY, cbind(X1, X2, 0, 1))$predictions
    muY00 <- predict(grfY, cbind(X1, X2, 0, 0))$predictions
    
    grfA <- probability_forest(X = cbind(X1, X2, T, Z), Y = as.factor(A))
    pi.A11 <- predict(grfA, cbind(X1, X2, 1, 1))$predictions[, "1"]
    pi.A10 <- predict(grfA, cbind(X1, X2, 1, 0))$predictions[, "1"]
    pi.A01 <- predict(grfA, cbind(X1, X2, 0, 1))$predictions[, "1"]
    pi.A00 <- predict(grfA, cbind(X1, X2, 0, 0))$predictions[, "1"]
    
    grfZ <- probability_forest(X = cbind(X1, X2), Y = as.factor(Z))
    pi.TZ <- (Z * grfZ$predictions[, 2] + (1 - Z) * grfZ$predictions[, 1]) / 2
  } else {
    stop("model must be set to LR or RF")
  }
  
  
  # IPW1
  Vhat <- function(d) {
    Ad <- as.numeric(c(cbind(1, X1, X2) %*% d) > 0)
    mean(as.numeric(Ad == A) * (2 * Z - 1) * (2 * T - 1) * (2 * A - 1) * Y / pi.TZ / (pi.A11 - pi.A01 - pi.A10 + pi.A00))
  }
  Vopt <- genoud(Vhat, default.domains = 1000, nvars = 3,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  
  policy["IPW1", ] <- Vopt$par
  PCD["IPW1"] <- compute.PCD(Vopt$par)
  
  
  # IPW2
  Vhat <- function(d) {
    Ad <- as.numeric(c(cbind(1, X1, X2) %*% d) > 0)
    mean(as.numeric(Ad == Z) * (2 * T - 1) * Y / pi.TZ / (pi.A11 - pi.A01 - pi.A10 + pi.A00))
  }
  Vopt <- genoud(Vhat, default.domains = 1000, nvars = 3,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  
  policy["IPW2", ] <- Vopt$par
  PCD["IPW2"] <- compute.PCD(Vopt$par)
  
  
  # Wald
  Vhat <- function(d) {
    Ad <- as.numeric(c(cbind(1, X1, X2) %*% d) > 0)
    mean(Ad * (muY11 - muY01 - muY10 + muY00) / (pi.A11 - pi.A01 - pi.A10 + pi.A00))
  }
  Vopt <- genoud(Vhat, default.domains = 1000, nvars = 3,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  
  policy["Wald", ] <- Vopt$par
  PCD["Wald"] <- compute.PCD(Vopt$par)
  
  
  # cross-fitted nuisance parameters
  if (model == "RF") {
    crossfit <- sample(1:V, length(Y), replace = TRUE)
    muY11 <- rep(NA, length(Y))
    muY10 <- rep(NA, length(Y))
    muY01 <- rep(NA, length(Y))
    muY00 <- rep(NA, length(Y))
    pi.A11 <- rep(NA, length(Y))
    pi.A10 <- rep(NA, length(Y))
    pi.A01 <- rep(NA, length(Y))
    pi.A00 <- rep(NA, length(Y))
    pi.TZ <- rep(NA, length(Y))
    
    XTZ <- cbind(X1, X2, T, Z)
    X11 <- cbind(X1, X2, T = 1, Z = 1)
    X10 <- cbind(X1, X2, T = 1, Z = 0)
    X01 <- cbind(X1, X2, T = 0, Z = 1)
    X00 <- cbind(X1, X2, T = 0, Z = 0)
    X <- cbind(X1, X2)
    
    for (fold in 1:V) {
      grfY <- regression_forest(X = XTZ[crossfit != fold, ], Y = Y[crossfit != fold])
      muY11[crossfit == fold] <- predict(grfY, X11[crossfit == fold, ])$predictions
      muY10[crossfit == fold] <- predict(grfY, X10[crossfit == fold, ])$predictions
      muY01[crossfit == fold] <- predict(grfY, X01[crossfit == fold, ])$predictions
      muY00[crossfit == fold] <- predict(grfY, X00[crossfit == fold, ])$predictions
      
      grfA <- probability_forest(X = XTZ[crossfit != fold, ], Y = as.factor(A[crossfit != fold]))
      pi.A11[crossfit == fold] <- predict(grfA, X11[crossfit == fold, ])$predictions[, "1"]
      pi.A10[crossfit == fold] <- predict(grfA, X10[crossfit == fold, ])$predictions[, "1"]
      pi.A01[crossfit == fold] <- predict(grfA, X01[crossfit == fold, ])$predictions[, "1"]
      pi.A00[crossfit == fold] <- predict(grfA, X00[crossfit == fold, ])$predictions[, "1"]
      
      grfZ <- probability_forest(X = X[crossfit != fold, ], Y = as.factor(Z[crossfit != fold]))
      pi.TZ[crossfit == fold] <- (Z[crossfit == fold] * predict(grfZ, X[crossfit == fold, ])$predictions[, "1"] + (1 - Z[crossfit == fold]) * predict(grfZ, X[crossfit == fold, ])$predictions[, "0"]) / 2
    }
  }
  
  # MR1
  delta.Y <- muY11 - muY01 - muY10 + muY00
  delta.A <- pi.A11 - pi.A01 - pi.A10 + pi.A00
  
  muY.TZ <- muY00 + Z * (muY01 - muY00) + T * (muY10 - muY00) + T * Z * (muY11 - muY01 - muY10 + muY00)
  muA.TZ <- pi.A00 + Z * (pi.A01 - pi.A00) + T * (pi.A10 - pi.A00) + T * Z * (pi.A11 - pi.A01 - pi.A10 + pi.A00)
  
  Vhat <- function(d) {
    Ad <- as.numeric(c(cbind(1, X1, X2) %*% d) > 0)
    mean(Ad * (delta.Y / delta.A + (2 * Z - 1) * (2 * T - 1) / pi.TZ / delta.A * (Y - muY.TZ - delta.Y / delta.A * (A - muA.TZ))))
  }
  Vopt <- genoud(Vhat, default.domains = 1000, nvars = 3,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  
  policy["MR1", ] <- Vopt$par
  PCD["MR1"] <- compute.PCD(Vopt$par)
  
  
  # MR2
  Vhat <- function(d) {
    Ad <- as.numeric(c(cbind(1, X1, X2) %*% d) > 0)
    mean(as.numeric(Ad == Z) * ((2 * Z - 1) * delta.Y / delta.A + (2 * T - 1) / pi.TZ / delta.A * (Y - muY.TZ - delta.Y / delta.A * (A - muA.TZ))))
  }
  Vopt <- genoud(Vhat, default.domains = 1000, nvars = 3,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  
  policy["MR2", ] <- Vopt$par
  PCD["MR2"] <- compute.PCD(Vopt$par)
  
  out <- list(policy = policy,
              PCD = PCD)
  return(out)
}


#######################
## Parametric models ##
#######################

N <- 5000
num_boot <- 500

PCD.para <- matrix(nrow = num_boot, ncol = 7)
colnames(PCD.para) <- c("IV.t0", "IV.t1", "IPW1", "IPW2", "Wald", "MR1", "MR2")

set.seed(9090)
for (i in 1:num_boot) {
  cat(i, as.character(Sys.time()), "\n")
  data <- simData(N)
  res <- learn_policy(data, model = "LR")
  PCD.para[i, ] <- res$PCD
}

save(PCD.para, file = "PCD_para_main.RData")


######################
## Machine learning ##
######################

N <- 1e4
num_boot <- 500

PCD.ml <- matrix(nrow = num_boot, ncol = 7)
colnames(PCD.ml) <- c("IV.t0", "IV.t1", "IPW1", "IPW2", "Wald", "MR1", "MR2")

set.seed(6580)
for (i in 1:num_boot) {
  cat(i, as.character(Sys.time()), "\n")
  data <- simData(N)
  res <- learn_policy(data, model = "RF", V = 4)
  PCD.ml[i, ] <- res$PCD
}

save(PCD.ml, file = "PCD_ml_main.RData")


