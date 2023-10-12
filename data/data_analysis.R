## A Semiparametric Instrumented Difference-in-Differences Approach to Policy Learning
## real data application on ALS data


library(rgenoud)

load("ALS_DATA.RData")

policy <- matrix(nrow = 7, ncol = 7)
rownames(policy) <- c("IV.t0", "IV.t1", "IPW1", "IPW2", "Wald", "MR1", "MR2")
colnames(policy) <- paste0("beta", 0:6)

set.seed(7795)

# IV using only data at time 0
glmZ <- glm(attitude ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
            family = binomial(link = "logit"),
            data = DATA[DATA$t == 0, ])
pi.Z <- DATA$attitude[DATA$t == 0] * glmZ$fitted.values + (1 - DATA$attitude[DATA$t == 0]) * (1 - glmZ$fitted.values)
  
glmA.z1 <- glm(year_edu ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
               family = binomial(link = "logit"),
               data = DATA[DATA$t == 0 & DATA$attitude == 1, ])
  
glmA.z0 <- glm(year_edu ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
               family = binomial(link = "logit"),
               data = DATA[DATA$t == 0 & DATA$attitude == 0, ])

delta.X <- predict(glmA.z1, DATA[DATA$t == 0, ], type = "response") - predict(glmA.z0, DATA[DATA$t == 0, ], type = "response")

Vhat <- function(d) {
  Ad <- as.numeric(c(cbind(1, as.matrix(DATA[DATA$t == 0, 5:10])) %*% d) > 0)
  mean(as.numeric(Ad == DATA$year_edu[DATA$t == 0]) * (2 * DATA$attitude[DATA$t == 0] - 1) * (2 * DATA$year_edu[DATA$t == 0] - 1) * DATA$wage_hour[DATA$t == 0] / pi.Z / delta.X)
}
Vopt <- genoud(Vhat, default.domains = 1000, nvars = 7,
               pop.size = 6000, max = TRUE, print.level = 0,
               wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
               P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
               P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
               starting.values = NULL, hard.generation.limit = FALSE,
               solution.tolerance = 1e-04, optim.method = "Nelder-Mead")

policy["IV.t0", ] <- Vopt$par / sqrt(sum(Vopt$par^2))


# IV using only data at time 1
glmZ <- glm(attitude ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
            family = binomial(link = "logit"),
            data = DATA[DATA$t == 1, ])
pi.Z <- DATA$attitude[DATA$t == 1] * glmZ$fitted.values + (1 - DATA$attitude[DATA$t == 1]) * (1 - glmZ$fitted.values)

glmA.z1 <- glm(year_edu ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
               family = binomial(link = "logit"),
               data = DATA[DATA$t == 1 & DATA$attitude == 1, ])

glmA.z0 <- glm(year_edu ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
               family = binomial(link = "logit"),
               data = DATA[DATA$t == 1 & DATA$attitude == 0, ])

delta.X <- predict(glmA.z1, DATA[DATA$t == 1, ], type = "response") - predict(glmA.z0, DATA[DATA$t == 1, ], type = "response")

Vhat <- function(d) {
  Ad <- as.numeric(c(cbind(1, as.matrix(DATA[DATA$t == 1, 5:10])) %*% d) > 0)
  mean(as.numeric(Ad == DATA$year_edu[DATA$t == 1]) * (2 * DATA$attitude[DATA$t == 1] - 1) * (2 * DATA$year_edu[DATA$t == 1] - 1) * DATA$wage_hour[DATA$t == 1] / pi.Z / delta.X)
}
Vopt <- genoud(Vhat, default.domains = 1000, nvars = 7,
               pop.size = 6000, max = TRUE, print.level = 0,
               wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
               P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
               P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
               starting.values = NULL, hard.generation.limit = FALSE,
               solution.tolerance = 1e-04, optim.method = "Nelder-Mead")

policy["IV.t1", ] <- Vopt$par / sqrt(sum(Vopt$par^2))


#########################################
## Instrumented Difference-in-Differences

# learn nuisance parameters
glmY11 <- glm(wage_hour ~ born_australia + married + uni_mem + gov_emp + age + year_expe,
              data = DATA[DATA$t == 1 & DATA$attitude == 1, ])
muY11 <- predict(glmY11, DATA)

glmY01 <- glm(wage_hour ~ born_australia + married + uni_mem + gov_emp + age + year_expe,
              data = DATA[DATA$t == 0 & DATA$attitude == 1, ])
muY01 <- predict(glmY01, DATA)

glmY10 <- glm(wage_hour ~ born_australia + married + uni_mem + gov_emp + age + year_expe,
              data = DATA[DATA$t == 1 & DATA$attitude == 0, ])
muY10 <- predict(glmY10, DATA)

glmY00 <- glm(wage_hour ~ born_australia + married + uni_mem + gov_emp + age + year_expe,
              data = DATA[DATA$t == 0 & DATA$attitude == 0, ])
muY00 <- predict(glmY00, DATA)
  

glmA11 <- glm(year_edu ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
              family = binomial(link = "logit"),
              data = DATA[DATA$t == 1 & DATA$attitude == 1, ])
pi.A11 <- predict(glmA11, DATA, type = "response")

glmA01 <- glm(year_edu ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
              family = binomial(link = "logit"),
              data = DATA[DATA$t == 0 & DATA$attitude == 1, ])
pi.A01 <- predict(glmA01, DATA, type = "response")

glmA10 <- glm(year_edu ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
              family = binomial(link = "logit"),
              data = DATA[DATA$t == 1 & DATA$attitude == 0, ])
pi.A10 <- predict(glmA10, DATA, type = "response")

glmA00 <- glm(year_edu ~ born_australia + married + uni_mem + gov_emp + age + year_expe, 
              family = binomial(link = "logit"),
              data = DATA[DATA$t == 0 & DATA$attitude == 0, ])
pi.A00 <- predict(glmA00, DATA, type = "response")
  
glmZ <- glm(attitude ~ born_australia + married + uni_mem + gov_emp + age + year_expe,
            family = binomial(link = "logit"),
            data = DATA)
pi.TZ <- (DATA$attitude * glmZ$fitted.values + (1 - DATA$attitude) * (1 - glmZ$fitted.values)) * (DATA$t * mean(DATA$t) + (1 - DATA$t) * (1 - mean(DATA$t)))


# IPW1
Vhat <- function(d) {
  Ad <- as.numeric(c(cbind(1, as.matrix(DATA[, 5:10])) %*% d) > 0)
  mean(as.numeric(Ad == DATA$year_edu) * (2 * DATA$attitude - 1) * (2 * DATA$t - 1) * (2 * DATA$year_edu - 1) * DATA$wage_hour / pi.TZ / (pi.A11 - pi.A01 - pi.A10 + pi.A00))
}
Vopt <- genoud(Vhat, default.domains = 1000, nvars = 7,
               pop.size = 6000, max = TRUE, print.level = 0,
               wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
               P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
               P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
               starting.values = NULL, hard.generation.limit = FALSE,
               solution.tolerance = 1e-04, optim.method = "Nelder-Mead")

policy["IPW1", ] <- Vopt$par / sqrt(sum(Vopt$par^2))


# IPW2
Vhat <- function(d) {
  Ad <- as.numeric(c(cbind(1, as.matrix(DATA[, 5:10])) %*% d) > 0)
  mean(as.numeric(Ad == DATA$attitude) * (2 * DATA$t - 1) * DATA$wage_hour / pi.TZ / (pi.A11 - pi.A01 - pi.A10 + pi.A00))
}
Vopt <- genoud(Vhat, default.domains = 1000, nvars = 7,
               pop.size = 6000, max = TRUE, print.level = 0,
               wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
               P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
               P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
               starting.values = NULL, hard.generation.limit = FALSE,
               solution.tolerance = 1e-04, optim.method = "Nelder-Mead")

policy["IPW2", ] <- Vopt$par / sqrt(sum(Vopt$par^2))


# Wald
Vhat <- function(d) {
  Ad <- as.numeric(c(cbind(1, as.matrix(DATA[, 5:10])) %*% d) > 0)
  mean(Ad * (muY11 - muY01 - muY10 + muY00) / (pi.A11 - pi.A01 - pi.A10 + pi.A00))
}
Vopt <- genoud(Vhat, default.domains = 1000, nvars = 7,
               pop.size = 6000, max = TRUE, print.level = 0,
               wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
               P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
               P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
               starting.values = NULL, hard.generation.limit = FALSE,
               solution.tolerance = 1e-04, optim.method = "Nelder-Mead")

policy["Wald", ] <- Vopt$par / sqrt(sum(Vopt$par^2))


# MR1
delta.Y <- muY11 - muY01 - muY10 + muY00
delta.A <- pi.A11 - pi.A01 - pi.A10 + pi.A00

muY.TZ <- muY00 + DATA$attitude * (muY01 - muY00) + DATA$t * (muY10 - muY00) + DATA$t * DATA$attitude * (muY11 - muY01 - muY10 + muY00)
muA.TZ <- pi.A00 + DATA$attitude * (pi.A01 - pi.A00) + DATA$t * (pi.A10 - pi.A00) + DATA$t * DATA$attitude * (pi.A11 - pi.A01 - pi.A10 + pi.A00)

Vhat <- function(d) {
  Ad <- as.numeric(c(cbind(1, as.matrix(DATA[, 5:10])) %*% d) > 0)
  mean(Ad * (delta.Y / delta.A + (2 * DATA$attitude - 1) * (2 * DATA$t - 1) / pi.TZ / delta.A * (DATA$wage_hour - muY.TZ - delta.Y / delta.A * (DATA$year_edu - muA.TZ))))
}
Vopt <- genoud(Vhat, default.domains = 1000, nvars = 7,
               pop.size = 6000, max = TRUE, print.level = 0,
               wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
               P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
               P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
               starting.values = NULL, hard.generation.limit = FALSE,
               solution.tolerance = 1e-04, optim.method = "Nelder-Mead")

policy["MR1", ] <- Vopt$par / sqrt(sum(Vopt$par^2))


# MR2
Vhat <- function(d) {
  Ad <- as.numeric(c(cbind(1, as.matrix(DATA[, 5:10])) %*% d) > 0)
  mean(as.numeric(Ad == DATA$attitude) * ((2 * DATA$attitude - 1) * delta.Y / delta.A + (2 * DATA$t - 1) / pi.TZ / delta.A * (DATA$wage_hour - muY.TZ - delta.Y / delta.A * (DATA$year_edu - muA.TZ))))
}
Vopt <- genoud(Vhat, default.domains = 1000, nvars = 7,
               pop.size = 6000, max = TRUE, print.level = 0,
               wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
               P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
               P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
               starting.values = NULL, hard.generation.limit = FALSE,
               solution.tolerance = 1e-04, optim.method = "Nelder-Mead")

policy["MR2", ] <- Vopt$par / sqrt(sum(Vopt$par^2))


policy

