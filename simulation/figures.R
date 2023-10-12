library(ggplot2)

estimator.names <- c("IV.t0", "IV.t1", "IPW1", "IPW2", "Wald", "MR1", "MR2")

#############
## main paper

load("PCD_para_main.RData")
load("PCD_ml_main.RData")

data <- data.frame(PCD = c(c(PCD.para), c(PCD.ml)),
                   Estimators = rep(rep(estimator.names, each = nrow(PCD.para)), times = 2),
                   Models = rep(c("Parametric", "Machine Learning"), each = length(c(PCD.para))))

data$Estimators <- factor(data$Estimators, levels = estimator.names)
data$Models <- factor(data$Models, levels = c("Parametric", "Machine Learning"))
p <- ggplot(data, aes(Estimators, PCD, fill = Estimators)) +
  geom_boxplot() + 
  facet_grid(. ~ Models) +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x=NULL) + 
  guides(fill = guide_legend(nrow = 1))
p

ggsave(plot = p,
       filename = "main.pdf",
       device = "pdf",
       width = 7.5,
       height = 4)


##############
## IV strength

# parametric models
load("PCD_para_weak.RData")
load("PCD_para_strong.RData")

data <- data.frame(PCD = c(c(PCD.para.weak), c(PCD.para.strong)),
                   Estimators = rep(rep(estimator.names, each = nrow(PCD.para.weak)), times = 2),
                   Strength = rep(c("Weak", "Strong"), each = length(c(PCD.para.weak))))

data$Estimators <- factor(data$Estimators, levels = estimator.names)
data$Strength <- factor(data$Strength, levels = c("Weak", "Strong"))
p <- ggplot(data, aes(Estimators, PCD, fill = Estimators)) +
  geom_boxplot() + 
  facet_grid(. ~ Strength) +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x=NULL) + 
  guides(fill = guide_legend(nrow = 1))
p

ggsave(plot = p,
       filename = "strength_para.pdf",
       device = "pdf",
       width = 7.5,
       height = 4)

# machine learning
load("PCD_ml_weak.RData")
load("PCD_ml_strong.RData")

data <- data.frame(PCD = c(c(PCD.ml.weak), c(PCD.ml.strong)),
                   Estimators = rep(rep(estimator.names, each = nrow(PCD.ml.weak)), times = 2),
                   Strength = rep(c("Weak", "Strong"), each = length(c(PCD.ml.weak))))

data$Estimators <- factor(data$Estimators, levels = estimator.names)
data$Strength <- factor(data$Strength, levels = c("Weak", "Strong"))
p <- ggplot(data, aes(Estimators, PCD, fill = Estimators)) +
  geom_boxplot() + 
  facet_grid(. ~ Strength) +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x=NULL) + 
  guides(fill = guide_legend(nrow = 1))
p

ggsave(plot = p,
       filename = "strength_ml.pdf",
       device = "pdf",
       width = 7.5,
       height = 4)

##############
## sample size

# parametric models
load("PCD_para_2500.RData")
load("PCD_para_1e4.RData")

data <- data.frame(PCD = c(c(PCD.para2500), c(PCD.para1e4)),
                   Estimators = rep(rep(estimator.names, each = nrow(PCD.para2500)), times = 2),
                   Sample.size = rep(c("n = 2500", "n = 10000"), each = length(c(PCD.para2500))))

data$Estimators <- factor(data$Estimators, levels = estimator.names)
data$Sample.size <- factor(data$Sample.size, levels = c("n = 2500", "n = 10000"))
p <- ggplot(data, aes(Estimators, PCD, fill = Estimators)) +
  geom_boxplot() + 
  facet_grid(. ~ Sample.size) +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x=NULL) + 
  guides(fill = guide_legend(nrow = 1))
p

ggsave(plot = p,
       filename = "sample_size_para.pdf",
       device = "pdf",
       width = 7.5,
       height = 4)


# machine learning
load("PCD_ml_5000.RData")
load("PCD_ml_2e4.RData")

data <- data.frame(PCD = c(c(PCD.ml5000), c(PCD.ml2e4)),
                   Estimators = rep(rep(estimator.names, each = nrow(PCD.ml5000)), times = 2),
                   Sample.size = rep(c("n = 5000", "n = 20000"), each = length(c(PCD.ml5000))))

data$Estimators <- factor(data$Estimators, levels = estimator.names)
data$Sample.size <- factor(data$Sample.size, levels = c("n = 5000", "n = 20000"))
p <- ggplot(data, aes(Estimators, PCD, fill = Estimators)) +
  geom_boxplot() + 
  facet_grid(. ~ Sample.size) +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x=NULL) + 
  guides(fill = guide_legend(nrow = 1))
p

ggsave(plot = p,
       filename = "sample_size_ml.pdf",
       device = "pdf",
       width = 7.5,
       height = 4)


