## Lianne Siegel
## Generate Data for Scenarios 1-3
## Also Generates Figure 2 for aggregation bias manuscript
## Binary outcome and binary covariate
library(mixmeta)
library(tidyverse)
library(ggpubr)
library(unix)

#######
## Generate Dataset with Binary Covariate (Aggregation Bias)
## Scenario 2
#######
#### Need to generate dataset where AB an issue
#### Start with 6 studies (common effect model)
#### no interaction between x and treatment effect within study 
#### but association between prop x and treatment effect

## settings
K <- 5
OR_t <- seq(0.5, 1, length.out = K)
logOR_t <- log(OR_t)

OR_x <- rep(1/2, K)
logOR_x <- log(OR_x)
pc <- rep(0.2,K)
n <- rep(100000, K)
n_x <- seq(0.1*n[1], 0.7*n[1], length.out = K) 
n_x

## generate data for each study
set.seed(123)
dat <- lapply(c(1:K), function(s){
  x <- c(rep(1, n_x[s]), rep(0, n[s] - n_x[s]))
  t <- rbinom(n[s], size = 1, prob = 0.5)
  logit_p <- log(pc[s]/(1-pc[s])) + logOR_t[s]*t + logOR_x[s]*x
  p <- 1/(1 + exp(-logit_p))
  y <- rbinom(n[s], size = 1, prob = p)
  dat <- data.frame(s, x, t, y)
})

dat_pooled <- do.call("rbind", dat)

#########
## "Naive" (Pooled) approach with study fixed effect
#########

## -------- fit pooled GLM ("naive" approach with study as fixed effect) ---------#
res_glm <- glm(y ~ factor(x)*t + factor(s), data = dat_pooled, family = "binomial")
summary(res_glm)


########
## Figure showing data/aggregation bias (panel b)
####### 
# separate analyses by subgroup
glm_study_neg <- lapply(dat, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ t, data = subset(d, x == 0), family = "binomial")
})

neg_coefs <- lapply(glm_study_neg, function(x){
  x$coef[2]
})

neg_vars <- lapply(glm_study_neg, function(x){
  vcov(x)[2,2]
})

glm_study_pos <- lapply(dat, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ t, data = subset(d, x == 1), family = "binomial")
})

pos_coefs <- lapply(glm_study_pos, function(x){
  x$coef[2]
})

pos_vars <- lapply(glm_study_pos, function(x){
  vcov(x)[2,2]
})

# pool separately
# neg_pooled <- mixmeta(unlist(neg_coefs), unlist(neg_vars), method = "fixed")
# pos_pooled <- mixmeta(unlist(pos_coefs), unlist(pos_vars), method = "fixed")

# make long dataset with coefs and weights
synth_ex <- as.data.frame(rbind(cbind(unlist(neg_coefs),
                                      unlist(neg_vars)),
                                cbind(unlist(pos_coefs),
                                      unlist(pos_vars))))
names(synth_ex) <- c("coefs", "vars")
synth_ex$weights <- 1/synth_ex$vars
synth_ex$x <- c(rep("Negative", K),
                rep("Positive", K))
synth_ex$study <- c(1:K, 1:K)
synth_ex$type <- "Study"

# add in pooled results
neg_results <- cbind.data.frame(summary(res_glm)$coef["t", "Estimate"], 
                                 vcov(res_glm)["t","t"],
                                 1/vcov(res_glm)["t","t"],
                                 "Negative",
                                 7,
                                 "Pooled")

pos_results <- cbind.data.frame(summary(res_glm)$coef["t", "Estimate"] + 
                                   summary(res_glm)$coef["factor(x)1:t", "Estimate"], 
                                 vcov(res_glm)["t","t"] +
                                   vcov(res_glm)["factor(x)1:t","factor(x)1:t"] + 
                                   2*vcov(res_glm)["t","factor(x)1:t"],
                                 1/(vcov(res_glm)["t","t"] +
                                      vcov(res_glm)["factor(x)1:t","factor(x)1:t"] + 
                                      2*vcov(res_glm)["t","factor(x)1:t"]),
                                 "Positive",
                                 7,
                                 "Pooled")

names(neg_results) <- names(synth_ex)
names(pos_results) <- names(synth_ex)

both_results <- bind_rows(neg_results, pos_results)
all_results <- bind_rows(synth_ex, both_results)

ab_plot1 <- ggplot(all_results, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type),
             fill = "blue") + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 20)) + 
  scale_alpha_manual(values = c(0.7, 1)) + 
  ylim(-1, 0.1) + 
  theme(legend.title = element_blank()) + 
  labs(x = "Baseline Serostatus",
       y = "Estimated logOR \n (Active/Placebo)") + 
  scale_size(range = c(1,10))
ab_plot1

#########
## Figure showing data/AB (panel a)
## Scenario 1
#########
## use same settings as (a) but change covariate dist to be equal
n2 <- rep(100000, K)

n_x2 <- rep(0.5*n2[1], K)

OR_x2 <- rep(1/2, K)
#OR_x <- seq((1/2), 1, length.out = K)
logOR_x2 <- log(OR_x2)

## generate data for each study
dat2 <- lapply(c(1:K), function(s){
  x <- c(rep(1, n_x2[s]), rep(0, n2[s] - n_x2[s]))
  t <- rbinom(n2[s], size = 1, prob = 0.5)
  logit_p <- log(pc[s]/(1-pc[s])) + logOR_t[s]*t + logOR_x2[s]*x
  p <- exp(logit_p)/(1 + exp(logit_p))
  y <- rbinom(n2[s], size = 1, prob = p)
  dat <- data.frame(s, x, t, y)
})

dat_pooled2 <- do.call("rbind", dat2)

# separate analyses by subgroup
glm_study_neg2 <- lapply(dat2, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ t, data = subset(d, x == 0), family = "binomial")
})

neg_coefs2 <- lapply(glm_study_neg2, function(x){
  x$coef[2]
})

neg_vars2 <- lapply(glm_study_neg2, function(x){
  vcov(x)[2,2]
})

glm_study_pos2 <- lapply(dat2, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ t, data = subset(d, x == 1), family = "binomial")
})

pos_coefs2 <- lapply(glm_study_pos2, function(x){
  x$coef[2]
})

pos_vars2 <- lapply(glm_study_pos2, function(x){
  vcov(x)[2,2]
})

# pool separately
# neg_pooled2 <- mixmeta(unlist(neg_coefs2), unlist(neg_vars2), method = "fixed")
# pos_pooled2 <- mixmeta(unlist(pos_coefs2), unlist(pos_vars2), method = "fixed")

# make long dataset with coefs and weights
synth_ex2 <- as.data.frame(rbind(cbind(unlist(neg_coefs2),
                                      unlist(neg_vars2)),
                                cbind(unlist(pos_coefs2),
                                      unlist(pos_vars2))))
names(synth_ex2) <- c("coefs", "vars")
synth_ex2$weights <- 1/synth_ex2$vars
synth_ex2$x <- c(rep("Negative", K),
                rep("Positive", K))
synth_ex2$study <- c(1:K, 1:K)
synth_ex2$type <- "Study"

# add in pooled results
# fit model
res_glm2 <-  glm(y ~ factor(x)*t + factor(s), data = dat_pooled2, family = "binomial")


neg_results2 <- cbind.data.frame(summary(res_glm2)$coef["t", "Estimate"], 
                                vcov(res_glm2)["t","t"],
                                1/vcov(res_glm2)["t","t"],
                                "Negative",
                                7,
                                "Pooled")

pos_results2 <- cbind.data.frame(summary(res_glm2)$coef["t", "Estimate"] + 
                                   summary(res_glm2)$coef["factor(x)1:t", "Estimate"], 
                                vcov(res_glm2)["t","t"] +
                                  vcov(res_glm2)["factor(x)1:t","factor(x)1:t"] + 
                                  2*vcov(res_glm2)["t","factor(x)1:t"],
                                1/(vcov(res_glm2)["t","t"] +
                                     vcov(res_glm2)["factor(x)1:t","factor(x)1:t"] + 
                                     2*vcov(res_glm2)["t","factor(x)1:t"]),
                                "Positive",
                                7,
                                "Pooled")

names(neg_results2) <- names(synth_ex2)
names(pos_results2) <- names(synth_ex2)

both_results2 <- bind_rows(neg_results2, pos_results2)
all_results2 <- bind_rows(synth_ex2, both_results2)

ab_plot2 <- ggplot(all_results2, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type),
             fill = "blue") + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 20)) + 
  scale_alpha_manual(values = c(0.7, 1)) + 
  ylim(-1, 0.1) + 
  theme(legend.title = element_blank()) + 
  labs(x = "Baseline Serostatus",
       y = "Estimated logOR \n (Active/Placebo)") +
  scale_size(range = c(1,10))

ab_plot2 


#########
## Figure showing data/AB (panel d - allocation ratio)
## Scenario 3
#########
## similar settings as in panel a but change allocation ratio 
## to be related to effect size
## also need heterogeneity in assoc between x and outcome
## generate data for each study
alloc_ratio <- seq(0.1, 0.7, length.out = K)

# no treatment effect
OR_t2 <- rep(1, K)
logOR_t2 <- log(OR_t2)

OR_x3 <- seq((1/2), 1, length.out = K)
logOR_x3 <- log(OR_x3)

# use same number of x for each study as in (a)
dat3 <- lapply(c(1:K), function(s){
  x <- c(rep(1, n_x2[s]), rep(0, n2[s] - n_x2[s]))
  t <- rbinom(n2[s], size = 1, prob = alloc_ratio[s])
  logit_p <- log(pc[s]/(1-pc[s])) + logOR_t2[s]*t + logOR_x3[s]*x
  p <- exp(logit_p)/(1 + exp(logit_p))
  y <- rbinom(n2[s], size = 1, prob = p)
  dat <- data.frame(s, x, t, y)
})

dat_pooled3 <- do.call("rbind", dat3)

## analyze assoc between x and outcome for trt and control separately
glm_study_ctrl <- lapply(dat3, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ x, data = subset(d, t == 0), family = "binomial")
})

ctrl_coefs <- lapply(glm_study_ctrl, function(x){
  x$coef[2]
})

ctrl_vars <- lapply(glm_study_ctrl, function(x){
  vcov(x)[2,2]
})

glm_study_trt <- lapply(dat3, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ x, data = subset(d, t == 1), family = "binomial")
})

trt_coefs <- lapply(glm_study_trt, function(x){
  x$coef[2]
})

trt_vars <- lapply(glm_study_trt, function(x){
  vcov(x)[2,2]
})

# pool separately
# ctrl_pooled <- mixmeta(unlist(ctrl_coefs), unlist(ctrl_vars), method = "fixed")
# trt_pooled <- mixmeta(unlist(trt_coefs), unlist(trt_vars), method = "fixed")

# make long dataset with coefs and weights
synth_ex3 <- as.data.frame(rbind(cbind(unlist(ctrl_coefs),
                                       unlist(ctrl_vars)),
                                 cbind(unlist(trt_coefs),
                                       unlist(trt_vars))))
names(synth_ex3) <- c("coefs", "vars")
synth_ex3$weights <- 1/synth_ex3$vars
synth_ex3$x <- c(rep("Control", K),
                 rep("Treatment", K))
synth_ex3$study <- c(1:K, 1:K)
synth_ex3$type <- "Study"

# add in pooled results
res_glm3 <-  glm(y ~ factor(x)*t + factor(s), data = dat_pooled3, family = "binomial")

ctrl_results <- cbind.data.frame(summary(res_glm3)$coef["factor(x)1", "Estimate"], 
                                 vcov(res_glm3)["factor(x)1","factor(x)1"],
                                 1/vcov(res_glm3)["factor(x)1","factor(x)1"],
                                 "Control",
                                 7,
                                 "Pooled")



trt_results <- cbind.data.frame(summary(res_glm3)$coef["factor(x)1", "Estimate"] + 
                                  summary(res_glm3)$coef["factor(x)1:t", "Estimate"] , 
                                vcov(res_glm3)["factor(x)1","factor(x)1"] + 
                                  vcov(res_glm3)["factor(x)1:t","factor(x)1:t"] + 
                                  2*vcov(res_glm3)["factor(x)1:t", "factor(x)1"],
                                1/vcov(res_glm3)["factor(x)1","factor(x)1"],
                                "Treatment",
                                7,
                                "Pooled")

names(ctrl_results) <- names(synth_ex3)
names(trt_results) <- names(synth_ex3)

both_results3 <- bind_rows(ctrl_results, trt_results)
all_results3 <- bind_rows(synth_ex3, both_results3)

ab_plot3 <- ggplot(all_results3, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type),
             fill = "blue") + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 20)) + 
  scale_alpha_manual(values = c(0.7, 1)) + 
  ylim(-1, 0.1) + 
  theme(legend.title = element_blank()) + 
  labs(x = "Treatment Group",
       y = "Estimated logOR \n (Seropositive/Seronegative)") +
  scale_size(range = c(1,10))
ab_plot3

#########
## Figure showing data/AB (panel c - allocation ratio; Scenario 1)
#########
## equal allocation ratios
alloc_ratio2 <- rep(0.5, K)


# use same number of x for each study as in (a)
# dat4 <- lapply(c(1:K), function(s){
#   x <- c(rep(1, n_x2[s]), rep(0, n2[s] - n_x2[s]))
#   t <- rbinom(n2[s], size = 1, prob = alloc_ratio2[s])
#   logit_p <- log(pc[s]/(1-pc[s])) + logOR_t2[s]*t + logOR_x2[s]*x
#   p <- exp(logit_p)/(1 + exp(logit_p))
#   y <- rbinom(n2[s], size = 1, prob = p)
#   dat <- data.frame(s, x, t, y)
# })
# 
# dat4_pooled <- do.call("rbind", dat4)
dat4 <- dat2
dat4_pooled <- dat_pooled2


## analyze assoc between x and outcome for trt and control separately
glm_study_ctrl2 <- lapply(dat4, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ factor(x), data = subset(d, t == 0), family = "binomial")
})

ctrl2_coefs <- lapply(glm_study_ctrl2, function(x){
  x$coef[2]
})

ctrl2_vars <- lapply(glm_study_ctrl2, function(x){
  vcov(x)[2,2]
})

glm_study_trt2 <- lapply(dat4, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ x, data = subset(d, t == 1), family = "binomial")
})

trt2_coefs <- lapply(glm_study_trt2, function(x){
  x$coef[2]
})

trt2_vars <- lapply(glm_study_trt2, function(x){
  vcov(x)[2,2]
})

# pool separately
# ctrl2_pooled <- mixmeta(unlist(ctrl2_coefs), unlist(ctrl2_vars), method = "fixed")
# trt2_pooled <- mixmeta(unlist(trt2_coefs), unlist(trt2_vars), method = "fixed")

# make long dataset with coefs and weights
synth_ex4 <- as.data.frame(rbind(cbind(unlist(ctrl2_coefs),
                                       unlist(ctrl2_vars)),
                                 cbind(unlist(trt2_coefs),
                                       unlist(trt2_vars))))
names(synth_ex4) <- c("coefs", "vars")
synth_ex4$weights <- 1/synth_ex4$vars
synth_ex4$x <- c(rep("Control", K),
                 rep("Treatment", K))
synth_ex4$study <- c(1:K, 1:K)
synth_ex4$type <- "Study"

# add in pooled results
#res_glm4 <-  glm(y ~ factor(x)*t + factor(s), data = dat4_pooled)
res_glm4 <- res_glm2

ctrl2_results <- cbind.data.frame(summary(res_glm4)$coef["factor(x)1", "Estimate"], 
                                 vcov(res_glm4)["factor(x)1","factor(x)1"],
                                 1/vcov(res_glm4)["factor(x)1","factor(x)1"],
                                 "Control",
                                 7,
                                 "Pooled")



trt2_results <- cbind.data.frame(summary(res_glm4)$coef["factor(x)1", "Estimate"] + 
                                  summary(res_glm4)$coef["factor(x)1:t", "Estimate"] , 
                                vcov(res_glm4)["factor(x)1","factor(x)1"] + 
                                  vcov(res_glm4)["factor(x)1:t","factor(x)1:t"] + 
                                  2*vcov(res_glm4)["factor(x)1:t", "factor(x)1"],
                                1/vcov(res_glm4)["factor(x)1","factor(x)1"],
                                "Treatment",
                                7,
                                "Pooled")
names(ctrl2_results) <- names(synth_ex4)
names(trt2_results) <- names(synth_ex4)

both_results4 <- bind_rows(ctrl2_results, trt2_results)
all_results4 <- bind_rows(synth_ex4, both_results4)

ab_plot4 <- ggplot(all_results4, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type),
             fill = "blue") + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 20)) + 
  scale_alpha_manual(values = c(0.7, 1)) + 
  theme(legend.title = element_blank()) + 
  labs(x = "Treatment Group",
       y = "Estimated logOR \n (Seropositive/Seronegative)") + 
  ylim(-1, 0.1) + 
  scale_size(range = c(1,10))

ab_plot4



#########
## Put panels together into one figure 
#########
ab_plot1_figure1 <- ab_plot1 + 
  theme(plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"))
ab_plot2 <- ab_plot2 + 
  theme(plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"))
ab_plot3 <- ab_plot3 + 
  theme(plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"))
ab_plot4 <- ab_plot4 + 
  theme(plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"))


ggarrange(ab_plot2, ab_plot1_figure1, ab_plot4, ab_plot3,
          nrow = 2,
          ncol = 2,
          labels = c("(a) Equal Covariate Distribution and Allocation Ratio (Scenario 1)",
                     "(b) Covariate Distribution Related to logOR (Scenario 2)",
                     "(c) Equal Covariate Distribution and Allocation Ratio (Scenario 1)",
                     "(d) Allocation Ratio Related to logOR (Scenario 3)"),
          hjust = 0,
          common.legend = TRUE,
          font.label = list(size = 10))

#ggsave("Manuscript Tables and Figures/Figure 2 - Scenarios_v5.eps", width = 10, height = 7, dpi = 600, device = cairo_pdf)


#########
## Figure showing data/AB 
## Scenario 4 (Supplemental: same as Scenario 1 but with true effect modification)
#########
## effect of interaction
logrOR <- log(0.8)

## generate data for each study
dat_scen4 <- lapply(c(1:K), function(s){
  x <- c(rep(1, n_x2[s]), rep(0, n2[s] - n_x2[s]))
  t <- rbinom(n2[s], size = 1, prob = 0.5)
  logit_p <- log(pc[s]/(1-pc[s])) + logOR_t[s]*t + logOR_x2[s]*x + logrOR*x*t
  p <- exp(logit_p)/(1 + exp(logit_p))
  y <- rbinom(n2[s], size = 1, prob = p)
  dat <- data.frame(s, x, t, y)
})

dat_pooled_scen4 <- do.call("rbind", dat_scen4)

# separate analyses by subgroup
glm_study_neg_scen4 <- lapply(dat_scen4, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ t, data = subset(d, x == 0), family = "binomial")
})

neg_coefs_scen4 <- lapply(glm_study_neg_scen4, function(x){
  x$coef[2]
})

neg_vars_scen4 <- lapply(glm_study_neg_scen4, function(x){
  vcov(x)[2,2]
})

glm_study_pos_scen4 <- lapply(dat_scen4, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ t, data = subset(d, x == 1), family = "binomial")
})

pos_coefs_scen4 <- lapply(glm_study_pos_scen4, function(x){
  x$coef[2]
})

pos_vars_scen4 <- lapply(glm_study_pos_scen4, function(x){
  vcov(x)[2,2]
})

# make long dataset with coefs and weights
synth_ex_scen4 <- as.data.frame(rbind(cbind(unlist(neg_coefs_scen4),
                                            unlist(neg_vars_scen4)),
                                      cbind(unlist(pos_coefs_scen4),
                                            unlist(pos_vars_scen4))))
names(synth_ex_scen4) <- c("coefs", "vars")
synth_ex_scen4$weights <- 1/synth_ex_scen4$vars
synth_ex_scen4$x <- c(rep("Negative", K),
                      rep("Positive", K))
synth_ex_scen4$study <- c(1:K, 1:K)
synth_ex_scen4$type <- "Study"

# add in pooled results
# fit model
res_glm_scen4 <-  glm(y ~ factor(x)*t + factor(s), data = dat_pooled_scen4, family = "binomial")


neg_results_scen4 <- cbind.data.frame(summary(res_glm_scen4)$coef["t", "Estimate"], 
                                      vcov(res_glm_scen4)["t","t"],
                                      1/vcov(res_glm_scen4)["t","t"],
                                      "Negative",
                                      7,
                                      "Pooled")

pos_results_scen4 <- cbind.data.frame(summary(res_glm_scen4)$coef["t", "Estimate"] + 
                                        summary(res_glm_scen4)$coef["factor(x)1:t", "Estimate"], 
                                      vcov(res_glm_scen4)["t","t"] +
                                        vcov(res_glm_scen4)["factor(x)1:t","factor(x)1:t"] + 
                                        2*vcov(res_glm_scen4)["t","factor(x)1:t"],
                                      1/(vcov(res_glm_scen4)["t","t"] +
                                           vcov(res_glm_scen4)["factor(x)1:t","factor(x)1:t"] + 
                                           2*vcov(res_glm_scen4)["t","factor(x)1:t"]),
                                      "Positive",
                                      7,
                                      "Pooled")

names(neg_results_scen4) <- names(synth_ex_scen4)
names(pos_results_scen4) <- names(synth_ex_scen4)

both_results_scen4 <- bind_rows(neg_results_scen4, pos_results_scen4)
all_results_scen4 <- bind_rows(synth_ex_scen4, both_results_scen4)

ab_plot_scen4 <- ggplot(all_results_scen4, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type),
             fill = "blue") + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 20)) + 
  scale_alpha_manual(values = c(0.7, 1)) + 
  theme(legend.title = element_blank()) + 
  ylim(-1, 0.1) + 
  labs(x = "Baseline Serostatus",
       y = "Estimated logOR \n (Active/Placebo)") +
  scale_size(range = c(1,10))

ab_plot_scen4

#ggsave("Manuscript Tables and Figures/Figure S1 - Scenario4.eps", width = 6, height = 4, dpi = 600, device = cairo_pdf)


###################
# Response to reviewer: figure showing treatment effects in each group for Scenario 3
###################

## -------- fit pooled GLM ("naive" approach with study as fixed effect) ---------#
res_glm_scen3 <- glm(y ~ factor(x)*t + factor(s), data = dat_pooled3, family = "binomial")
summary(res_glm_scen3)


########
## Figure showing data/aggregation bias (panel b)
####### 
# separate analyses by subgroup
glm_study_neg_scen3 <- lapply(dat3, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ t, data = subset(d, x == 0), family = "binomial")
})

neg_coefs_scen3 <- lapply(glm_study_neg_scen3, function(x){
  x$coef[2]
})

neg_vars_scen3 <- lapply(glm_study_neg_scen3, function(x){
  vcov(x)[2,2]
})

glm_study_pos_scen3 <- lapply(dat3, function(x){
  d <- as.data.frame(x)
  res <- glm(y ~ t, data = subset(d, x == 1), family = "binomial")
})

pos_coefs_scen3 <- lapply(glm_study_pos_scen3, function(x){
  x$coef[2]
})

pos_vars_scen3 <- lapply(glm_study_pos_scen3, function(x){
  vcov(x)[2,2]
})


# make long dataset with coefs and weights
synth_ex_scen3 <- as.data.frame(rbind(cbind(unlist(neg_coefs_scen3),
                                      unlist(neg_vars_scen3)),
                                cbind(unlist(pos_coefs_scen3),
                                      unlist(pos_vars_scen3))))
names(synth_ex_scen3) <- c("coefs", "vars")
synth_ex_scen3$weights <- 1/synth_ex_scen3$vars
synth_ex_scen3$x <- c(rep("Negative", K),
                rep("Positive", K))
synth_ex_scen3$study <- c(1:K, 1:K)
synth_ex_scen3$type <- "Study"

# add in pooled results
neg_results_scen3 <- cbind.data.frame(summary(res_glm_scen3)$coef["t", "Estimate"], 
                                vcov(res_glm_scen3)["t","t"],
                                1/vcov(res_glm_scen3)["t","t"],
                                "Negative",
                                7,
                                "Pooled")

pos_results_scen3 <- cbind.data.frame(summary(res_glm_scen3)$coef["t", "Estimate"] + 
                                  summary(res_glm_scen3)$coef["factor(x)1:t", "Estimate"], 
                                vcov(res_glm_scen3)["t","t"] +
                                  vcov(res_glm_scen3)["factor(x)1:t","factor(x)1:t"] + 
                                  2*vcov(res_glm_scen3)["t","factor(x)1:t"],
                                1/(vcov(res_glm_scen3)["t","t"] +
                                     vcov(res_glm_scen3)["factor(x)1:t","factor(x)1:t"] + 
                                     2*vcov(res_glm_scen3)["t","factor(x)1:t"]),
                                "Positive",
                                7,
                                "Pooled")

names(neg_results_scen3) <- names(synth_ex_scen3)
names(pos_results_scen3) <- names(synth_ex_scen3)

both_results_scen3 <- bind_rows(neg_results_scen3, pos_results_scen3)
all_results_scen3 <- bind_rows(synth_ex_scen3, both_results_scen3)

ab_plot_scen3 <- ggplot(all_results_scen3, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type),
             fill = "blue") + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 20)) + 
  scale_alpha_manual(values = c(0.7, 1)) + 
  ylim(-1, 0.15) + 
  theme(legend.title = element_blank()) + 
  labs(x = "Baseline Serostatus",
       y = "Estimated logOR \n (Active/Placebo)") + 
  scale_size(range = c(1,10))
ab_plot_scen3

#ggsave("Manuscript Tables and Figures/SupplementalScenario3.eps", width = 6, height = 4, dpi = 600, device = cairo_pdf)



