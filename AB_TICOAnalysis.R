### Lianne Siegel
### Aggregation Bias Manuscript
### TICO/ACTIV-3 results using publicly available data

library(haven)
library(survival)
library(cmprsk)
library(broom)
library(coxphf)
library(tidyverse)
library(crrSC)
library(meta)
library(mixmeta)

### Import data
tico <- read_sas("./TICO Data/ticopublic.sas7bdat")

### Variables we will need: 
###   - active (active vs. placebo)
###   - Treatment assigned (treatment)
###      - "Bamlanivimab", "Sotrovimab", "Amubarvimab-romlusevimab", "Tixagevimab-cilgavimab", "Ensovibep"
###   - gs_pos_0 (GenScript positive on Day 0/baseline)

table(tico$active, useNA = "ifany")
table(tico$treatment, useNA = "ifany")

# 82 NAs for baseline GS
table(tico$gs_pos_0, useNA = "ifany")

### Remove Pfizer data and participants without baseline GS
tico <- filter(tico, treatment != "F1",
               treatment != "F2",
               is.na(gs_pos_0) == FALSE)
table(tico$treatment)
table(tico$gs_pos_0, useNA = "ifany")
nrow(tico)

### Make list of data from each study 
tico <- tico %>%
  mutate(trial = fct_collapse(treatment,
                              "Bamlanivimab" = c("A1", "A2"),
                              "Sotrovimab" = c("B1", "B2"),
                              "Amubarvimab-romlusevimab" = c("C1", "C2"),
                              "Tixagevimab-cilgavimab" = c("D1", "D2"),
                              "Ensovibep" = c("E1", "E2")))
table(tico$trial)

df_list <- split(tico, tico$trial)
table(tico$trial, tico$gs_pos_0)
table(tico$trial, tico$recovery)
summary(tico$t2recovery)

###################################
## Overall results
####################################
sr_list <- lapply(df_list, FUN = function(x){
  crr(ftime = x$t2recovery, fstatus = x$recovery,
      cov1 = model.matrix(~ active, data = x)[, -1], cencode = 0)})
sr_est_list <- lapply(sr_list, function(x){x$coef})
sr_se_list <- lapply(sr_list, function(x){sqrt(x$var)})
sr_ci_lower_list <- lapply(sr_list, function(x){x$coef - qnorm(0.975)*sqrt(x$var)})
sr_ci_upper_list <- lapply(sr_list, function(x){x$coef + qnorm(0.975)*sqrt(x$var)})
sr_study_res <- data.frame(cbind(unlist(sr_est_list), unlist(sr_se_list), unlist(sr_ci_lower_list), unlist(sr_ci_upper_list)))

## add sample sizes and event rates
sr_study_res$rec_events <- unlist(lapply(df_list, function(x){table(x$recovery)[2]}))
sr_study_res$rec_n <- unlist(lapply(sr_list, function(x){x$n}))

### pooled results (one-stage Fine-Gray stratified by study)
sr_pooled <- crrs(ftime = tico$t2recovery, fstatus = tico$recovery,
                  cov1 = model.matrix(~ active, data = tico)[, -1], strata = tico$trial, cencode = 0)
pooled_te <- sr_pooled$coef
pooled_se <- sqrt(sr_pooled$var)
pooled_lower <- pooled_te - qnorm(0.975)*pooled_se
pooled_upper <- pooled_te + qnorm(0.975)*pooled_se

# overall results
ov_sum <- c(pooled_te, pooled_se, pooled_lower, pooled_upper)
ov_events <- sum(sr_study_res$rec_events)
ov_n <- sum(sr_study_res$rec_n)

# put results in one dataframe
sr_all <- as.data.frame(rbind(sr_study_res, c(ov_sum, ov_events, ov_n)))
names(sr_all) <- c("Estimate", "SE", "Lower", "Upper", "Events", "N")
sr_all$Study <- c(names(df_list), "Overall")
sr_all$GS <- rep("Overall", nrow(sr_all))



###########################################
## Study results for interaction with GS
##########################################
# separate analyses by trial
sr_gs_list <- lapply(df_list, FUN = function(x){
  crr(ftime = x$t2recovery, fstatus = x$recovery,
      cov1 = model.matrix(~ active*gs_pos_0, data = x)[, -1], cencode = 0)
  })

# separate log RRR's for pos/neg
est_gs_list <- lapply(sr_gs_list, function(x){
  ests <- c(x$coef[1], x$coef[1] + x$coef[3])
  names(ests) <- c("neg", "pos")
  return(ests)})

#confidence intervals for positive and negative
ci_gs_lower_list <- lapply(sr_gs_list, function(x){
  se_pos <- sqrt(x$var[1,1] + x$var[3,3] + 2*x$var[1,3])
  se_neg <- sqrt(x$var[1,1])
  lower <- c(x$coef[1] - qnorm(0.975)*sqrt(x$var[1,1]), 
             x$coef[1] + x$coef[3] - qnorm(0.975)*se_pos)
  ses <- c(se_neg, se_pos)
  names(ses) <- c("neg", "pos")
  names(lower) <- c("neg", "pos")
  return(list(lower, ses))})

ci_gs_upper_list <- lapply(sr_gs_list, function(x){
  se_pos <- sqrt(x$var[1,1] + x$var[3,3] + 2*x$var[1,3])
  se_neg <- sqrt(x$var[1,1])
  upper <- c(x$coef[1] + qnorm(0.975)*sqrt(x$var[1,1]), 
             x$coef[1] + x$coef[3] + qnorm(0.975)*se_pos)
  ses <- c(se_neg, se_pos)
  names(upper) <- c("neg", "pos")
  names(ses) <- c("neg", "pos")
  return(list(upper, ses))})

n_gs_list <- lapply(df_list, function(x){
  res <- table(x$gs_pos_0)
  names(res) <- c("neg", "pos")
  return(res)
})

event_gs_list <- lapply(df_list, function(x){
  res <- table(x$gs_pos_0, x$recovery)[1:2,2]
  names(res) <- c("neg", "pos")
  return(res)
})

## put results into matrix
gs_study_res <- as.data.frame(cbind(unlist(est_gs_list),
                                    unlist(lapply(ci_gs_lower_list, function(l) l[[2]])),
                                    unlist(lapply(ci_gs_lower_list, function(l) l[[1]])),
                                    unlist(lapply(ci_gs_upper_list, function(l) l[[1]])),
                                    unlist(event_gs_list),
                                    unlist(n_gs_list)))
gs_study_res$Study <- rep(names(df_list), each = 2)
gs_study_res$GS <- rep(c("Neg", "Pos"), 5)
names(gs_study_res)[1:6] <- c("Estimate", "SE", "Lower", "Upper", "Events", "N")


###############################################
### Pooled Interaction Results - Sustained Recovery (one-stage F-G)
###############################################
sr_pooled_gs_within <- crrs(ftime = tico$t2recovery, fstatus = tico$recovery,
                     cov1 = model.matrix(~ active:trial + gs_pos_0:trial + active*gs_pos_0, data = tico)[, -1], strata = tico$trial, cencode = 0)

sr_pooled_gs_sharedTE <- crrs(ftime = tico$t2recovery, fstatus = tico$recovery,
                              cov1 = model.matrix(~ gs_pos_0:trial + active*gs_pos_0, data = tico)[, -1], strata = tico$trial, cencode = 0)

# get GS-specific estimates for treatment effect
sr_ests <- sr_pooled_gs_sharedTE$coef[c(1,7)]
sr_vcov <- sr_pooled_gs_sharedTE$var[c(1,7), c(1,7)]

# treatment effect for GS negative
gs_neg_te <- sr_ests[1]
gs_neg_se <- sqrt(sr_pooled_gs_sharedTE$var[1, 1])
gs_neg_lower <- sr_ests[1] - qnorm(0.975)*gs_neg_se
gs_neg_upper <- sr_ests[1] + qnorm(0.975)*gs_neg_se

# treatment effect for GS positive
gs_pos_te <- sum(sr_ests)
X <- c(1,1)
gs_pos_se <- sqrt(diag(t(X) %*% (sr_vcov %*% X)))
gs_pos_lower <- gs_pos_te - qnorm(0.975)*gs_pos_se
gs_pos_upper <- gs_pos_te + qnorm(0.975)*gs_pos_se

# overall results for gs pos/neg
sr_gs_both <- as.data.frame(rbind(c(gs_neg_te, gs_neg_se, gs_neg_lower, gs_neg_upper),
                                  c(gs_pos_te, gs_pos_se, gs_pos_lower, gs_pos_upper)))
sr_gs_both$Study <- rep("Overall",2)
sr_gs_both$GS <- c("Neg", "Pos")
sr_gs_both$Events <- (gs_study_res %>% group_by(GS)%>% summarize(Events = sum(Events)))$Events
sr_gs_both$N <- (gs_study_res %>% group_by(GS)%>% summarize(N = sum(N)))$N

names(sr_gs_both)[1:4] <- c("Estimate", "SE", "Lower", "Upper")

# add summary results to matrix of treatment effects
gs_res <- rbind(gs_study_res, sr_gs_both)

# add overall (no interaction results) to data frame
sr_full <- rbind(sr_all, gs_res)

########################
## Subgroup-specific funnel plots
########################
# seropositive
tico_meta_pos <- metagen(TE = subset(gs_res, GS == "Pos")$Estimate,
                          seTE = subset(gs_res, GS == "Pos")$SE)


#seronegative
tico_meta_neg <- metagen(TE = subset(gs_res, GS == "Neg")$Estimate,
                         seTE = subset(gs_res, GS == "Neg")$SE)


pdf("Manuscript Tables and Figures/Supplemental Figure 3 - TICO Funnel Plots.pdf",
    height = 4, width = 10)
par(mfrow = c(1,2))

# seronegative
funnel(tico_meta_neg, xlab = "log(RRR)")
mtext("(a) TICO/ACTIV-3: Seronegative", side = 3, line = 1, cex = 1)

# seropositive
funnel(tico_meta_pos, xlab = "log(RRR)")
mtext("(b) TICO/ACTIV-3: Seropositive", side = 3, line = 1, cex = 1)

dev.off()


###########################
## Full set of models for comparison
###########################

# ----------------------------- one stage ------------------------------------------------------------------------------------------#
# within-trial information only (stratified)
sr_pooled_gs_within <- crrs(ftime = tico$t2recovery, fstatus = tico$recovery,
                            cov1 = model.matrix(~ active:trial + gs_pos_0:trial + active*gs_pos_0, data = tico)[, -1], strata = tico$trial, cencode = 0)
ind1 <- which(colnames(model.matrix(~ active:trial + gs_pos_0:trial + active*gs_pos_0, data = tico)[, -1]) == "active:gs_pos_0")
sr_within_res <- cbind(sr_pooled_gs_within$coef, sqrt(diag(sr_pooled_gs_within$var)))[ind1,]

# shared treatment effect term (primary analysis)
sr_pooled_gs_sharedTE <- crrs(ftime = tico$t2recovery, fstatus = tico$recovery,
                              cov1 = model.matrix(~ gs_pos_0:trial + active*gs_pos_0, data = tico)[, -1], strata = tico$trial, cencode = 0)
ind2 <- which(colnames(model.matrix(~ gs_pos_0:trial + active*gs_pos_0, data = tico)[, -1]) == "gs_pos_0:active" )
sr_sharedTE_res <- cbind(sr_pooled_gs_sharedTE$coef, sqrt(diag(sr_pooled_gs_sharedTE$var)))[ind2,]

# between-study information (model 1)
sr_pooled_gs_btw <-crrs(ftime = tico$t2recovery, fstatus = tico$recovery,
                        cov1 = model.matrix(~ active*gs_pos_0, data = tico)[, -1], strata = tico$trial, cencode = 0)
ind3 <- which(colnames(model.matrix(~ active*gs_pos_0, data = tico)[, -1]) == "active:gs_pos_0" )
sr_btw_res <- cbind(sr_pooled_gs_btw$coef, sqrt(diag(sr_pooled_gs_btw$var)))[ind3,]

### within-study information only (centering covariate/model 3)
tico <- tico %>%
  group_by(trial) %>%
  mutate(mean_gs = mean(gs_pos_0),
         gs_cent = gs_pos_0 - mean_gs)

sr_pooled_cent <-crrs(ftime = tico$t2recovery, fstatus = tico$recovery,
                      cov1 = model.matrix(~ gs_pos_0:trial + active:mean_gs + active:gs_cent, data = tico)[, -1], strata = tico$trial, cencode = 0)
ind4 <-  which(colnames(model.matrix(~ gs_pos_0:trial + active:mean_gs + active:gs_cent, data = tico)[, -1]) == "active:gs_cent" )
sr_cent_res <- cbind(sr_pooled_cent$coef, sqrt(diag(sr_pooled_cent$var)))[ind4,]

one_stage_res <- rbind(sr_within_res, sr_cent_res, sr_sharedTE_res, sr_btw_res)

### ---------------------------- two-stage ---------------------------------------------------------------###
### first fit fine-gray model within each study
tico_fg_list <- lapply(df_list, function(dat){
  res <- crr(ftime = dat$t2recovery, fstatus = dat$recovery,
              cov1 = model.matrix(~ active*gs_pos_0, data = dat)[, -1], cencode = 0)
  return(res)
})

### between-study information (multivariate meta-analysis)
### coefs
tico_coefs_ls <- lapply(tico_fg_list, function(res){
  return(res$coef)
})

tico_coefs <- do.call("rbind", tico_coefs_ls)

### covariance matrices
tico_vcov <- lapply(tico_fg_list, function(res){
  return(res$var)
})

tico_mmeta <- mixmeta(tico_coefs, tico_vcov, method = "fixed")
mmeta_res <- cbind(tico_mmeta$coefficients, sqrt(diag(tico_mmeta$vcov)))["active:gs_pos_0",]

### pool only interaction terms
tico_int_vars <- lapply(tico_fg_list, function(res){
  return(res$var[3,3])
})

tico_umeta <- mixmeta(tico_coefs[,3], tico_int_vars, method = "fixed")
umeta_res <- cbind(tico_umeta$coefficients, sqrt(diag(tico_umeta$vcov)))

### combine all results into matrix (one-stage and two-stage)
all_tico_res <- rbind(sr_within_res, sr_cent_res, umeta_res, sr_sharedTE_res, sr_btw_res, mmeta_res)
colnames(all_tico_res) <- c("Estimate", "SE")
all_tico_res <- data.frame(all_tico_res)
all_tico_res$lowerCI <- all_tico_res$Estimate + qnorm(0.975)*all_tico_res$SE
all_tico_res$upperCI <- all_tico_res$Estimate - qnorm(0.975)*all_tico_res$SE

# flip reference group (seronegative/seropositive) and exponentiate 
all_tico_exp <- all_tico_res %>%
  select(Estimate, lowerCI, upperCI) %>%
  mutate_all(.funs = function(x){round(exp(-x), 2)})

all_tico_exp
rownames(all_tico_exp) <- c("Fully stratified (model 2)", "Centered (model 3)", "Univariate meta-analysis", "Shared treatment effect",
                            "No stratification (model 1)", "Multivariate meta-analysis")

all_tico_exp$`95% CI` <- paste0("[", all_tico_exp$lowerCI, ", ", all_tico_exp$upperCI, "]")

all_tico_exp <- all_tico_exp %>% select(Estimate, `95% CI`)
all_tico_exp
