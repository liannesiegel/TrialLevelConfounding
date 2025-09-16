### Lianne Siegel 
### Table 2 - Aggregation Bias Paper
### Results for Scenarios 1-3 with and without stratification
### Includes one-stage and two-stage approaches

library(mixmeta)
library(lme4)
source("R Code/SynthEx_Figure2.R")

#########################################
### Scenario 1 - no confounding by study
#########################################
### ------- one-stage ----------------------------------------###
### between-study information
dat_pooledscen1 <- dat_pooled2
res_onestage_bs1 <- glm(y ~ factor(x)*t + factor(s), data = dat_pooledscen1, family = "binomial")
summary(res_onestage_bs1)

### within-study information only (stratification)
res_onestage_ws1 <- glm(y ~ factor(x)*factor(s) + factor(s)*t + factor(x):t + factor(s), data = dat_pooledscen1, family = "binomial")
summary(res_onestage_ws1)

### within-study information only (center covariate)
# mean of covariate within trial
dat_pooledscen1 <- dat_pooledscen1 %>%
  group_by(s) %>%
  mutate(meanx = mean(x), x_cent = x - meanx)

# don't include term for mean x because same in all studies but study-specific terms for x
res_onestage_ws1_cent <- glm(y ~ factor(s)*x + t + x_cent:t, data = dat_pooledscen1, family = "binomial")
summary(res_onestage_ws1_cent)


# one-stage stratifying only main effect of x 
res_onestage_xonly <- glm(y ~ factor(x)*factor(s) + t + factor(x):t + factor(s), data = dat_pooledscen1, family = "binomial")
summary(res_onestage_xonly)

# one-stage stratifying only main effect of trt
res_onestage_tonly <- glm(y ~ factor(x) + factor(s)*t + factor(x):t + factor(s), data = dat_pooledscen1, family = "binomial")
summary(res_onestage_tonly)


### ------ two-stage -----------------------------------------###
### first fit logistic regression within each study
scen1dat_list <- dat2
scen1_glm_list <- lapply(scen1dat_list, function(dat){
  res <- glm(y ~ factor(x)*t, data = dat, family = "binomial")
  return(res)
})

### between-study information (multivariate meta-analysis)
### coefs
scen1_coefs_ls <- lapply(scen1_glm_list, function(res){
  return(coef(res)[2:4])
})

scen1_coefs <- do.call("rbind",scen1_coefs_ls)

### covariance matrices
scen1_vcov <- lapply(scen1_glm_list, function(res){
  return(vcov(res)[2:4, 2:4])
})

scen1_mmeta <- mixmeta(scen1_coefs, scen1_vcov, method = "fixed")

### pool only interaction terms
scen1_int_vars <- lapply(scen1_glm_list, function(res){
  return(vcov(res)[4,4])
})

scen1_umeta <- mixmeta(scen1_coefs[,3], scen1_int_vars, method = "fixed")

### ----------- collate results ------------------------------------- ###
### want interaction coef (diff in log odds ratio), [95% CI]; p-value 
### for each method (b-s one stage, two-stage; ws one stage, two-stage)
bs_1stage_scen1_ci <- confint.default(res_onestage_bs1)
bs_1stage_scen1 <- paste0(sprintf(coef(res_onestage_bs1)["factor(x)1:t"], fmt = '%#.2f'), 
                          " [",
                          sprintf(bs_1stage_scen1_ci["factor(x)1:t", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(bs_1stage_scen1_ci["factor(x)1:t", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_bs1)$coef["factor(x)1:t","Pr(>|z|)"], fmt = '%#.2f'))

ws_1stage_scen1_ci <- confint.default(res_onestage_ws1)
ws_1stage_scen1 <- paste0(sprintf(coef(res_onestage_ws1)["factor(x)1:t"], fmt = '%#.2f'), 
                          " [",
                          sprintf(ws_1stage_scen1_ci["factor(x)1:t", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(ws_1stage_scen1_ci["factor(x)1:t", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_ws1)$coef["factor(x)1:t","Pr(>|z|)"], fmt = '%#.2f'))


cent_1stage_scen1_ci <- confint.default(res_onestage_ws1_cent)
cent_1stage_scen1 <- paste0(sprintf(coef(res_onestage_ws1_cent)["t:x_cent"], fmt = '%#.2f'), 
                          " [",
                          sprintf(cent_1stage_scen1_ci ["t:x_cent", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(cent_1stage_scen1_ci ["t:x_cent", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_ws1_cent)$coef["t:x_cent","Pr(>|z|)"], fmt = '%#.2f'))


bs_2stage_scen1 <- paste0(sprintf(summary(scen1_mmeta)$coef["factor(x)1:t", "Estimate"], fmt = '%#.2f'), 
                          " [",
                          sprintf(summary(scen1_mmeta)$coef["factor(x)1:t", "95%ci.lb"], fmt = '%#.2f'),
                          ", ",
                          sprintf(summary(scen1_mmeta)$coef["factor(x)1:t", "95%ci.ub"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(scen1_mmeta)$coef["factor(x)1:t", "Pr(>|z|)"], fmt = '%#.2f'))


ws_2stage_scen1 <- paste0(sprintf(summary(scen1_umeta)$coef[1, "Estimate"], fmt = '%#.2f'), 
                          " [",
                          sprintf(summary(scen1_umeta)$coef[1, "95%ci.lb"], fmt = '%#.2f'),
                          ", ",
                          sprintf(summary(scen1_umeta)$coef[1, "95%ci.ub"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(scen1_umeta)$coef[1, "Pr(>|z|)"], fmt = '%#.2f'))

scenario1_all <- c(bs_1stage_scen1, bs_2stage_scen1, ws_1stage_scen1, cent_1stage_scen1, ws_2stage_scen1)


###############################################
### Scenario 2 - confounding by covariate dist
###############################################
### ------- one-stage ----------------------------------------###
### between-study information
dat_pooledscen2 <- dat_pooled
res_onestage_bs2 <- glm(y ~ factor(x)*t + factor(s), data = dat_pooledscen2, family = "binomial")
summary(res_onestage_bs2)

### within-study information only (stratification)
res_onestage_ws2 <- glm(y ~ factor(x)*factor(s) + factor(s)*t + factor(x):t + factor(s), data = dat_pooledscen2, family = "binomial")
summary(res_onestage_ws2)

### within-study information only (centering covariate)
dat_pooledscen2 <- dat_pooledscen2 %>%
  group_by(s) %>%
  mutate(meanx = mean(x), x_cent = x - meanx)

res_onestage_ws2_cent <- glm(y ~ factor(s)*factor(x) + t + t:meanx + x_cent:t , data = dat_pooledscen2, family = "binomial")
summary(res_onestage_ws2_cent)

### ------ two-stage -----------------------------------------###
### first fit logistic regression within each study
scen2dat_list <- dat
scen2_glm_list <- lapply(scen2dat_list, function(dat){
  res <- glm(y ~ factor(x)*t, data = dat, family = "binomial")
  return(res)
})

### between-study information (multivariate meta-analysis)
### coefs
scen2_coefs_ls <- lapply(scen2_glm_list, function(res){
  return(coef(res)[2:4])
})

scen2_coefs <- do.call("rbind",scen2_coefs_ls)

### covariance matrices
scen2_vcov <- lapply(scen2_glm_list, function(res){
  return(vcov(res)[2:4, 2:4])
})

scen2_mmeta <- mixmeta(scen2_coefs, scen2_vcov, method = "fixed")

### pool only interaction terms
scen2_int_vars <- lapply(scen2_glm_list, function(res){
  return(vcov(res)[4,4])
})

scen2_umeta <- mixmeta(scen2_coefs[,3], scen2_int_vars, method = "fixed")

### ----------- collate results ------------------------------------- ###
### want interaction coef (diff in log odds ratio), [95% CI]; p-value 
### for each method (b-s one stage, two-stage; ws one stage, two-stage)
bs_1stage_scen2_ci <- confint.default(res_onestage_bs2)
bs_1stage_scen2 <- paste0(sprintf(coef(res_onestage_bs2)["factor(x)1:t"], fmt = '%#.2f'), 
                          " [",
                          sprintf(bs_1stage_scen2_ci["factor(x)1:t", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(bs_1stage_scen2_ci["factor(x)1:t", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_bs2)$coef["factor(x)1:t","Pr(>|z|)"], fmt = '%#.2f'))

ws_1stage_scen2_ci <- confint.default(res_onestage_ws2)
ws_1stage_scen2 <- paste0(sprintf(coef(res_onestage_ws2)["factor(x)1:t"], fmt = '%#.2f'), 
                          " [",
                          sprintf(ws_1stage_scen2_ci["factor(x)1:t", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(ws_1stage_scen2_ci["factor(x)1:t", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_ws2)$coef["factor(x)1:t","Pr(>|z|)"], fmt = '%#.2f'))

cent_1stage_scen2_ci <- confint.default(res_onestage_ws2_cent)
cent_1stage_scen2 <- paste0(sprintf(coef(res_onestage_ws2_cent)["t:x_cent"], fmt = '%#.2f'), 
                            " [",
                            sprintf(cent_1stage_scen2_ci ["t:x_cent", "2.5 %"], fmt = '%#.2f'),
                            ", ",
                            sprintf(cent_1stage_scen2_ci ["t:x_cent", "97.5 %"], fmt = '%#.2f'),
                            "]; p = ",
                            sprintf(summary(res_onestage_ws2_cent)$coef["t:x_cent","Pr(>|z|)"], fmt = '%#.2f'))


bs_2stage_scen2 <- paste0(sprintf(summary(scen2_mmeta)$coef["factor(x)1:t", "Estimate"], fmt = '%#.2f'), 
                          " [",
                          sprintf(summary(scen2_mmeta)$coef["factor(x)1:t", "95%ci.lb"], fmt = '%#.2f'),
                          ", ",
                          sprintf(summary(scen2_mmeta)$coef["factor(x)1:t", "95%ci.ub"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(scen2_mmeta)$coef["factor(x)1:t", "Pr(>|z|)"], fmt = '%#.2f'))


ws_2stage_scen2 <- paste0(sprintf(summary(scen2_umeta)$coef[1, "Estimate"], fmt = '%#.2f'), 
                          " [",
                          sprintf(summary(scen2_umeta)$coef[1, "95%ci.lb"], fmt = '%#.2f'),
                          ", ",
                          sprintf(summary(scen2_umeta)$coef[1, "95%ci.ub"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(scen2_umeta)$coef[1, "Pr(>|z|)"], fmt = '%#.2f'))

scenario2_all <- c(bs_1stage_scen2, bs_2stage_scen2, ws_1stage_scen2, cent_1stage_scen2, ws_2stage_scen2)

###############################################
### Scenario 3 - confounding by allocation ratio
###############################################
### ------- one-stage ----------------------------------------###
### between-study information
dat_pooledscen3 <- dat_pooled3
res_onestage_bs3 <- glm(y ~ factor(x)*t + factor(s), data = dat_pooledscen3, family = "binomial")
summary(res_onestage_bs3)

### within-study information only (stratification)
res_onestage_ws3 <- glm(y ~ factor(x)*factor(s) + factor(s)*t + factor(x):t + factor(s), data = dat_pooledscen3, family = "binomial")
summary(res_onestage_ws3)

### within-study information only (centering covariate)
dat_pooledscen3 <- dat_pooledscen3 %>%
  group_by(s) %>%
  mutate(meanx = mean(x), x_cent = x - meanx)

### mean(x) same in all studies
res_onestage_ws3_cent <- glm(y ~ factor(s) + x:factor(s) + t + t:meanx + x_cent:t , data = dat_pooledscen3, family = "binomial")
summary(res_onestage_ws3_cent)

res_onestage_ws3_cent <- glm(y ~ factor(s) + x:factor(s) + t + x_cent:t , data = dat_pooledscen3, family = "binomial")
summary(res_onestage_ws3_cent)



### ------ two-stage -----------------------------------------###
### first fit logistic regression within each study
scen3dat_list <- dat3
scen3_glm_list <- lapply(scen3dat_list, function(dat){
  res <- glm(y ~ factor(x)*t, data = dat, family = "binomial")
  return(res)
})

### between-study information (multivariate meta-analysis)
### coefs
scen3_coefs_ls <- lapply(scen3_glm_list, function(res){
  return(coef(res)[2:4])
})

scen3_coefs <- do.call("rbind",scen3_coefs_ls)

### covariance matrices
scen3_vcov <- lapply(scen3_glm_list, function(res){
  return(vcov(res)[2:4, 2:4])
})

scen3_mmeta <- mixmeta(scen3_coefs, scen3_vcov, method = "fixed")

### pool only interaction terms
scen3_int_vars <- lapply(scen3_glm_list, function(res){
  return(vcov(res)[4,4])
})

scen3_umeta <- mixmeta(scen3_coefs[,3], scen3_int_vars, method = "fixed")

### ----------- collate results ------------------------------------- ###
### want interaction coef (diff in log odds ratio), [95% CI]; p-value 
### for each method (b-s one stage, two-stage; ws one stage, two-stage)
bs_1stage_scen3_ci <- confint.default(res_onestage_bs3)
bs_1stage_scen3 <- paste0(sprintf(coef(res_onestage_bs3)["factor(x)1:t"], fmt = '%#.2f'), 
                          " [",
                          sprintf(bs_1stage_scen3_ci["factor(x)1:t", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(bs_1stage_scen3_ci["factor(x)1:t", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_bs3)$coef["factor(x)1:t","Pr(>|z|)"], fmt = '%#.2f'))

ws_1stage_scen3_ci <- confint.default(res_onestage_ws3)
ws_1stage_scen3 <- paste0(sprintf(coef(res_onestage_ws3)["factor(x)1:t"], fmt = '%#.2f'), 
                          " [",
                          sprintf(ws_1stage_scen3_ci["factor(x)1:t", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(ws_1stage_scen3_ci["factor(x)1:t", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_ws3)$coef["factor(x)1:t","Pr(>|z|)"], fmt = '%#.2f'))

cent_1stage_scen3_ci <- confint.default(res_onestage_ws3_cent)
cent_1stage_scen3 <- paste0(sprintf(coef(res_onestage_ws3_cent)["t:x_cent"], fmt = '%#.2f'), 
                            " [",
                            sprintf(cent_1stage_scen3_ci ["t:x_cent", "2.5 %"], fmt = '%#.2f'),
                            ", ",
                            sprintf(cent_1stage_scen3_ci ["t:x_cent", "97.5 %"], fmt = '%#.2f'),
                            "]; p = ",
                            sprintf(summary(res_onestage_ws3_cent)$coef["t:x_cent","Pr(>|z|)"], fmt = '%#.2f'))


bs_2stage_scen3 <- paste0(sprintf(summary(scen3_mmeta)$coef["factor(x)1:t", "Estimate"], fmt = '%#.2f'), 
                          " [",
                          sprintf(summary(scen3_mmeta)$coef["factor(x)1:t", "95%ci.lb"], fmt = '%#.2f'),
                          ", ",
                          sprintf(summary(scen3_mmeta)$coef["factor(x)1:t", "95%ci.ub"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(scen3_mmeta)$coef["factor(x)1:t", "Pr(>|z|)"], fmt = '%#.2f'))


ws_2stage_scen3 <- paste0(sprintf(summary(scen3_umeta)$coef[1, "Estimate"], fmt = '%#.2f'), 
                          " [",
                          sprintf(summary(scen3_umeta)$coef[1, "95%ci.lb"], fmt = '%#.2f'),
                          ", ",
                          sprintf(summary(scen3_umeta)$coef[1, "95%ci.ub"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(scen3_umeta)$coef[1, "Pr(>|z|)"], fmt = '%#.2f'))

scenario3_all <- c(bs_1stage_scen3, bs_2stage_scen3, ws_1stage_scen3, cent_1stage_scen3, ws_2stage_scen3)


###########################
# Scenario 4 (True effect modification, no confounding)
#############################
### ------- one-stage ----------------------------------------###
### between-study information
dat_pooledscen4 <- dat_pooled_scen4
res_onestage_bs4 <- glm(y ~ factor(x)*t + factor(s), data = dat_pooledscen4, family = "binomial")
summary(res_onestage_bs4)

### within-study information only (stratification)
res_onestage_ws4 <- glm(y ~ factor(x)*factor(s) + factor(s)*t + factor(x):t + factor(s), data = dat_pooledscen4, family = "binomial")
summary(res_onestage_ws4)

### within-study information only (centering covariate)
dat_pooledscen4 <- dat_pooledscen4 %>%
  group_by(s) %>%
  mutate(meanx = mean(x), x_cent = x - meanx)

res_onestage_ws4_cent <- glm(y ~ factor(s)*factor(x) + t + t:meanx + x_cent:t , data = dat_pooledscen4, family = "binomial")
summary(res_onestage_ws4_cent)

### ------ two-stage -----------------------------------------###
### first fit logistic regression within each study
scen4dat_list <- dat_scen4
scen4_glm_list <- lapply(scen4dat_list, function(dat){
  res <- glm(y ~ factor(x)*t, data = dat, family = "binomial")
  return(res)
})

### between-study information (multivariate meta-analysis)
### coefs
scen4_coefs_ls <- lapply(scen4_glm_list, function(res){
  return(coef(res)[2:4])
})

scen4_coefs <- do.call("rbind",scen4_coefs_ls)

### covariance matrices
scen4_vcov <- lapply(scen4_glm_list, function(res){
  return(vcov(res)[2:4, 2:4])
})

scen4_mmeta <- mixmeta(scen4_coefs, scen4_vcov, method = "fixed")

### pool only interaction terms
scen4_int_vars <- lapply(scen4_glm_list, function(res){
  return(vcov(res)[4,4])
})

scen4_umeta <- mixmeta(scen4_coefs[,3], scen4_int_vars, method = "fixed")

### ----------- collate results ------------------------------------- ###
### want interaction coef (diff in log odds ratio), [95% CI]; p-value 
### for each method (b-s one stage, two-stage; ws one stage, two-stage)
bs_1stage_scen4_ci <- confint.default(res_onestage_bs4)
bs_1stage_scen4 <- paste0(sprintf(coef(res_onestage_bs4)["factor(x)1:t"], fmt = '%#.2f'), 
                          " [",
                          sprintf(bs_1stage_scen4_ci["factor(x)1:t", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(bs_1stage_scen4_ci["factor(x)1:t", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_bs4)$coef["factor(x)1:t","Pr(>|z|)"], fmt = '%#.2f'))

ws_1stage_scen4_ci <- confint.default(res_onestage_ws4)
ws_1stage_scen4 <- paste0(sprintf(coef(res_onestage_ws4)["factor(x)1:t"], fmt = '%#.2f'), 
                          " [",
                          sprintf(ws_1stage_scen4_ci["factor(x)1:t", "2.5 %"], fmt = '%#.2f'),
                          ", ",
                          sprintf(ws_1stage_scen4_ci["factor(x)1:t", "97.5 %"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(res_onestage_ws4)$coef["factor(x)1:t","Pr(>|z|)"], fmt = '%#.2f'))

cent_1stage_scen4_ci <- confint.default(res_onestage_ws4_cent)
cent_1stage_scen4 <- paste0(sprintf(coef(res_onestage_ws4_cent)["t:x_cent"], fmt = '%#.2f'), 
                            " [",
                            sprintf(cent_1stage_scen4_ci ["t:x_cent", "2.5 %"], fmt = '%#.2f'),
                            ", ",
                            sprintf(cent_1stage_scen4_ci ["t:x_cent", "97.5 %"], fmt = '%#.2f'),
                            "]; p = ",
                            sprintf(summary(res_onestage_ws4_cent)$coef["t:x_cent","Pr(>|z|)"], fmt = '%#.2f'))


bs_2stage_scen4 <- paste0(sprintf(summary(scen4_mmeta)$coef["factor(x)1:t", "Estimate"], fmt = '%#.2f'), 
                          " [",
                          sprintf(summary(scen4_mmeta)$coef["factor(x)1:t", "95%ci.lb"], fmt = '%#.2f'),
                          ", ",
                          sprintf(summary(scen4_mmeta)$coef["factor(x)1:t", "95%ci.ub"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(scen4_mmeta)$coef["factor(x)1:t", "Pr(>|z|)"], fmt = '%#.2f'))


ws_2stage_scen4 <- paste0(sprintf(summary(scen4_umeta)$coef[1, "Estimate"], fmt = '%#.2f'), 
                          " [",
                          sprintf(summary(scen4_umeta)$coef[1, "95%ci.lb"], fmt = '%#.2f'),
                          ", ",
                          sprintf(summary(scen4_umeta)$coef[1, "95%ci.ub"], fmt = '%#.2f'),
                          "]; p = ",
                          sprintf(summary(scen4_umeta)$coef[1, "Pr(>|z|)"], fmt = '%#.2f'))

scenario4_all <- c(bs_1stage_scen4, bs_2stage_scen4, ws_1stage_scen4, cent_1stage_scen4, ws_2stage_scen4)



##################################
### Make results table
###################################
table_all <- cbind(scenario1_all, scenario2_all, scenario3_all)
colnames(table_all) <- c("Scenario 1", "Scenario 2", "Scenario 3")
table_all2 <- cbind(c("Between-Study Information",
                      "Between-Study Information",
                      "Within-Study Information Only",
                      "Within-Study Information Only",
                      "Within-Study Information Only"),
                    c("One-stage",
                      "Two-stage",
                      "One-stage (Stratification)",
                      "One-stage (Centering Covariate)",
                      "Two-stage"),
                    table_all)
## export
#write.table(table_all2, file = "table2_v2.txt", sep = "_", quote = FALSE, row.names = T)

#################
## Supplementary Table for Scenario 4
#################
table_scen4_labeled <- cbind(c("Between-Study Information",
                      "Between-Study Information",
                      "Within-Study Information Only",
                      "Within-Study Information Only",
                      "Within-Study Information Only"),
                    c("One-stage",
                      "Two-stage",
                      "One-stage (Stratification)",
                      "One-stage (Centering Covariate)",
                      "Two-stage"),
                    scenario4_all)
colnames(table_scen4_labeled) <- c("Model", "Stages", "Scenario 4")

