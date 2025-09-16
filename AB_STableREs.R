### Lianne Siegel
### Supplementary Table for Aggregation Bias paper
### Synthetic examples with random main effects 

library(mixmeta)
library(lme4)
source("R Code/SynthEx_Figure2.R")
source("R Code/AB_Table2.R")

#########################################
### Scenario 1 - no confounding by study
#########################################
### ------- one-stage random main effects  ----------------------------------------###
### between-study information
dat_pooledscen1 <- dat_pooled2
dat_pooledscen1$x_fac <- factor(dat_pooledscen1$x)

# only include random treatment effect for t as prognostic effect of x is 
# the same across all studies (including both results in singularity warning)
res_onestage_bs1_re <- glmer(y ~ x_fac*t + factor(s) + (0 + t | s), data = dat_pooledscen1, family = "binomial")
summary(res_onestage_bs1_re)
summary(res_onestage_bs1_re)$coef
summary(res_onestage_bs1_re)$coef[,"Estimate"] - qnorm(0.975)*summary(res_onestage_bs1_re)$coef[,"z value"]
summary(res_onestage_bs1_re)$coef[,"Estimate"] + qnorm(0.975)*summary(res_onestage_bs1_re)$coef[,"z value"]



###############################################
### Scenario 2 - confounding by covariate dist
###############################################
### ------- one-stage ----------------------------------------###
### between-study information
dat_pooledscen2$x_fac <- factor(dat_pooledscen2$x)
# only include random treatment effect for t as prognostic effect of x is 
# the same across all studies (including both results in singularity warning)
res_onestage_bs2_re <- glmer(y ~ x_fac*t + factor(s) + (0 + t | s), data = dat_pooledscen2, family = "binomial")
#res_onestage_bs2_re_xandt <- glmer(y ~ x_fac*t + factor(s) + (0 + t + x| s), data = dat_pooledscen2, family = "binomial")

summary(res_onestage_bs2_re)
summary(res_onestage_bs2_re)$coef
re_2_lowerci <- summary(res_onestage_bs2_re)$coef[,"Estimate"] - qnorm(0.975)*summary(res_onestage_bs2_re)$coef[,"Std. Error"]
re_2_upperci <- summary(res_onestage_bs2_re)$coef[,"Estimate"] + qnorm(0.975)*summary(res_onestage_bs2_re)$coef[,"Std. Error"]

paste0(sprintf(summary(res_onestage_bs2_re)$coef["x_fac1:t", "Estimate"], fmt = '%#.2f'), 
       " [",
       sprintf(re_2_lowerci["x_fac1:t"], fmt = '%#.2f'),
       ", ",
       sprintf(re_2_upperci["x_fac1:t"], fmt = '%#.2f'),
       "]; p = ",
       sprintf(summary(res_onestage_bs2_re)$coef["x_fac1:t","Pr(>|z|)"], fmt = '%#.2f'))

### random interaction effect
res_onestage_bs2_re_int <- glmer(y ~ x_fac*t + factor(s) + (0 + t + x_fac:t| s), data = dat_pooledscen2, family = "binomial")




###############################################
### Scenario 3 - confounding by allocation ratio
###############################################
### ------- one-stage ----------------------------------------###
### between-study information
dat_pooledscen3 <- dat_pooled3
dat_pooledscen3$x_fac <- factor(dat_pooledscen3$x)
# only include random treatment effect for x as treatment effect is
# the same across all studies (including both results in singularity warning)
res_onestage_bs3_re <- glmer(y ~ x_fac*t + factor(s) + (0 + x | s), data = dat_pooledscen3, family = "binomial")
#res_onestage_bs3_re_xandt <- glmer(y ~ x_fac*t + factor(s) + (0 + t + x| s), data = dat_pooledscen3, family = "binomial")

summary(res_onestage_bs3_re)
summary(res_onestage_bs3_re)$coef
re_3_lowerci <- summary(res_onestage_bs3_re)$coef[,"Estimate"] - qnorm(0.975)*summary(res_onestage_bs3_re)$coef[,"Std. Error"]
re_3_upperci <- summary(res_onestage_bs3_re)$coef[,"Estimate"] + qnorm(0.975)*summary(res_onestage_bs3_re)$coef[,"Std. Error"]

paste0(sprintf(summary(res_onestage_bs3_re)$coef["x_fac1:t", "Estimate"], fmt = '%#.2f'), 
       " [",
       sprintf(re_3_lowerci["x_fac1:t"], fmt = '%#.2f'),
       ", ",
       sprintf(re_3_upperci["x_fac1:t"], fmt = '%#.2f'),
       "]; p = ",
       sprintf(summary(res_onestage_bs3_re)$coef["x_fac1:t","Pr(>|z|)"], fmt = '%#.2f'))


##################################
## Explore impact of sample size shrinkage by sampling from datasets (Scenario 2)
###################################
set.seed(123)
# reducing sample size does not seem to change much
dat_pooledscen2_sub <- dat_pooledscen2[sample(1:nrow(dat_pooledscen2), 1000, replace = FALSE), ]
res_onestage_bs2_re_sub <- glmer(y ~ x_fac*t + factor(s) + (0 + t | s), data = dat_pooledscen2_sub, family = "binomial")

summary(res_onestage_bs2_re_sub)


###############################
# Scenario 4
##############################
### ------- one-stage ----------------------------------------###
### between-study information
dat_pooledscen4 <- dat_pooled4
dat_pooledscen4$x_fac <- factor(dat_pooledscen4$x)
# only include random treatment effect for t as prognostic effect of x is 
# the same across all studies (including both results in singularity warning)
res_onestage_bs4_re <- glmer(y ~ x_fac*t + factor(s) + (0 + t | s), data = dat_pooledscen4, family = "binomial")
#res_onestage_bs4_re_xandt <- glmer(y ~ x_fac*t + factor(s) + (0 + t + x| s), data = dat_pooledscen4, family = "binomial")

summary(res_onestage_bs4_re)
summary(res_onestage_bs4_re)$coef
re_4_lowerci <- summary(res_onestage_bs4_re)$coef[,"Estimate"] - qnorm(0.975)*summary(res_onestage_bs4_re)$coef[,"Std. Error"]
re_4_upperci <- summary(res_onestage_bs4_re)$coef[,"Estimate"] + qnorm(0.975)*summary(res_onestage_bs4_re)$coef[,"Std. Error"]

paste0(sprintf(summary(res_onestage_bs4_re)$coef["x_fac1:t", "Estimate"], fmt = '%#.2f'), 
       " [",
       sprintf(re_4_lowerci["x_fac1:t"], fmt = '%#.2f'),
       ", ",
       sprintf(re_4_upperci["x_fac1:t"], fmt = '%#.2f'),
       "]; p = ",
       sprintf(summary(res_onestage_bs4_re)$coef["x_fac1:t","Pr(>|z|)"], fmt = '%#.2f'))



