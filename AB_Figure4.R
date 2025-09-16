### Lianne Siegel
### Figure 4 - Aggregation Bias Manuscript
### Implementation of Richard Riley's and Peter Godolphin's approaches to 
### estimating subgroup-specific treatment effects 
### Two panels (one for Richard's approach and one for Peter's)
### Overlayed on Scenario 2 panel from Figure 1

library(mixmeta)
library(ggpubr)
library(marginaleffects)

source("R Code/SynthEx_Figure2.R")
source("R Code/AB_Table2.R")

#######################################################
## Subgroup-specific treatment effects - Riley et al.
#######################################################
# suggests multivariate meta-analysis ignoring correlation interaction terms
# same here as separately pooling main effect and interaction

# ------------ seronegative as reference group -------------------------------------#
# pooled interaction terms from Table 1 analysis ("AB_Table1.R")
scen2_umeta

# pool just main effect of treatment (treatment effect in seronegatives)
glm_t_coefs <- unlist(lapply(scen2_glm_list, function(x){
  t_coef <- summary(x)$coef["t","Estimate"]
  return(t_coef)
}))

glm_t_vars <- unlist(lapply(scen2_glm_list, function(x){
  t_vars <- vcov(x)["t","t"]
}))

# fit model 
t_meta <- mixmeta(glm_t_coefs, glm_t_vars, method = "fixed")
summary(t_meta)

# estimated TE in seronegatives
t_neg_rich <- summary(t_meta)$coef[1,"Estimate"]

# estimated TE in seropositives
t_pos_rich <- summary(t_meta)$coef[1,"Estimate"] + summary(scen2_umeta)$coef[1,"Estimate"]

# --------- seropositive as reference group -------------------------------------------------#
# refit glms in each study switching reference group
glm_study2 <- lapply(dat, function(x){
  x$x_rev <- fct_rev(factor(x$x))
  res <- glm(y ~ x_rev*t, data = x, family = "binomial")
  return(res)
})

## fit univariate meta-analysis model w/just interaction terms
glm_int_coefs2 <- lapply(glm_study2, function(x){
  int_coef <- summary(x)$coef[4,1]
})

glm_int_coefs2 <- unlist(glm_int_coefs2)

glm_int_vars2 <- lapply(glm_study2, function(x){
  int_vars <- vcov(x)[4,4]
})

int_meta2 <- mixmeta(glm_int_coefs2, glm_int_vars2, method = "fixed")
summary(int_meta2)

## pool main effects of treatment
glm_t_coefs2 <- unlist(lapply(glm_study2, function(x){
  int_coef <- summary(x)$coef[3,1]
}))

glm_t_vars2 <- unlist(lapply(glm_study2, function(x){
  t_vars <- vcov(x)[3,3]
}))

# fit model 
t_meta2 <- mixmeta(glm_t_coefs2, glm_t_vars2, method = "fixed")

# estimated treatment effect in seropositives
t_pos_rich2 <- summary(t_meta2)$coef[1, "Estimate"]
t_neg_rich2 <- t_pos_rich2 + summary(int_meta2)$coef[1, "Estimate"]

#######################################################
## Subgroup-specific treatment effects - Godolphin et al.
#######################################################
# need subgroup effects and their variances
# from SynthEx.R
neg_coefs_vec <- unlist(neg_coefs)
neg_vars_vec <- unlist(neg_vars)
pos_coefs_vec <- unlist(pos_coefs)
pos_vars_vec <- unlist(pos_vars)

muhat <- (1/(sum((1/neg_vars_vec) + (1/pos_vars_vec))))*sum((1/neg_vars_vec)*neg_coefs_vec + (1/pos_vars_vec)*pos_coefs_vec)
A1 <- -1*(1/(sum((1/neg_vars_vec) + (1/pos_vars_vec))))*sum(1/pos_vars_vec)
gammahat <- summary(scen2_umeta)$coef[1, "Estimate"]
beta1 <- muhat + A1*gammahat
beta2 <- beta1 + gammahat

##############################################################################
## Subgroup-specific treatment effects - Marginalizing over stratified model
##############################################################################

# res_onestage_ws2 has model results for fully stratified model (Scenario 2)
# estimated interaction (log odds ratio): -0.015784  
res_onestage_ws2 

#seronegative
newdat_neg <- datagrid(model = res_onestage_ws2,
                       t = c(0,1),
                       x = 0,
                       grid_type = "counterfactual")

pred_neg <- avg_comparisons(res_onestage_ws2,
                            newdata = newdat_neg,
                            variable = "t",
                            type = "response",
                            comparison = "lnor")

# pred_neg_balanced <- avg_comparisons(res_onestage_ws2,
#                                      newdata = "balanced",
#                                      variable = "t",
#                                      type = "response",
#                                      comparison = "lnor")

#seropositive
newdat_pos <- datagrid(model = res_onestage_ws2,
                       t = c(0,1),
                       x = 1,
                       grid_type = "counterfactual")

pred_pos <- avg_comparisons(res_onestage_ws2,
                            newdata = newdat_pos,
                            variable = "t",
                            type = "response",
                            comparison = "lnor")

# pred_pos_balanced <- avg_comparisons(res_onestage_ws2,
#                                      newdata = newdat_pos,
#                                      variable = "t",
#                                      type = "response",
#                                      comparison = "lnor")

#pred_pos$estimate - pred_neg$estimate


#######################################################
## Add to plot similar to Figure 1
#######################################################
# data from Scenario 2 for Figure 1
all_results

# -------------- Marginalizing ---------------------------------#
all_results_margin <- all_results
all_results_margin[11:12,1] <- c(pred_neg$estimate, pred_pos$estimate)
all_results_margin[11:12, "type"] <- c("Pooled")

ab_plot_margin <- ggplot(all_results_margin, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type,
                 fill = type)) + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 20)) + 
  scale_alpha_manual(values = c(0.7, 1)) + 
  scale_linetype_manual(values = c(1,2)) + 
  scale_fill_manual(values = c("pink", "black")) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "vertical") + 
  labs(x = "Baseline Serostatus",
       y = "Estimated logOR") +
  ggtitle("(c) ATE Marginalizing Across Studies") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape = guide_legend(nrow=2,byrow=TRUE),
         alpha = guide_legend(nrow=2,byrow=TRUE),
         color = guide_legend(nrow=2,byrow=TRUE)) +
  scale_size(range = c(1,10))

ab_plot_margin


# -------------- Riley et al. ---------------------------------#
# swap out pooled subgroup effects for new ones
# riley et al. - first reference group
all_results_riley <- all_results
all_results_riley[11:12,1] <- c(t_neg_rich, t_pos_rich)
all_results_riley2 <- rbind(all_results_riley, all_results_riley[11:12,])
all_results_riley2[13:14,1] <- c(t_neg_rich2, t_pos_rich2)
all_results_riley2[13:14,"study"] <- rep(8,2)
all_results_riley2[11:12, "type"] <- c("Pooled - Seronegative as Reference")
all_results_riley2[13:14, "type"] <- c("Pooled - Seropositive as Reference")

all_results_riley3 <- cbind(all_results_riley2,
                           c(rep("Study", 10),
                             rep("Seronegative", 2),
                             rep("Seropositive", 2)))

ab_plot_riley <- ggplot(all_results_riley2, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type,
                 fill = type)) + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 23, 20)) + 
  scale_alpha_manual(values = c(0.7, 0.7, 1)) + 
  scale_linetype_manual(values = c(1,1,2)) + 
  scale_fill_manual(values = c("green", "darkblue", "black")) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "vertical") + 
  labs(x = "Baseline Serostatus",
       y = "Estimated logOR") +
  ggtitle("(b) Pooling Main Effects and Interaction Separately") + 
  guides(fill=guide_legend(nrow=3,byrow=TRUE),
         shape = guide_legend(nrow=3,byrow=TRUE),
         alpha = guide_legend(nrow=3,byrow=TRUE),
         color = guide_legend(nrow=3,byrow=TRUE)) +
  scale_size(range = c(1,10))

ab_plot_riley

# -------- Godolphin et al. --------------------------------#
all_results_godolphin <- all_results
all_results_godolphin[11:12,1] <- c(beta1, beta2)
all_results_godolphin[11:12, "type"] <- c("Pooled")

ab_plot_godolphin <- ggplot(all_results_godolphin, aes(x = x, y = coefs)) + 
  geom_point(color = "black",
             aes(size = weights,
                 shape = type,
                 alpha = type,
                 fill = type)) + 
  geom_line(aes(group = study,
                linetype = type,
                alpha = type)) + 
  guides(size = FALSE) + 
  scale_shape_manual(values = c(23, 20)) + 
  scale_alpha_manual(values = c(0.7, 1)) + 
  scale_linetype_manual(values = c(1,2)) + 
  scale_fill_manual(values = c("purple", "black")) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "vertical") + 
  labs(x = "Baseline Serostatus",
       y = "Estimated logOR") +
  ggtitle("(d) Godolphin et al. (2023)") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape = guide_legend(nrow=2,byrow=TRUE),
         alpha = guide_legend(nrow=2,byrow=TRUE),
         color = guide_legend(nrow=2,byrow=TRUE)) +
  scale_size(range = c(1,10))

ab_plot_godolphin


# ------------ three panels --------------------------------------#
ab_plot1_title <- ab_plot1 +
  ggtitle("(a) Between-Study Information (Scenario 2)") + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "vertical") + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE),
         shape = guide_legend(nrow=2,byrow=TRUE),
         alpha = guide_legend(nrow=2,byrow=TRUE),
         color = guide_legend(nrow=2,byrow=TRUE)) +
  scale_linetype_manual(values = c(1,2)) 
  

ab_plot1_use <- ab_plot1_title +
  theme(plot.margin = margin(0.5, 0.5, 1.3, 0.5, "cm"))
ab_plot_margin_use <- ab_plot_margin + 
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
ab_plot_riley_use <- ab_plot_riley +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
ab_plot_godolphin_use <- ab_plot_godolphin +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))


ggarrange(ab_plot1_use, ab_plot_riley_use, ab_plot_margin_use, ab_plot_godolphin_use,
          nrow = 2,
          ncol = 2,
          common.legend = FALSE,
          #font.label = list(size = 10),
          legend = "bottom")


#ggsave("Manuscript Tables and Figures/Figure4 - SubgroupSpecificTEs_v4.eps", width = 10, height = 8, dpi = 600, device = cairo_pdf)


