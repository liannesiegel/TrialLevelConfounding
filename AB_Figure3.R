## Lianne Siegel
## Aggregation Bias Manuscript
## Figure 3 - funnel plots from hypothetical example
library(meta)
source("R Code/SynthEx_Figure2.R")

#################
### Scenario 1 (no aggregation bias)
#################
# odds ratios and variances from each study
scenario1_sum <- synth_ex2

# seropositive
scenario1_meta_pos <- metagen(TE = subset(scenario1_sum, x == "Positive")$coefs,
                          seTE = sqrt(subset(scenario1_sum, x == "Positive")$vars)*10)
funnel(scenario1_meta_pos, xlab = "logOR", ylab = "Standard Error x 10")

#seronegative
scenario1_meta_neg <- metagen(TE = subset(scenario1_sum, x == "Negative")$coefs,
                              seTE = sqrt(subset(scenario1_sum, x == "Negative")$vars)*10)
scen1_neg_fun <- funnel(scenario1_meta_neg, xlab = "logOR", ylab = "Standard Error x 10")


#################
### Scenario 2
### (confounding by study - covariate dist related to trt effect)
#################
scenario2_sum <- synth_ex
scenario2_meta_pos <- metagen(TE = subset(scenario2_sum, x == "Positive")$coefs,
                          seTE = sqrt(subset(scenario2_sum, x == "Positive")$vars)*10)
scen2_pos_fun <- funnel(scenario2_meta_pos, xlab = "logOR", ylab = "Standard Error x 10")

scenario2_meta_neg <- metagen(TE = subset(scenario2_sum, x == "Negative")$coefs,
                              seTE = sqrt(subset(scenario2_sum, x == "Negative")$vars)*10)
scen2_neg_fun <- funnel(scenario2_meta_neg, xlab = "logOR", ylab = "Standard Error x 10")


############
## Figure with 4 panels
#############

pdf("Figure 2 - Funnel Plots_Final.pdf", height = 5, width = 8)
par(mfrow = c(2,2))
par(mar = c(3, 3, 2, 1), oma = c(2, 2, 1, 1), xaxs = "i", yaxs = "i")

funnel(scenario1_meta_pos, xlab = "logOR", ylab = "Standard Error x 10",
       xlim = c(-1, 0.4), ylim = c(0.4, 0), axes = FALSE)
axis(side = 1, at = seq(-1,0.4,by=0.2), pos = 0.4)
axis(side = 2, at = seq(0, 0.4, by = 0.1), las = 2, pos = -1)
mtext("(a) Scenario 1: Seropositive", side = 3, line = 1, cex = 1)


funnel(scenario1_meta_neg, xlab = "logOR", ylab = "Standard Error x 10",
            xlim = c(-1, 0.4), ylim = c(0.3, 0), axes = FALSE)
axis(side = 1, at = seq(-1,0.4,by=0.2), pos = 0.3)
axis(side = 2, at = seq(0, 0.4, by = 0.1), las = 2, pos = -1)
mtext("(b) Scenario 1: Seronegative", side = 3, line = 1, cex = 1)

funnel(scenario2_meta_pos, xlab = "logOR", ylab = "Standard Error x 10",
            xlim = c(-1.3, 1), ylim = c(0.8, 0), axes = FALSE)
axis(side = 1, at = seq(-1.5, 1.2, by = 0.2), pos = 0.8)
axis(side = 2, at = seq(0, 0.8, by = 0.1), las = 2, pos = -1.3)
mtext("(c) Scenario 2: Seropositive", side = 3, line = 1, cex = 1)

mtext("Estimated Treatment Effect (logOR)", side = 1, outer = TRUE, cex = 1)
mtext("Standard Error x 10", side = 2, outer = TRUE, cex = 1)

funnel(scenario2_meta_neg, xlab = "logOR", ylab = "Standard Error x 10",
       xlim = c(-0.8, 0.1), ylim = c(0.35, 0), axes = FALSE)
axis(side = 1, at = seq(-0.8, 0.1, by = 0.2), pos = 0.35)
axis(side = 2, at = seq(0, 0.4, by = 0.1), las = 2, pos = -0.8)
mtext("(d) Scenario 2: Seronegative", side = 3, line = 1, cex = 1)

dev.off()


##############################
## Supplemental Scenario 4 (no trial-level confounding and true effect modification)
##############################

# odds ratios and variances from each study
scenario4_sum <- synth_ex_scen4

# seropositive
scenario4_meta_pos <- metagen(TE = subset(scenario4_sum, x == "Positive")$coefs,
                              seTE = sqrt(subset(scenario4_sum, x == "Positive")$vars)*10)

#seronegative
scenario4_meta_neg <- metagen(TE = subset(scenario4_sum, x == "Negative")$coefs,
                              seTE = sqrt(subset(scenario4_sum, x == "Negative")$vars)*10)

pdf("Manuscript Tables and Figures/Supplemental Figure 2 - Scenario 4 Funnel Plots.pdf",
    height = 4, width = 10)
par(mfrow = c(1,2))

# seronegative
funnel(scenario4_meta_neg, xlab = "logOR", ylab = "Standard Error x 10")
mtext("(a) Scenario 4: Seronegative", side = 3, line = 1, cex = 1)

# seropositive
funnel(scenario4_meta_pos, xlab = "logOR", ylab = "Standard Error x 10")
mtext("(b) Scenario 4: Seropositive", side = 3, line = 1, cex = 1)

dev.off()


