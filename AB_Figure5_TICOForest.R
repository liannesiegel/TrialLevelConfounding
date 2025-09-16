### Lianne Siegel
### Aggregation Bias Mansucript
### Figure 5 - Forest Plots of TICO/ACTIV-3 Results
source("./R Code/AB_TICOAnalysis.R")
library(forestploter)

### need to reorder levels of Study and GS variables
sr_full$Study2 <- factor(sr_full$Study, levels = c(names(df_list), "Overall"))
sr_full$GS2 <- factor(sr_full$GS, levels = c("Neg", "Pos", "Overall"),
                      labels = c("Negative", "Positive", "Overall"))

#######################################################
### Forest plot of Sustained Recovery Results
#######################################################

# ----------------- with forestploter ---------------------------------------#
# create empty row for plotted RRR's
sr_full$` ` <- paste(rep(" ", 18), collapse = " ")

sr_full$`RRR (95% CI)` <- ifelse(is.na(sr_full$Estimate), "",
                                 sprintf("%.2f (%.2f to %.2f)",
                                         exp(sr_full$Estimate), 
                                         exp(sr_full$Lower), exp(sr_full$Upper)))

# sort dataset by study
sr_full2 <- arrange(sr_full, Study2, GS2) %>%
  dplyr::select(Study2, GS2, ` `, `RRR (95% CI)`, Estimate, SE, Lower, Upper, Events, N)

# add blank rows 
sr_full3 <- sr_full2
for (i in c(1, 5, 9, 13, 17, 21)){
  sr_full3 <- sr_full3%>% add_row(.before = i)
}

sr_full3$Study3 <- c(names(df_list)[1], rep(" ", 3), names(df_list)[2], rep(" ", 3),
                     names(df_list)[3], rep(" ", 3), names(df_list)[4], rep(" ", 3),
                     names(df_list)[5], rep(" ", 3),
                     "Overall", rep(" ", 3))

sr_full3$`GenScript Antibody` <- ifelse(is.na(sr_full3$GS2), "", sr_full3$GS2)

sr_full3$`GenScript Antibody` <- factor(factor(sr_full3$`GenScript Antibody`), levels = c("", "1", "2", "3"),
                                        labels = c("", "Negative", "Positive", "Overall"))


sr_full3$`RRR (95% CI)` <- ifelse(is.na(sr_full3$`RRR (95% CI)`), "", 
                                  sr_full3$`RRR (95% CI)`)

sr_full3$Trial <- sr_full3$Study3

sr_full4 <- sr_full3 %>% select(Trial, `GenScript Antibody`, ` `, `RRR (95% CI)`, Estimate,
                                SE, Lower, Upper, Events, N)
sr_full4$` ` <- paste(rep(" ", 18), collapse = " ")

sr_full4 <- sr_full4 %>% mutate(N = replace_na(as.character(N), ""),
                                Events = replace_na(as.character(Events), ""))

# make forest plot
fig4a <- forest(sr_full4[,c(1:4,9:10)],
             est = exp(sr_full4$Estimate),
             lower = exp(sr_full4$Lower),
             upper = exp(sr_full4$Upper),
             size = sr_full4$SE,
             ci_column= 3,
             ref_line = 1,
             is_summary = c(rep(FALSE, nrow(sr_full4)-3), rep(TRUE, 3)),
             arrow_lab = c("Placebo Better", "Treatment Better"),
             xlim = c(0.5, 2))

plot(fig4a)


#######################################################
### Forest Plot of Interaction Coefficients
#######################################################

# list of interaction coefficients
int_gs_list <- lapply(sr_gs_list, function(x){
  ests <- c(x$coef[3])
  return(ests)})

# lower bound of coefficient
int_gs_list <- lapply(sr_gs_list, function(x){
  ests <- log(summary(x)$conf.int["active:gs_pos_0",c("2.5%", "97.5%")])
  ses <- summary(x)$coef["active:gs_pos_0", c("coef", "se(coef)")]
  ests <- c(ses, ests)
  return(ests)})

int_gs_res <- data.frame(t(data.frame(int_gs_list)))

### add overall results
### within-study information only (stratification of all main effects)
overall_int_with_int  <- sr_pooled_gs_within$coef[11]
overall_int_with_se  <- sqrt(sr_pooled_gs_within$var[11,11])
overall_int_with_lower <- overall_int_with_int - qnorm(0.975)*overall_int_with_se
overall_int_with_upper <- overall_int_with_int + qnorm(0.975)*overall_int_with_se
overall_int_with <- c(overall_int_with_int,
                      overall_int_with_se,
                      overall_int_with_lower,
                      overall_int_with_upper)

### no stratification of main effect of treatment
overall_int_sharedTE_int  <- sr_pooled_gs_sharedTE$coef[7]
overall_int_sharedTE_se  <- sqrt(sr_pooled_gs_sharedTE$var[7,7])
overall_int_sharedTE_lower <- overall_int_sharedTE_int - qnorm(0.975)*overall_int_sharedTE_se
overall_int_sharedTE_upper <- overall_int_sharedTE_int + qnorm(0.975)*overall_int_sharedTE_se
overall_int_sharedTE <- c(overall_int_sharedTE_int,
                          overall_int_sharedTE_se,
                          overall_int_sharedTE_lower,
                          overall_int_sharedTE_upper)

int_gs_res <- rbind(int_gs_res, overall_int_with, overall_int_sharedTE)

int_gs_res$Study <- c(names(sr_gs_list),
                      "Overall - Within-Study Information Only",
                      "Overall - Common Main Effect of Treatment")
names(int_gs_res) <- c("Estimate", "SE", "Lower", "Upper", "Trial")


# create empty row for plotted lines
int_gs_res$` ` <- paste(rep(" ", 22), collapse = " ")

# exponentiate results and flip
int_gs_res$`rRRR* (95% CI)` <- ifelse(is.na(int_gs_res$Estimate), "",
                                      sprintf("%.2f (%.2f to %.2f)",
                                              1/exp(int_gs_res$Estimate), 
                                              1/exp(int_gs_res$Upper),
                                              1/exp(int_gs_res$Lower)))


# make forest plot
fig4b <- forest(int_gs_res[,c(5:7)],
                est = 1/exp(int_gs_res$Estimate),
                lower =  1/exp(int_gs_res$Upper),
                upper = 1/exp(int_gs_res$Lower),
                size = int_gs_res$SE,
                ci_column= 2,
                ref_line = 1,
                is_summary = c(rep(FALSE, nrow(int_gs_res)-2), rep(TRUE, 2)),
                arrow_lab = c("Treatment Better in GS Positive", "Treatment Better in GS Negative"),
                xlim = c(0, 3),
                footnote = "*Ratio of RRR in seronegative participants \
             to RRR in seropositive participants")
plot(fig4b)

###################################
### Export Figure w/Two Panels
### This code no longer works w/updated version of forestploter
### Can export each figure above separately instead (e.g. using ggsave())
#####################################
# pdf("Figure 5 - TICO.pdf", height = 10, width = 10)
# par(mfrow = c(1,2), mar = c(3, 3, 2, 1), oma = c(2, 2, 1, 1))
# plot(fig4a)
# mtext("(a) ", side = 3, line = 1, cex = 1)
# plot(fig4b)
# mtext("(b) ", side = 3, line = 1, cex = 1)
# 
# dev.off()









