library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(gplots)
library(stringr)
library(splitstackshape)
library(forcats)
library(ggridges)
library(DescTools)
library(pheatmap)
library(broom)
library(rstatix)
library(ggsignif)


##*##
## This script produces all of the plots with a generalized input.

##Specify the input/output folder 
input.folder<-"~/Documents/StudentData/Elias/Elias_smFISH_paper2026/Github_101526/smFISH_analysis/wnt1.holdopen/"
setwd(input.folder)

#Specify an input spot count datafile with names of each sample "treatment" are column headers and each row is spot counts for each counted cell
input.datafile<-"307_output_all_cells_including_zero.csv"

#Specify a unique name for the dataset analyzed so that plots will save with unique names
dataset.name<-"wnt1.holdopen.307"

#Specify the order of conditions that you want to plot (conditions are defined to be groups of treatments)
conditions.order<-c("Incision_24h","Agarose_24h","Amputation_24h")


#Generalized Analysis
setwd(input.folder)
data.spots<-read.table(input.datafile, header=TRUE, sep=",")
str(data.spots)

#sample.order<-colnames(data.spots)
treatment.order<-c("Incision_24h_1","Incision_24h_2","Incision_24h_3","Agarose_24h_1", "Agarose_24h_2","Agarose_24h_3","Agarose_24h_4","Amputation_24h_1","Amputation_24h_2","Amputation_24h_3","Amputation_24h_4")
sample.order<-treatment.order

data.spots.long<-data.spots %>% pivot_longer(1:ncol(data.spots),names_to="treatment", values_to="mRNA.per.cell")
#wnt1.306b.spots.long<-factor(wnt1.306b.spots.long$treatment,levels=levels(wnt1.306b.spots.long$treatment))
str(data.spots.long)

data.spots.long$treatment <- factor(data.spots.long$treatment, levels=sample.order)

#sample.order<-levels(wnt1.306b.spots.long$treatment)
#sample.order<-c("Incision_24h_1","Incision_24h_2","Incision_24h_3","Agarose_24h_1","Agarose_24h_2","Agarose_24h_3","Agarose_24h_4","Amputation_24h_1","Amputation_24h_2","Amputation_24h_3","Amputation_24h_4")

#wnt1.306b.spots.long$treatment = factor(wnt1.306b.spots.long$treatment, levels=sample.order)

data.spots.long <-data.spots.long %>% mutate(condition=treatment) %>% separate_wider_delim(condition, delim="_", names= c("type","time","sample"))
data.spots.long <-data.spots.long %>% mutate(conditions=paste0(type,"_",time))


#wnt1.holdopen.spots.long.average$type = factor(wnt1.holdopen.spots.long.average$type, levels=c("Incision", "Agarose", "Amputation"))
#conditions.order<-c("Incision_0h","Incision_6h","Incision_12h","Amputation_0h","Amputation_6h","Amputation_12h")
#conditions.order<-c("Incision_0h","Amputation_0h","Incision_6h","Amputation_6h","Incision_12h","Amputation_12h")

data.spots.long %>% ggplot(aes(x=factor(conditions,levels=conditions.order),y=mRNA.per.cell,color=type)) +
  geom_jitter(alpha=0.3)+
  scale_color_manual(values = c("Incision" = "#4bb254", "Amputation" = "#1E90FF","Agarose"="orange"))+
  #  geom_violin()+
  theme_classic()+
  ylim(0,200)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 15))+
  ggtitle(paste0(dataset.name," mRNA/cell (sum replicates)"))+
  theme(plot.title = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")+
  theme(panel.border = element_blank())+
  theme(
    plot.background = element_rect(fill = "white", colour=NA), # Entire plot background
    panel.background = element_rect(fill = "white"), # Plot panel background
    panel.border = element_blank(),
    axis.text = element_text(color = "black"), # Axis tick labels
    axis.title = element_text(color = "black"), # Axis titles
    title= element_text(color = "black"), # plot title
    #    panel.grid.major = element_line(color = "grey30"), # Major grid lines
    #    panel.grid.minor = element_line(color = "grey20") # Minor grid lines
    axis.line = element_line(color = "black"), # Axis lines
    axis.ticks = element_line(colour = "black")
  )
ggsave(paste0(dataset.name,".mRNApercell.repsum.png"),width=3,height=5)


data.spots.long %>% ggplot(aes(x=factor(treatment,levels=treatment.order),y=mRNA.per.cell,color=type)) +
  geom_jitter(alpha=0.3)+
  scale_color_manual(values = c("Incision" = "#4bb254", "Amputation" = "#1E90FF","Agarose"="orange"))+
  #  geom_violin()+
  theme_classic()+
  ylim(0,200)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text= element_text(face = "bold"),text = element_text(size = 15))+
  ggtitle(paste0(dataset.name," mRNA/cell (each replicate)"))+
  theme(plot.title = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")+
  theme(panel.border = element_blank())+
  theme(
    plot.background = element_rect(fill = "white", colour=NA), # Entire plot background
    panel.background = element_rect(fill = "white"), # Plot panel background
    panel.border = element_blank(),
    axis.text = element_text(color = "black"), # Axis tick labels
    axis.title = element_text(color = "black"), # Axis titles
    title= element_text(color = "black"), # plot title
    #    panel.grid.major = element_line(color = "grey30"), # Major grid lines
    #    panel.grid.minor = element_line(color = "grey20") # Minor grid lines
    axis.line = element_line(color = "black"), # Axis lines
    axis.ticks = element_line(colour = "black")
  )
ggsave(paste0(dataset.name,".mRNApercell.eachrep.png"),width=6,height=8)


library(ggridges)
data.spots.long %>% filter(mRNA.per.cell > 0) %>% ggplot(aes(x=log(mRNA.per.cell+1,10), y=treatment,fill=type))+
  geom_density_ridges(scale=2,alpha=0.6,jittered_points = TRUE, position = position_raincloud(height = 0.2),
                      point_alpha = 0.3, point_size = 0.5)+
  scale_fill_manual(values = c("Incision" = "#4bb254", "Amputation" = "#1E90FF","Agarose"="orange"))+
#  xlim(0.1,max(log(data.spots.long$mRNA.per.cell+1,10)))+
  theme(legend.position = "none")+
  theme_classic()+theme(legend.position = "none")
ggsave(paste0(dataset.name,".ggridges.jittered.png"),width=3,height=5)


data.spots.long %>% filter(mRNA.per.cell > 0) %>% ggplot(aes(x=log(mRNA.per.cell+1,10), y=treatment,fill=type))+
  geom_density_ridges(scale=3,alpha=0.6)+
  scale_fill_manual(values = c("Incision" = "#4bb254", "Amputation" = "#1E90FF","Agarose"="orange"))+
#  xlim(0.1,max(log(wnt1.306b.spots.long$mRNA.per.cell+1,10)))+
  theme(legend.position = "none")+
  theme_classic()+theme(legend.position = "none")
ggsave(paste0(dataset.name,".ggridges.png"),width=3,height=5)











###*
## ============================================================
## Hurdle NB GLMM for smFISH data with heavy zero-inflation from
## rare responsive cells (wnt1/notum induction at wound sites)
##
## Hurdle negative binomial GLMM, which models:
##   (1) P(cell is responsive), via a logistic hurdle component
##   (2) intensity | responsive, via a zero-truncated NB
## both with the same condition:replicate random-effect structure
## used throughout this analysis. Zeros are NOT manually removed
## from the data -- glmmTMB routes each row to the correct part of
## the likelihood automatically based on whether count == 0.
## ============================================================

library(glmmTMB)
library(DHARMa)
library(emmeans)
library(dplyr)
library(ggplot2)

## ------------------------------------------------------------
## 0. Expected data
## ------------------------------------------------------------
## df: one row per cell, columns:
##   condition : factor, e.g. "Amputation_12h", "Incision_12h", ...
##   replicate : unique ID per biological replicate
##   count     : mRNA spot count (INCLUDING zeros -- do not filter)
##
## df <- read.csv("your_data.csv")   # reshape from wide (as uploaded)
##                                    # to long format first if needed
## df$condition <- factor(df$condition)
## df$replicate <- factor(df$replicate)


#format my data to Claude's generalized solution
df<-na.omit(data.spots.long)
colnames(df)<-c("replicate","count","type","time","sample","condition")
df$condition<-factor(df$condition)
#model <- glm.nb(count ~ condition, data = df)

## Set your reference/control level explicitly so contrasts make sense:
df$condition <- relevel(df$condition, ref = "Incision_24h")


## ------------------------------------------------------------
## 1. Sanity check: confirm overdispersion persists among nonzero
##    cells before committing to NB over Poisson for the intensity part
## ------------------------------------------------------------
nonzero_check <- df %>%
  filter(count > 0) %>%
  group_by(condition) %>%
  summarise(
    n_nonzero = n(),
    mean_nz = mean(count),
    var_nz = var(count),
    var_mean_ratio = var_nz / mean_nz,
    .groups = "drop"
  )
print(nonzero_check)
## var_mean_ratio >> 1 (especially well above ~2-3) is strong evidence
## against Poisson for the intensity-given-responsive component.

## ------------------------------------------------------------
## 2. Fit hurdle models: truncated_nbinom2 vs truncated_poisson,
##    compare via AIC to formally confirm NB is needed even after
##    conditioning on response
## ------------------------------------------------------------
fit_hurdle_nb <- glmmTMB(
  count ~ condition + (1 | condition:replicate),
  ziformula = ~ condition + (1 | condition:replicate),
  family = truncated_nbinom2,
  data = df
)

fit_hurdle_pois <- glmmTMB(
  count ~ condition + (1 | condition:replicate),
  ziformula = ~ condition + (1 | condition:replicate),
  family = truncated_poisson,
  data = df
)

AIC(fit_hurdle_nb, fit_hurdle_pois)   # expect NB to win by a wide margin given the var/mean ratios above

model <- fit_hurdle_nb
summary(model)

## ------------------------------------------------------------
## 2b. CONVERGENCE CHECK -- do this BEFORE trusting any emmeans
##     output. A response-scale estimate near 0 paired with an
##     astronomically wide CI (e.g. upper bound ~1e+289) is NOT a
##     real confidence interval -- it's a symptom of non-convergence,
##     most likely because the (1|condition:replicate) random effect
##     in the cond (intensity) formula is poorly identified for
##     conditions with very few nonzero cells (your 0h baselines).
## ------------------------------------------------------------

# Positive-definite Hessian is the basic convergence sanity check --
# FALSE means the optimizer did not find a proper solution
model$sdr$pdHess

# Look for "Model convergence problem" / "NA/NaN function evaluation"
# warnings when the model was fit, and check for absurdly large SEs
# in the cond part of the summary output above
summary(model)$coefficients$cond

# Check random-effect variance estimates -- a variance component that
# is exactly 0, or implausibly huge, both indicate estimation trouble
VarCorr(model)

## If pdHess is FALSE or you see huge SEs / degenerate variances,
## the fix is almost always to SIMPLIFY the random-effect structure
## for the cond (intensity) formula specifically, since that's the
## part being asked to estimate a per-condition-replicate variance
## from very sparse nonzero data at 0h. Try dropping the
## condition:replicate nesting to a simpler (1|replicate) term for
## the cond part (keeping the richer structure for ziformula, which
## has ALL cells -- including zeros -- to work with, so it's much
## less data-starved):

fit_hurdle_nb_simplified <- glmmTMB(
  count ~ condition + (1 | replicate),                       # simplified: no condition nesting here
  ziformula = ~ condition + (1 | condition:replicate),        # unchanged -- zi has plenty of data
  family = truncated_nbinom2,
  data = df
)

fit_hurdle_nb_simplified$sdr$pdHess   # check convergence again
AIC(fit_hurdle_nb, fit_hurdle_nb_simplified)

## If fit_hurdle_nb_simplified converges cleanly (pdHess = TRUE, no
## warnings, reasonable SEs) while the original doesn't, use the
## simplified version from here on -- a converged simpler model beats
## an unconverged more complex one every time. If BOTH still show
## convergence issues, the sparse 0h data may not support a per-
## replicate random effect on the intensity component at all; as a
## further fallback, try condition as the only intensity predictor
## with no random effect there, while keeping it in ziformula:
# fit_hurdle_nb_norand <- glmmTMB(
#   count ~ condition,
#   ziformula = ~ condition + (1 | condition:replicate),
#   family = truncated_nbinom2,
#   data = df
# )

model <- fit_hurdle_nb_simplified   # or fit_hurdle_nb_norand, whichever converges cleanly
summary(model)

## ------------------------------------------------------------
## 3. Diagnostics
## ------------------------------------------------------------
sim_res <- simulateResiduals(model, n = 1000)
plot(sim_res)
testDispersion(sim_res)
plotResiduals(sim_res, form = df$condition)

## ------------------------------------------------------------
## 4. Two separate, biologically meaningful comparisons
## ------------------------------------------------------------

## (a) PROBABILITY OF RESPONSE: does condition change the fraction
##     of cells that express any transcript at all? This pulls from
##     the hurdle (zero) component of the model.
emm_hurdle <- emmeans(model, ~ condition, component = "zi", type = "response")
print(emm_hurdle)   # response scale = probability of being a STRUCTURAL zero (non-responsive)
## Note: glmmTMB's zi component models P(non-response), so smaller
## values here mean a HIGHER fraction of cells respond. Flip sign/
## framing as needed when reporting ("1 - zi probability" = % responsive).

pairs(emm_hurdle, adjust = "tukey")


#Perform testing of hurdle conditions only to determine whether conditions cause activation vs 0h
emm_hurdle_rawpairs<-as.data.frame(pairs(emm_hurdle, adjust = "none"))
p.adjust(emm_hurdle_rawpairs %>% filter(contrast %in% c("Incision_24h / Agarose_24h","Incision_24h / Amputation_24h","Agarose_24h / Amputation_24h")) %>% pull(p.value), method = "bonferroni", n = 3)
activation.comparisons<-c("Incision_24h / Agarose_24h","Incision_24h / Amputation_24h","Agarose_24h / Amputation_24h")
emm_hurdle_activationtest<-p.adjust(emm_hurdle_rawpairs %>% filter(contrast %in% activation.comparisons) %>% pull(p.value), method = "BH", n = length(activation.comparisons))
emm_hurdle_activationtest<-as.data.frame(emm_hurdle_activationtest)
emm_hurdle_activationtest<-cbind(activation.comparisons,emm_hurdle_activationtest)
emm_hurdle_activationtest.ft<-flextable(emm_hurdle_activationtest)
save_as_image(emm_hurdle_activationtest.ft, path = paste0(dataset.name,".HurdleNB2.emm_ActivationOnly_response.png"),res=300)

#Try test of both activation and difference in activation at 6 and 12
#activation.comparisons2<-c("Incision_0h / Incision_6h","Amputation_0h / Amputation_6h","Incision_0h / Incision_12h","Amputation_0h / Amputation_12h", "Amputation_6h / Incision_6h", "Amputation_12h / Incision_12h")
#emm_hurdle_activationtest2<-p.adjust(emm_hurdle_rawpairs %>% filter(contrast %in% activation.comparisons2) %>% pull(p.value), method = "BH", n = length(activation.comparisons2))
#emm_hurdle_activationtest2<-as.data.frame(emm_hurdle_activationtest2)
#emm_hurdle_activationtest2<-cbind(activation.comparisons2,emm_hurdle_activationtest2)
#emm_hurdle_activationtest2.ft<-flextable(emm_hurdle_activationtest2)
#save_as_image(emm_hurdle_activationtest2.ft, path = paste0(dataset.name,".HurdleNB2.emm_ActivationAndCompare_response.png"),res=300)



## (b) INTENSITY GIVEN RESPONSE: among cells that DO respond, does
##     condition change how many transcripts they make? This pulls
##     from the truncated-NB count component (conditional mean).
## IMPORTANT FIX: component = "cond" gives the mean of the UNTRUNCATED
## NB distribution (mu), not the mean of cells conditional on being
## responsive. Per glmmTMB's own documentation: "In hurdle models,
## component = 'cmean' produces means of the truncated conditional
## distribution, while component = 'cond', type = 'response' produces
## means of the untruncated conditional distribution." Since this
## dataset is so zero-heavy, mu can be tiny even when the true
## per-responsive-cell mean (what your raw replicate averages and the
## violin plot show) is much larger -- component="cmean" is the one
## that matches the raw data.
emm_intensity <- emmeans(model, ~ condition, component = "cmean", type = "response")
print(emm_intensity)

pairs(emm_intensity, adjust = "tukey")

library(flextable)
library(magrittr)
emm_intensity_pairs_df<-as.data.frame(pairs(emm_intensity, adjust = "tukey"))
emm_intensity_pairs.ft <- flextable(emm_intensity_pairs_df) %>% set_caption(caption = "HurdleNB model intensity Tukey-adjusted pvalues")
save_as_image(emm_intensity_pairs.ft, path = paste0(dataset.name,".HurdleNB2.emm_intensity.png"),res=300)

emm_intensity.ft2 <- flextable(emm_intensity_pairs_df %>% filter(contrast %in% c("Incision_0h - Incision_6h","Amputation_0h - Amputation_6h","Amputation_6h - Incision_6h","Amputation_12h - Incision_12h")))
save_as_image(emm_intensity.ft2, path = paste0(dataset.name,".HurdleNB2.emm_intensity_filtered.png"),res=300)


#Perform testing of hurdle conditions only to determine whether conditions cause activation vs 0h
emm_intensity_rawpairs<-as.data.frame(pairs(emm_intensity, adjust = "none"))
#p.adjust(emm_intensity_rawpairs %>% filter(contrast %in% c("Incision_0h - Incision_6h","Amputation_0h - Amputation_6h","Amputation_6h - Incision_6h","Amputation_12h - Incision_12h")) %>% pull(p.value), method = "bonferroni", n = 4)
intensity.comparisons<-c("Incision_24h - Agarose_24h","Incision_24h - Amputation_24h","Agarose_24h - Amputation_24h")
emm_intensity_activationtest<-p.adjust(emm_intensity_rawpairs %>% filter(contrast %in% intensity.comparisons) %>% pull(p.value), method = "BH", n = length(intensity.comparisons))
emm_intensity_activationtest<-as.data.frame(emm_intensity_activationtest)
emm_intensity_activationtest<-cbind(activation.comparisons,emm_intensity_activationtest)
emm_intensity_activationtest.ft<-flextable(emm_intensity_activationtest)
save_as_image(emm_intensity_activationtest.ft, path = paste0(dataset.name,".HurdleNB2.emm_IntensityLimited_response.png"),res=300)




## ------------------------------------------------------------
## 5. Focused comparisons: Incision vs Amputation at 6h and 12h,
##    for BOTH components (consistent with the rest of this analysis)
## ------------------------------------------------------------
#targeted_pairs <- list(
#  c("Amputation_6h",  "Incision_6h"),
#  c("Amputation_12h", "Incision_12h")
#)

#for (p in targeted_pairs) {
#  cat("\n---", p[1], "vs", p[2], "---\n")
#  cat("Responsiveness (hurdle component):\n")
#  print(pairs(emm_hurdle, adjust = "none")[
#    grepl(paste(p, collapse = "|"), as.character(pairs(emm_hurdle, adjust = "none")))
#  ])
#  cat("Intensity given response (truncated-NB component):\n")
#  print(pairs(emm_intensity, adjust = "none")[
#    grepl(paste(p, collapse = "|"), as.character(pairs(emm_intensity, adjust = "none")))
#  ])
#}
## (If the string-matching filter above is finicky with your exact
## condition names, it's just as easy to read off the relevant rows
## directly from the full `pairs(emm_hurdle)` / `pairs(emm_intensity)`
## tables printed in section 4.)

## ------------------------------------------------------------
## 6. Visualization: split view of the two components
## ------------------------------------------------------------

## (a) Fraction of responsive cells per condition, with model-based
##     95% CIs (from the hurdle component's replicate-aware GLMM,
##     NOT a naive pooled-cell binomial CI -- consistent with
##     avoiding pseudoreplication throughout this analysis) and
##     individual replicate datapoints overlaid for raw visual context

hurdle_emm_df <- as.data.frame(emm_hurdle) %>%
  rename(p_nonresponse = response, LCL_nonresponse = asymp.LCL, UCL_nonresponse = asymp.UCL) %>%
  mutate(
    pct_responsive = 100 * (1 - p_nonresponse),
    # CI bounds flip when converting P(non-response) -> P(response)
    pct_responsive_LCL = 100 * (1 - UCL_nonresponse),
    pct_responsive_UCL = 100 * (1 - LCL_nonresponse),
    # SE is unaffected by the 1-p flip (derivative of 1-p w.r.t. p is
    # -1, so magnitude is unchanged); just rescale to the 0-100% axis
    pct_responsive_SE = 100 * SE
  )

hurdle_emm_df <- hurdle_emm_df %>% 
  arrange(factor(condition, levels = conditions.order))

# raw per-replicate proportions (n=3 points per condition), no modeling
rep_responsive <- df %>%
  group_by(condition, replicate) %>%
  summarise(pct_responsive = mean(count > 0) * 100, .groups = "drop")
rep_responsive <- rep_responsive %>% 
  arrange(factor(condition, levels = conditions.order))


p_responsive <- ggplot() +
  geom_col(data = hurdle_emm_df, aes(x = fct_relevel(condition, conditions.order), y = pct_responsive),
           fill = "steelblue", alpha = 0.55) +
  geom_errorbar(data = hurdle_emm_df,
                aes(x = fct_relevel(condition, conditions.order), ymin = pct_responsive_LCL, ymax = pct_responsive_UCL),
                width = 0.2, linewidth = 0.8, color = "black") +
  geom_jitter(data = rep_responsive, aes(x = condition, y = pct_responsive),
              width = 0.12, size = 2.8, shape = 21, color = "black", fill = "white", stroke = 1) +
  labs(
    title = "Percent of cells with detectable transcript",
    subtitle = "Bar + black CI = hurdle-model estimate (replicate-aware); open circles = individual replicates",
    x = NULL, y = "% responsive cells"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_responsive)
ggsave("HurdleNB2_p_responsive.png")

## ---- SE version: same plot, but with +/- 1 SE error bars instead
##      of the 95% CI. SE bars are narrower than the CI (SE is
##      roughly CI_width / (2*1.96) for a Wald-type CI, i.e. SE bars
##      cover a smaller, less conservative range) -- useful when you
##      want a less visually dramatic error bar for a figure, but
##      keep in mind SE bars do NOT directly correspond to a 95%
##      coverage statement the way the CI version does.

p_responsive_se <- ggplot() +
  geom_col(data = hurdle_emm_df, aes(x = fct_relevel(condition, conditions.order), y = pct_responsive),
           fill = "steelblue", alpha = 0.55) +
  geom_errorbar(data = hurdle_emm_df,
                aes(x = fct_relevel(condition, conditions.order),
                    ymin = pct_responsive - pct_responsive_SE,
                    ymax = pct_responsive + pct_responsive_SE),
                width = 0.2, linewidth = 0.8, color = "black") +
  geom_jitter(data = rep_responsive, aes(x = condition, y = pct_responsive),
              width = 0.12, size = 2.8, shape = 21, color = "black", fill = "white", stroke = 1) +
  labs(
    title = "Percent of cells with detectable transcript",
    subtitle = "Bar + black error bar = hurdle-model estimate +/- SE (replicate-aware); open circles = individual replicates",
    x = NULL, y = "% responsive cells"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_responsive_se)
ggsave("HurdleNB2_p_responsive_se.png")


hurdle_emm_df<-hurdle_emm_df %>% mutate(condition_copy=condition) %>% separate_wider_delim(condition_copy, delim ="_" ,names=c("type","time"))
hurdle_emm_df$type<-factor(hurdle_emm_df$type)

p_responsive_se <- ggplot() +
  geom_col(data = hurdle_emm_df, aes(x = factor(condition,levels=conditions.order), y = pct_responsive, fill=type), 
           alpha = 1) +
  scale_fill_manual(values = c("Incision" = "#4bb254", "Amputation" = "#1E90FF","Agarose"="orange"))+
  geom_errorbar(data = hurdle_emm_df,
                aes(x = condition,
                    ymin = pct_responsive - pct_responsive_SE,
                    ymax = pct_responsive + pct_responsive_SE),
                width = 0.2, linewidth = 0.8, color = "black") +
  geom_jitter(data = rep_responsive, aes(x = condition, y = pct_responsive),
              width = 0.12, size = 2.8, shape = 21, color = "black", fill = "white", stroke = 1) +
  labs(
    title = "Percent of cells with detectable transcript",
    subtitle = "Bar + black error bar = hurdle-model estimate +/- SE (replicate-aware); open circles = individual replicates",
    x = NULL, y = "% responsive cells"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  #  theme_minimal(base_size = 12) +
  #  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.1))+
  #  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 12, face = "bold"))+
  theme(axis.text.x = element_text(face = "bold", size = 12))+
  theme(axis.text.y = element_text(face = "bold", size = 12))+
  theme(strip.text = element_text(face = "bold"))+
  theme(legend.position = "none") 
print(p_responsive_se)
ggsave("HurdleNB2_p_responsive_se_color.png")


## Note: bars/CIs come from the model's hurdle component (proper
## replicate random-effect structure), not from pooling all ~5000
## cells per replicate into one binomial proportion -- pooling would
## make the CI artificially tight and repeat the same pseudo-
## replication issue flagged earlier in this analysis for the count
## data itself. If any replicate's open circle sits visibly outside
## the CI band, that's worth a second look (same logic as checking
## individual points against a GLMM estimate elsewhere in this script).

## (b) Distribution of counts AMONG RESPONSIVE CELLS ONLY (raw data,
##     visual complement to the intensity-component test) -- this is
##     the "remove zeros" view you were already looking at, now paired
##     with a model that handles it rigorously rather than an ad hoc fit
p_intensity <- df %>%
  filter(count > 0) %>%
  ggplot(aes(x = fct_relevel(condition, conditions.order), y = count)) +
  geom_violin(fill = "grey85", color = "grey60", scale = "width", trim = TRUE) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
  scale_y_continuous(trans = "log1p") +
  labs(title = "mRNA count among responsive cells only",
       x = NULL, y = "mRNA spots per cell (log1p scale)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_intensity)
ggsave("HurdleNB2_p_intensity.png")

## ------------------------------------------------------------
## 6b.5. Density version of the intensity plot: one density curve
##       per condition (log1p x-axis, same scale as p_intensity),
##       with replicate means overlaid as rug marks and the hurdle
##       model's estimated mean +/- 95% CI (from emm_intensity,
##       section 4) shown as a vertical line + shaded band.
## ------------------------------------------------------------

# Replicate-level means among responsive cells only (raw, no modeling)
rep_intensity_means <- df %>%
  filter(count > 0) %>%
  group_by(condition, replicate) %>%
  summarise(mean_count = mean(count), .groups = "drop")

rep_intensity_means <- rep_intensity_means %>% 
  arrange(factor(condition, levels = conditions.order))


# Model-estimated mean (response/raw-count scale) +/- 95% CI, from
# emm_intensity (already computed in section 4) -- these are on the
# same raw-count units as the x-axis, before the log1p display transform
intensity_emm_df <- as.data.frame(emm_intensity) %>%
  rename(model_mean = emmean, model_LCL = asymp.LCL, model_UCL = asymp.UCL)

intensity_emm_df <- intensity_emm_df %>% 
  arrange(factor(condition, levels = conditions.order))


p_intensity_density <- df %>%
  filter(count > 0) %>%
  ggplot(aes(x = count)) +
  geom_density(fill = "grey80", color = "grey40", alpha = 0.6) +
  # model estimate: shaded 95% CI band across the full panel height,
  # plus a solid line at the point estimate
  geom_rect(data = intensity_emm_df,
            aes(xmin = model_LCL, xmax = model_UCL, ymin = 0, ymax = Inf),
            inherit.aes = FALSE, fill = "firebrick", alpha = 0.15) +
  geom_vline(data = intensity_emm_df, aes(xintercept = model_mean),
             color = "firebrick", linewidth = 1) +
  # replicate means as rug marks along the bottom (n=3 per condition)
  geom_rug(data = rep_intensity_means, aes(x = mean_count),
           color = "black", linewidth = 1, sides = "b", length = unit(0.06, "npc")) +
  scale_x_continuous(trans = "log1p") +
  facet_wrap(~ fct_relevel(condition, c("Incision_0h","Incision_6h","Incision_12h","Amputation_0h","Amputation_6h","Amputation_12h")), ncol = 3, scales = "free_y") +
  # facet_wrap(~ fct_relevel(condition, c("Incision_0h","Incision_6h","Incision_12h","Amputation_0h","Amputation_6h","Amputation_12h")), ncol = 3, scales = "fixed") +
  labs(
    title = "Density of mRNA counts among responsive cells",
    subtitle = "Grey density = raw distribution; black ticks = replicate means; red line/band = hurdle-model estimate +/- 95% CI",
    x = "mRNA spots per cell (log1p scale)", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 12, face = "bold"))+
  theme(axis.text = element_text(face = "bold", size = 8))+
  theme(strip.text = element_text(face = "bold"))

print(p_intensity_density)
ggsave("HurdleNB2_p_intensity_density.png")

## Note: the model's estimated mean (red line) is the hurdle model's
## conditional mean for the truncated-NB component, which can differ
## noticeably from the raw density's visual peak/mode -- NB
## distributions are right-skewed, so the mean sits to the right of
## the mode, and the CI band reflects uncertainty in that mean
## specifically (driven by between-replicate variability via the
## random effect), not the spread of individual cells.

## ------------------------------------------------------------
## 6b.6. Violin-plot version of the same comparison: x = condition,
##       violin = per-cell density (nonzero cells only), open blue
##       circles = jittered replicate means, red point + error bar =
##       hurdle-model estimate +/- 95% CI. Same visual pattern used
##       earlier in this analysis for the raw-count GLMM plots.
## ------------------------------------------------------------

# Jittered x-positions for the replicate mean points, so the 3 points
# per condition don't overlap -- same jitter approach used earlier
# (nudged left of center, small random jitter within that offset)
rep_intensity_means <- rep_intensity_means %>%
  mutate(
    x_num = as.numeric(factor(condition, levels = levels(df$condition))),
    x_jit = x_num - 0.15 + runif(n(), min = -0.05, max = 0.05)
  )


rep_intensity_means <- rep_intensity_means %>% arrange(factor(condition, levels = conditions.order))
intensity_emm_df<-intensity_emm_df %>% arrange(factor(condition, levels = conditions.order))


#fct_relevel(condition, conditions.order)

p_intensity_violin <- ggplot() +
  geom_violin(data = df %>% filter(count > 0),
              aes(x = fct_relevel(condition, conditions.order), y = count),
              fill = "grey85", color = "grey60", scale = "width", trim = TRUE) +
  geom_point(data = rep_intensity_means,
             aes(x = fct_relevel(condition, conditions.order), y = mean_count),
             #             aes(x = x_jit, y = mean_count),
             color = "steelblue", size = 2.5, shape = 1, position=position_nudge(x=-0.15)) +
  geom_point(data = intensity_emm_df, aes(x = fct_relevel(condition, conditions.order), y = model_mean),
             color = "firebrick", size = 3, position = position_nudge(x = 0.15)) +
  geom_errorbar(data = intensity_emm_df,
                aes(x = fct_relevel(condition, conditions.order), ymin = model_LCL, ymax = model_UCL),
                width = 0.12, color = "firebrick", linewidth = 0.8,
                position = position_nudge(x = 0.15)) +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "mRNA count among responsive cells",
    subtitle = "Grey = per-cell distribution (nonzero cells); Blue = replicate means (n=3); Red = hurdle-model estimate \u00b1 95% CI",
    x = "Condition", y = "mRNA spots per cell (log1p scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size=11),
        plot.subtitle=element_text(size=10),
        axis.text.x = element_text(angle = 30, hjust = 1))+
  #  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.1))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 12, face = "bold"))+
  theme(axis.text.x = element_text(face = "bold", size = 12))+
  theme(axis.text.y = element_text(face = "bold", size = 12))+
  theme(strip.text = element_text(face = "bold"))

print(p_intensity_violin)
ggsave("HurdleNB2_p_intensity_violin.png")


p_intensity <- df %>%
  filter(count > 0) %>%
  ggplot(aes(x = factor(condition,levels=conditions.order), y = count,color=type)) +
  scale_color_manual(values = c("Incision" = "#4bb254", "Amputation" = "#1E90FF","Agarose"="orange"))+
  geom_violin(fill = "grey85", color = "grey60", scale = "width", trim = TRUE) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
  scale_y_continuous(trans = "log1p") +
  labs(title = "mRNA count among responsive cells only",
       x = NULL, y = "mRNA spots per cell (log1p scale)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(size = 12))+
  #  theme(axis.text = element_text(face = "bold", size = 12))+
  theme(axis.text.x = element_text(face = "bold", size = 13))+
  theme(axis.text.y = element_text(face = "bold", size = 12))+
  theme(legend.position = "none")

print(p_intensity)
ggsave("HurdleNB2_p_intensity1.png")



p_intensity <- df %>%
  filter(count > 0) %>%
  ggplot(aes(x = factor(condition,levels=conditions.order), y = count,color=type)) +
  scale_color_manual(values = c("Incision" = "#4bb254", "Amputation" = "#1E90FF","Agarose"="orange"))+
  geom_violin(fill = "grey85", color = "grey60", scale = "width", trim = TRUE) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
  scale_y_continuous(trans = "log1p") +
  labs(title = "mRNA count among responsive cells only",
       x = NULL, y = "mRNA spots per cell (log1p scale)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(size = 12))+
  #  theme(axis.text = element_text(face = "bold", size = 12))+
  theme(axis.text.x = element_text(face = "bold", size = 13))+
  theme(axis.text.y = element_text(face = "bold", size = 15))+
  theme(legend.position = "none")

print(p_intensity)
ggsave("HurdleNB2_p_intensity2.png",width=4, height=8)



## Note: x-axis order follows levels(df$condition) -- if you want a
## specific order (e.g. Incision_0h, Amputation_0h, Incision_6h, ...)
## rather than whatever order the factor currently has, set it
## explicitly beforehand with:
##   df$condition <- factor(df$condition, levels = conditions.order)
## and rerun this section so rep_intensity_means' x_num/x_jit and the
## violin's x-axis both follow the same order.

## ------------------------------------------------------------
## 6c. P-VALUE-ANNOTATED versions of the intensity and responsiveness
##     plots, for just the two comparisons of interest:
##     Incision_6h vs Amputation_6h, Incision_12h vs Amputation_12h.
##
##     P-values come from the SAME hurdle-model contrasts already
##     computed in section 4 (emm_hurdle for responsiveness,
##     emm_intensity for intensity-given-response) -- not a separate
##     ad hoc test -- so the brackets shown here are consistent with
##     everything else in this analysis.
## ------------------------------------------------------------

# install.packages("ggpubr")
library(ggpubr)

comparisons <- list(
  c("Incision_6h",  "Amputation_6h"),
  c("Incision_12h", "Amputation_12h")
)

## Helper: pull a p-value out of an emmeans pairwise-contrast data
## frame (e.g. from as.data.frame(pairs(emm_obj, adjust = "tukey")))
## for a specific pair of conditions.
##
## IMPORTANT: when emmeans() is called with type = "response" on a
## log-link model (as emm_intensity and emm_hurdle both are), pairs()
## produces RATIO contrasts, formatted as "A / B" -- not "A - B".
## This must match that separator or nothing will be found.
get_emm_contrast_p <- function(pw_df, condA, condB) {
  contrast_str <- as.character(pw_df$contrast)
  parts <- strsplit(contrast_str, " / ", fixed = TRUE)
  match_idx <- vapply(parts, function(x) setequal(trimws(x), c(condA, condB)), logical(1))
  if (!any(match_idx)) stop("Contrast not found: ", condA, " vs ", condB,
                            ". Available contrasts: ", paste(contrast_str, collapse = "; "))
  pw_df$p.value[which(match_idx)[1]]
}

## Defensive parser: emmeans sometimes formats very small p-values as
## the string "<0.0001" rather than a plain number (this can persist
## even after as.data.frame(), depending on emmeans version). Handle
## both cases so star-thresholding and the printed label always work.
format_p_label <- function(p) {
  if (is.character(p) && grepl("^<", p)) {
    numeric_approx <- as.numeric(sub("^<", "", p))  # e.g. 0.0001, used only for star thresholding
    display <- paste0("p ", p)                       # e.g. "p <0.0001"
  } else {
    numeric_approx <- as.numeric(p)
    display <- paste0("p = ", formatC(numeric_approx, format = "e", digits = 2))
  }
  stars <- dplyr::case_when(
    numeric_approx < 0.0001 ~ "****",
    numeric_approx < 0.001  ~ "***",
    numeric_approx < 0.01   ~ "**",
    numeric_approx < 0.05   ~ "*",
    TRUE                    ~ "ns"
  )
  paste0(stars, "\n", display)
}

build_pval_df <- function(emm_obj, comparisons, adjust = "tukey") {
  pw_df <- as.data.frame(pairs(emm_obj, adjust = adjust))
  tibble::tibble(
    group1 = vapply(comparisons, `[`, character(1), 1),
    group2 = vapply(comparisons, `[`, character(1), 2)
  ) %>%
    rowwise() %>%
    mutate(p_value_raw = list(get_emm_contrast_p(pw_df, group1, group2))) %>%
    ungroup() %>%
    mutate(
      p_value = vapply(p_value_raw, function(p) if (is.character(p) && grepl("^<", p))
        as.numeric(sub("^<", "", p)) else as.numeric(p), numeric(1)),
      label = vapply(p_value_raw, format_p_label, character(1))
    ) %>%
    select(-p_value_raw)
}

