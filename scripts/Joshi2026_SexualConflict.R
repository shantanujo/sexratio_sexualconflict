setwd("C:/Users/shant/Documents/SexualConflict_Manuscript/latest/final_ms")

library(data.table)
library(tidyverse)
library(lme4)
library(glmmTMB)
library(emmeans)
library(ggeffects)
setwd("~/data")
##########       Sex  Ratio    ##########
newcount = fread("sex_ratio_counts.csv")
agg <- newcount %>%
  group_by(sp) %>%
  summarise(   males   = sum(males),    females = sum(females),    total   = males + females)
agg
mat <- matrix(
  c(agg$males[1], agg$females[1], agg$males[2],   agg$females[2]),
  nrow = 2,
  byrow = TRUE
)
rownames(mat) <- levels(newcount$sp)
colnames(mat) <- c("Male", "Female")

res_fisher <- fisher.test(mat, alternative = "less")  # Run Fisher's exact test (one‐sided)
res_fisher
res_prop <- with(agg,    # Perform two‐sample proportion test
                 prop.test(
                   x = males, n = total, alternative = "less", correct     = FALSE  ))
datemod= glmer(data=newcount, cbind(males,females) ~ sp *(1|date_j), family=binomial)
exp(confint(datemod, parm="beta_", method="profile"))
newcount %>%group_by(sp) %>%    # adding columns for Fig 1b
  summarise(n       = n(),
            mean    = mean(ratio, na.rm = TRUE),
            sd      = sd(ratio, na.rm = TRUE),
            se      = sd / sqrt(n),
            ci_lower = mean - qt(0.975, df = n - 1) * se,
            ci_upper = mean + qt(0.975, df = n - 1) * se ) 


##########       SINGLE  SPECIES TRIALS   ##########
caged <- fread("single_species_trials.csv") %>%
  mutate(
    sp      = factor(sp),no_male = factor(no_male)  )

# Drop unused levels within species
caged <- caged %>% group_by(sp) %>%mutate(no_male = droplevels(no_male))%>% ungroup()
levels(droplevels(caged$no_male[caged$sp == "TRA"]))  # should be "1" "4"

#####  Survival   ######
mod_surv <- glmer(female_surv ~ no_male * sp + (1 | date_j),
                  data    = caged, family  = binomial,control = glmerControl(optimizer= "bobyqa",
                    optCtrl= list(maxfun = 2e5)
                  ))
mod_surv1 <- glmer(female_surv ~ no_male * sp +mat_th+ (1 | date_j),
                  data    = caged,
                  family  = binomial,
                  control = glmerControl(
                    optimizer        = "bobyqa",
                    optCtrl          = list(maxfun = 2e5)
                  ))
anova(mod_surv, mod_surv1, test = "Chisq")

raw_surv <- caged %>%
  group_by(no_male, sp) %>%
  summarise(n_total = n(),
    n_survived = sum(female_surv == 1),
    prop_surv  = n_survived / n_total) %>%ungroup()
emm <- emmeans(mod_surv, ~ sp | no_male,type = "response")
emm_df <- as.data.frame(emm) %>%
  rename( prop_surv = prob,
    lower     = asymp.LCL,
    upper     = asymp.UCL)

######     Fecundity    #######
cage2 <- subset(caged, female_surv== '1') #keeping only females that survived
eggmod <- glmmTMB(
  eggs ~ no_male*sp + (1|date_j),
  family = nbinom2,ziformula = ~sp,
  data   = cage2)
summary(eggmod)$coefficients$cond
summary(eggmod)$coefficients$zi
#comparing models with and without mating attempts
eggmod1 <- glmmTMB(
  eggs ~ no_male*sp + mat_th+(1|date_j),
  family = nbinom2,ziformula = ~sp,
  data   = cage2)
anova(eggmod, eggmod1, test = "Chisq")

preds <- ggpredict(eggmod, terms = c( "no_male", "sp"))
# Filter out the unwanted combination
preds_filtered <- subset(preds, !(x == "6" & group == "TRA"))
colnames(preds_filtered)[1] = "no_male"
preds_filtered$sp = preds_filtered$group

##### Fertility   ###
cage3 <- subset(caged, eggs >0) # for fertility analysis

cage3$failures = cage3$eggs - cage3$larvae

larmod <- glmer(cbind(larvae, failures) ~ no_male * sp + (1 | date_j), data = cage3, family = binomial)
larmoda <- glmer(cbind(larvae, failures) ~ no_male * sp +mat_th + (1 | date_j), data = cage3, family = binomial)
anova(larmoda, larmod, test="LRT")
summary(larmod)
preds1 <- ggpredict(larmod, terms = c( "no_male", "sp"))
# Filter out the unsupported combination
preds_filtered1 <- subset(preds1, !(x == "6" & group == "TRA"))
colnames(preds_filtered1)[1] = "no_male"
preds_filtered1$sp = preds_filtered$group

###########     'MIXED'   SPECIES   ############
mix = fread("mixed_species_trials.csv") 
mix <- mix %>%
  mutate(  female_sp    = factor(female_sp),
    con_male =factor(con_male),
    ratio =factor(ratio2) )
#Survival
mod_surv2 <- glmer(female_surv ~ ratio * female_sp + (1 | date_j),
                data    = mix,
                family  = binomial,
                control = glmerControl(
                  optimizer        = "bobyqa",
                  optCtrl          = list(maxfun = 2e5)
                ))
emm2 <- emmeans(mod_surv2, ~ female_sp | ratio,type = "response")
summary(mod_surv2)
# convert to data frame for ggplot
emm_df2 <- as.data.frame(emm2) %>%
  rename(  prop_surv = prob,
    lower     = asymp.LCL,
    upper     = asymp.UCL)

#### FECUNDITY ####
mix2 <- subset(mix, female_surv== '1')
eggmod_mix <- glmmTMB(
  eggs ~ ratio*female_sp + (1|date_j),
  family = nbinom2,ziformula = ~female_sp,
  data   = mix2)
summary(eggmod_mix)$coefficients$cond
summary(eggmod_mix)$coefficients$zi
predsa <- ggpredict(eggmod_mix, terms = c("ratio", "female_sp"))
colnames(predsa)[1] = "ratio"
predsa$female_sp = predsa$group

#### FERTILIITY ####
mix3 <- subset(mix, eggs > 0) 
mix3$failures = mix3$eggs - mix3$larvae
larmod1 <- glmer(cbind(larvae, failures) ~ ratio * female_sp + (1 | date_j), data = mix3, family = binomial)
larmod2 <- glmer(cbind(larvae, failures) ~ ratio * female_sp+tot_th + (1 | date_j), data = mix3, family = binomial)
anova(larmod1, larmod2)
summary(larmod1)
predsb <- ggpredict(larmod1, terms = c("ratio", "female_sp"))
colnames(predsb)[1] = "ratio"
predsb$female_sp = predsb$group

emm <- emmeans(larmod, ~ no_male*sp)
pairs(emm)

###### FIGURES  #####
species_palette <- c("EXS" = "#cc9200", "TRA" = "#9200cc")

# Fig 1a
ggplot(newcount, aes(x = sp, y = ratio)) +
  geom_jitter(width = 0.08, alpha = 0.6, color="gray50")+
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0,
    linewidth = 1.2,
    aes(color = sp)
  )+stat_summary(fun = mean, geom = "point",
                 size = 6.4, aes(color = sp), shape=18)  +
  labs( x = "Species",  y = "Sex Ratio (M/F)") +
  geom_hline(
    yintercept = 1,linetype   = "dashed",
    color = "gray50" ) +
  theme_classic(base_size  =14) +
  scale_y_continuous(limits = c(0, 7), breaks = seq(0, 7, by = 1))+
  scale_color_manual(values = species_palette) +guides(color = "none")

### Fig. 1b
ggplot(newcount, aes(x=date_j, y=ratio, col=sp)) + geom_jitter(width = 0.1, alpha = 0.5, size=2.2) +
  geom_smooth(method="loess", span=0.9, linewidth=1.4) +
  scale_color_manual(values = species_palette) + 
  theme_classic(base_size=14) +labs(color="Species")+
  xlab("Day of Year") + ylab("Sex Ratio (M/F)") +
  guides(color = "none")

###Fig 2a
# dodge width for overlapping points
pd <- position_dodge(width = 0.6)
fig2a=ggplot()+
  geom_point(data    = emm_df,
    aes(x = factor(no_male), y = prop_surv, color = sp),
    shape    = 18,
    size     = 8,
    position = pd ) + # error bars for emmeans CIs
  geom_errorbar(
    data    = emm_df,
    aes(x = factor(no_male), ymin = lower, ymax = upper, color = sp),
    width    = 0,linewidth=1.2,
    position = pd) +
  scale_y_continuous(
    name    = "Survival Probability",
    limits  = c(0,1),
    breaks  = seq(0,1,0.2)) +
  xlab("Number of Males") +
  labs(color = "Species") +
  theme_classic(base_size = 14) + scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color = "none")
fig2a

#Fig2b
fig2b=ggplot()+
  geom_point( data    = emm_df2,
    aes(x = factor(ratio), y = prop_surv, color = female_sp),
    shape    = 18,
    size     = 8,
    position = pd) + # error bars for emmeans CIs
  geom_errorbar(
    data    = emm_df2,
    aes(x = factor(ratio), ymin = lower, ymax = upper, color = female_sp),
    width    = 0,linewidth=1.2,
    position = pd) +
  scale_y_continuous(
    name    = "Survival Probability",
    limits  = c(0,1),
    breaks  = seq(0,1,0.2) ) +
  xlab("Conspecific: Heterospecific Male Ratio") +
  labs(color = "Species") +
  theme_classic(base_size = 14) + scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color = "none")
fig2b

#Fig 3a
fig3a <- ggplot() +
  geom_jitter(
    data = cage2,
    aes(x = no_male, y = eggs),
    width = 0.15, size = 2, alpha = 0.2  ) +
  geom_point(
    data = preds_filtered,
    aes(x = no_male, y = predicted, color = factor(sp)),
    position = position_dodge(width = 0.4),
    size = 7, shape = 18
  ) +
  geom_errorbar(
    data = preds_filtered,
    aes(x = no_male, ymin = conf.low, ymax = conf.high, color = factor(sp)),
    position = position_dodge(width = 0.5),
    width = 0,
    linewidth = 1.3) +
  facet_grid(. ~ sp, scales = "free_x", space = "free_x") +
  labs(x = "No. of Males", y = "No. of Eggs") +
  theme_classic(base_size = 12) + scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color="none")
fig3a


# Fig 3b
fig3b= ggplot() + # Jittered raw data in the background
  geom_jitter(data = mix2, aes(x = ratio, y = eggs, color=factor(female_sp)), 
              width = 0.15, colour = "black", size = 2, alpha = 0.2) +
  # Prediction points and error bars in the foreground
  geom_point(data = predsa, 
             aes(x = ratio, y = predicted, color = factor(female_sp)), 
             position = position_dodge(width = 0.4), size = 7, shape=18)  + 
  geom_errorbar(data = predsa, 
                aes(x = ratio, ymin = conf.low, ymax = conf.high, color = factor(female_sp)), 
                position = position_dodge(width = 0.5), width = 0, linewidth=1.4) +
  labs(x = "Conspecific: Heterospecific Male Ratio", y = "No. of Eggs", color = "") +
  theme_classic(base_size = 12) +
  facet_wrap(~female_sp)+ scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color="none")
fig3b

# Fig 3c
fig3c <- ggplot() +
  geom_jitter(
    data = cage3,
    aes(x = no_male, y = ferti),
    width = 0.15, size = 2, alpha = 0.2 ) +
  geom_point(
    data = preds_filtered1,
    aes(x = no_male, y = predicted, color = factor(sp)),
    position = position_dodge(width = 0.5),
    size = 7, shape = 18) +
  geom_errorbar(
    data = preds_filtered1,
    aes(x = no_male, ymin = conf.low, ymax = conf.high, color = factor(sp)),
    position = position_dodge(width = 0.5),
    width = 0,
    linewidth = 1.4 ) +
  facet_grid(. ~ sp, scales = "free_x", space = "free_x") +
  labs(x = "No. of Males", y = "Fertility") +
  theme_classic(base_size = 12) + scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color="none")
fig3c

# Fig 3d
fig3d= ggplot() +
  # Jittered raw data in the background
  geom_jitter(data = mix3, aes(x = ratio, y = ferti), 
              width = 0.15, colour = "black", size = 2, alpha = 0.2) +
  # Prediction points and error bars in the foreground
  geom_point(data = predsb, 
             aes(x = ratio, y = predicted, color = factor(female_sp)), 
             position = position_dodge(width = 0.5), size = 7, shape=18)  + 
  geom_errorbar(data = predsb, 
                aes(x = ratio, ymin = conf.low, ymax = conf.high, color = factor(female_sp)), 
                position = position_dodge(width = 0.5), width = 0, linewidth=1.4) +
  labs(x = "Conspecific: Heterospecific Male Ratio", y = "Fertility", color = "") +
  theme_classic(base_size = 12) +
  facet_wrap(~female_sp)+ scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color="none")
fig3d


##### SUPPLEMENTARY FIGURE  #####
ten=  read.csv("~/SexualConflict_Manuscript/latest/final_ms/supplementary/tenerals_2024.csv")
levels(ten$Sp)= c("EXS", "TRA")
ten <- ten %>%
  mutate(  Sp = recode(Sp,    "ENEX" = "EXS",
                       "ENTR" = "TRA"))
ten2 <- ten %>%
  mutate(
    Sex = str_trim(Sex)   # remove leading/trailing spaces
  ) %>%
  group_by(Date, Sp, Sex) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  pivot_wider(
    names_from  = Sex,
    values_from = Count,
    values_fill = 0
  ) %>%
  mutate(ratio = M / F )

species_palette <- c("EXS" = "#cc9200", "TRA" = "#9200cc")
ggplot(ten2, aes(x=Sp, y=ratio, col=Sp)) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", size=1.7, linewidth=1.2, shape=18) +
  xlab("Species") + ylab("Immature Sex Ratio") + ylim(0,2)+
  theme(text = element_text(size=22)) + scale_color_brewer(palette="Dark2") + theme_classic()+
  geom_jitter(alpha=0.4, width=0.1)+
  scale_color_manual(values = species_palette)+
  guides(color = "none")
