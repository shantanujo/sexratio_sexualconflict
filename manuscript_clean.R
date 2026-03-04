setwd("C:/Users/shant/Documents/R/cage_final")
library(ggeffects)
library(tidyverse)
library(dplyr)
library(lme4)
library(DHARMa)
library(sjPlot)
library(mgcv)
library(sjPlot)
library(mgcViz)
library(glmmTMB)
library(broom)
library(car)    
library(emmeans) 
library(multcomp)

##########       SINGLE  SPECIES   ##########
##sex ratio counts
newcount = readRDS("finalcount.RDS")
colnames(newcount)[13]= "ratio"
## cage trial data, mix= both species
#cagef <- read.csv("singlesp_25.csv")
#cagef$ferti <- cagef$larvae / cagef$eggs
caged = readRDS("singlesp.RDS")
caged <- caged %>%
  mutate(
    sp       = factor(sp),
    date_cat = factor(date_cat),
    obs      = factor(obs),
    no_male =factor(no_male),
    date_cat = factor(date_cat)
  )
cagea <- subset(caged, !eggs== '0')
cage <- subset(caged, female_dead== 'N')

### Mix Species ###
mix = readRDS("mixsp.RDS")
mix <- mix %>%
  mutate(
    female_sp       = factor(female_sp),
    obs      = factor(obs),
    con_male =factor(con_male)
  )
mix$tot_th = mix$th_con + mix$th_het



######  Sex Ratio   ####

#Fishers exact
# Create contingency table (rows = species, cols = sex)
library(dplyr)
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

# Run Fisher's exact test (one‐sided)
res_fisher <- fisher.test(mat, alternative = "less")
res_fisher
# Perform two‐sample proportion test
res_prop <- with(agg,
                 prop.test(
                   x = males,
                   n = total,
                   alternative = "less",
                   correct     = FALSE  ))
res_prop

##Sex ratio final model
datemod= glmer(data=newcount, cbind(males,females) ~ sp *(1|date_j), family=binomial)
summary(datemod)
# profile likelihood CIs
exp(confint(datemod, parm="beta_", method="profile"))



#####  SURVIVAL ANALYSIS   ######
mod_surv <- glmer(female_surv ~ no_male * sp + (1 | date_j),
                data    = caged,
                family  = binomial,
                control = glmerControl(
                  optimizer        = "bobyqa",
                  optCtrl          = list(maxfun = 2e5)
                ))

# assume fem_death: 1 = died, 0 = survived
raw_surv <- caged %>%
  group_by(no_male, sp) %>%
  summarise(
    n_total    = n(),
    n_survived = sum(female_surv == 1),
    prop_surv  = n_survived / n_total
  ) %>%
  ungroup()
emm <- emmeans(mod_surv,
               ~ sp | no_male,
               type = "response")

# convert to data frame for ggplot
emm_df <- as.data.frame(emm) %>%
  rename(
    prop_surv = prob,
    lower     = asymp.LCL,
    upper     = asymp.UCL
  )

######     FECUNDITY    #######
eggmod <- glmmTMB(
  eggs ~ sp*no_male + (1|date_j),
  family = nbinom2,ziformula = ~no_male*sp,
  data   = caged)

# Generate prediction data
preds <- ggpredict(eggmod, terms = c( "no_male", "sp"))
# Filter out the unsupported combination
preds_filtered <- subset(preds, !(x == "6" & group == "TRA"))
colnames(preds_filtered)[1] = "no_male"
# Plot: point-range grouped by species, color by no_male
caged$group=caged$sp
preds_filtered$sp = preds_filtered$group

cagedf <- caged %>%
  filter(!(sp == "TRA" & no_male == 6))
preds_filtered <- preds_filtered %>%
  filter(!(sp == "TRA" & no_male == 6))



##### FERTILITY   ###
caged2=cagea
lar = caged2$larvae
fail = caged2$failures
## larvae Conspecific
larmod <- glmmTMB(
  cbind(larvae, failures) ~ no_male * sp + (1 | date_j),
  family = betabinomial(),
  data = caged2
)
preds1 <- ggpredict(larmod, terms = c("no_male", "sp"))
# Filter out unsupported combinations (e.g., TRA × 6)
preds_filtered1 <- subset(preds1, !(group == "TRA" & x == "6"))
caged2$group=caged2$sp
preds_filtered1$sp = preds_filtered1$group



###########     'MIXED'   SPECIES   ############

mix$sp <- as.factor(mix$sp)
mix$no_male <- as.factor(mix$no_male)
mix$female_sp <- as.factor(mix$female_sp)
mix$con_male <- as.factor(mix$con_male)
mix$ferti <- mix$larvae / mix$eggs
mix$date= as.Date(mix$date,format= "%m/%d/%Y")
mix$date_j = as.POSIXlt(mix$date,format = "%d%b%y")$yday
mix[is.na(mix$eggs)] <- 0
mix[is.na(mix$larvae)] <- 0
levels(mix$female_sp) <- c("EXS", "TRA")
mix <- mix |>
  mutate(ratio2 = factor(con_male, levels = c(1, 3), labels = c("1:3", "3:1")))

mix2 <- subset(mix, !eggs== '0')
mix3 <- subset(mix, female_d== '1')

god_s2 <- glmer(female_d ~ ratio2 * female_sp + (1 | date_j),
                data    = mix,
                family  = binomial,
                control = glmerControl(
                  optimizer        = "bobyqa",
                  optCtrl          = list(maxfun = 2e5)
                ))
emm2 <- emmeans(god_s2,
                ~ female_sp | ratio2,
                type = "response")

# convert to data frame for ggplot
emm_df2 <- as.data.frame(emm2) %>%
  rename(
    prop_surv = prob,
    lower     = asymp.LCL,
    upper     = asymp.UCL
  )
pd <- position_dodge(width = 0.6)



#### FECUNDITY ####
eggmod_nb_mix <- glmmTMB(
  eggs ~ ratio2*female_sp + (1|date_j),
  family = nbinom2,ziformula = ~con_male*female_sp,
  data   = mix)

predsa <- ggpredict(eggmod_nb_mix, terms = c("ratio2", "female_sp"))
predsa$female_sp = predsa$group

#### FERTILIITY ####
mix2 =filter(mix, eggs>0)
mix2$failures = mix2$eggs - mix2$larvae
larmod1 <- glmmTMB(
  cbind(larvae, failures) ~ ratio2 * female_sp + (1 | date_j),
  family = betabinomial(),
  data = mix2
)
preds2 <- ggpredict(
  larmod1,
  terms = c("ratio2", "female_sp"),
  type = "fixed"
)
preds2$sp <- preds2$group


###### FIGURES  #####

# FIG 1 a  -  ratio plot
species_palette <- c("EXS" = "#cc9200", "TRA" = "#9200cc")

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

### Fig. 1B
#summary stats
newcount %>%
  group_by(sp) %>%
  summarise(
    n       = n(),
    mean    = mean(ratio, na.rm = TRUE),
    sd      = sd(ratio, na.rm = TRUE),
    se      = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se
  ) 

ggplot(newcount, aes(x=date_j, y=ratio, col=sp)) + geom_jitter(width = 0.1, alpha = 0.5, size=2.2) +
  geom_smooth(method="loess", span=0.9, linewidth=1.4) +
  scale_color_manual(values = species_palette) + 
  theme_classic(base_size=14) +labs(color="Species")+
  xlab("Day of Year") + ylab("Sex Ratio (M/F)") +
  guides(color = "none")

###Fig 2a
# dodge width for side-by-side points
pd <- position_dodge(width = 0.6)

fig2a=ggplot()+
  geom_point(
    data    = emm_df,
    aes(x = factor(no_male), y = prop_surv, color = sp),
    shape    = 18,
    size     = 8,
    position = pd
  ) +
  # error bars for emmeans CIs
  geom_errorbar(
    data    = emm_df,
    aes(x = factor(no_male), ymin = lower, ymax = upper, color = sp),
    width    = 0,linewidth=1.2,
    position = pd
  ) +
  scale_y_continuous(
    name    = "Survival Probability",
    limits  = c(0,1),
    breaks  = seq(0,1,0.2)
  ) +
  xlab("Number of Males") +
  labs(color = "Species") +
  theme_classic(base_size = 14) + scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color = "none")
fig2a

#Fig2b
fig2b=ggplot()+
  geom_point(
    data    = emm_df2,
    aes(x = factor(ratio2), y = prop_surv, color = female_sp),
    shape    = 18,
    size     = 8,
    position = pd
  ) +
  # error bars for emmeans CIs
  geom_errorbar(
    data    = emm_df2,
    aes(x = factor(ratio2), ymin = lower, ymax = upper, color = female_sp),
    width    = 0,linewidth=1.2,
    position = pd
  ) +
  scale_y_continuous(
    name    = "Survival Probability",
    limits  = c(0,1),
    breaks  = seq(0,1,0.2)
  ) +
  xlab("Conspecific: Heterospecific Male Ratio") +
  labs(color = "Species") +
  theme_classic(base_size = 14) + scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color = "none")
fig2b

#Fig 3a
fig3a= ggplot() +
  # Jittered raw data in the background
  geom_jitter(data = cagedf, aes(x = no_male, y = eggs), 
              width = 0.15, colour = "black", size = 2, alpha = 0.2) +
  # Prediction points and error bars in the foreground
  geom_point(data = preds_filtered, 
             aes(x = no_male, y = predicted, color=sp), 
             position = position_dodge(width = 0.5), size = 8, shape=18)  + 
  geom_errorbar(
    data = preds_filtered,
    aes(x = no_male, ymin = conf.low, ymax = conf.high, color = sp),
    position = position_dodge(width = 0.5),
    width = 0, linewidth = 1.4)+
  labs(x = "No. of Males", y = "No. of Eggs") +
  guides(color = "none")+
  theme_classic(base_size = 13) +
  facet_wrap(~group, scales="free_x") + scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color="none")
fig3a

# Fig 3b
fig3b= ggplot() +
  # Jittered raw data in the background
  geom_jitter(data = mix, aes(x = ratio2, y = eggs), 
              width = 0.15, colour = "black", size = 2, alpha = 0.2) +
  # Prediction points and error bars in the foreground
  geom_point(data = predsa, 
             aes(x = x, y = predicted, color = factor(female_sp)), 
             position = position_dodge(width = 0.4), size = 8, shape=18)  + 
  geom_errorbar(data = predsa, 
                aes(x = x, ymin = conf.low, ymax = conf.high, color = factor(female_sp)), 
                position = position_dodge(width = 0.5), width = 0, linewidth=1.5) +
  labs(x = "Conspecific: Heterospecific Male Ratio", y = "No. of Eggs", color = "") +
  theme_classic(base_size = 13) +
  facet_wrap(~group)+ scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color="none")
fig3b
# Fig 3c
fig3c=ggplot() +
  geom_jitter(data = caged2, aes(x = no_male, y = ferti), 
              width = 0.15, colour = "black", size = 2, alpha = 0.2) +
  geom_point(data = preds_filtered1, 
             aes(x = x, y = predicted, color=factor(sp)), 
             position = position_dodge(width = 0.5), size = 8, shape=18) +
  geom_errorbar(
    data = preds_filtered1,
    aes(x = x, ymin = conf.low, ymax = conf.high, color = factor(sp)),
    position = position_dodge(width = 0.5),
    width = 0,
    linewidth = 1.3 )+
  labs(x = "No. of Males", y = "Fertility", color = "") +
  theme_classic(base_size = 13) +
  facet_wrap(~group, scales="free_x") +
  scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color="none")
fig3c
# Fig 3d
fig3d= ggplot() +
  # Jittered raw data in the background
  geom_jitter(data = mix2, aes(x = ratio2, y = ferti), 
              width = 0.15, colour = "black", size = 2, alpha = 0.2) +
  # Prediction points and error bars in the foreground
  geom_point(data = preds2, 
             aes(x = x, y = predicted, color = factor(sp)), 
             position = position_dodge(width = 0.5), size = 8, shape=18)  + 
  geom_errorbar(data = preds2, 
                aes(x = x, ymin = conf.low, ymax = conf.high, color = factor(sp)), 
                position = position_dodge(width = 0.5), width = 0, linewidth=1.4) +
  labs(x = "Conspecific: Heterospecific Male Ratio", y = "Fertility", color = "") +
  theme_classic(base_size = 13) +
  facet_wrap(~group)+ scale_color_manual(values = species_palette)+
  theme(strip.background = element_blank())+
  guides(color="none")
fig3d

###      ###SUPPLEMENTARY FIGURE     ###
ten=  read.csv("tenerals_2024.csv")
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
  mutate(
    ratio = M / F
  )

species_palette <- c("EXS" = "#cc9200", "TRA" = "#9200cc")
ggplot(ten2, aes(x=Sp, y=ratio, col=Sp)) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", size=1.5, linewidth=1.2) +
  xlab("Species") + ylab("Immature Sex Ratio") + ylim(0,2)+
  theme(text = element_text(size=22)) + scale_color_brewer(palette="Dark2") + theme_classic()+
  geom_jitter(alpha=0.4, width=0.1)+
  scale_color_manual(values = species_palette)+
  guides(color = "none")
