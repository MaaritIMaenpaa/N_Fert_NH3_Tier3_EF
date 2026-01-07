# Packages ----
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(interactions)
library(ggpubr)
library(car)
library(glmmTMB)
library(DHARMa)
library(tidyverse)
library(performance)

#
#
#

# This is repeat code from original script ----

# Data ----
load(here::here("data/Wrangled/NH3predict.Rda"))

# Remove potential NAs from base-model variables
NH3predict <- NH3predict %>% drop_na(c(Fertiliser, Application, Method))

NH3predict <- NH3predict[which(is.na(NH3predict$NH3loss)==FALSE &
                                 is.na(NH3predict$Fertiliser)==FALSE & 
                                 is.na(NH3predict$Application)==FALSE &
                                 is.na(NH3predict$Method)==FALSE &
                                 is.na(NH3predict$Ref)==FALSE),]

# Removing a study with large negative values for NH3loss
NH3predict <- NH3predict[which(NH3predict$Ref!="Velthof et al. 1990"),]

#
#
# 

## DEFINING VARIABLES ----


# Defining class of vectors in the analysis

#
# Model 1 (A1)

NH3predict$Fertiliser <- factor(NH3predict$Fertiliser)
NH3predict$Place <- factor(NH3predict$Place)
NH3predict$Method <- factor(NH3predict$Method)
NH3predict$Application <- factor(NH3predict$Application)
NH3predict$Ref <- factor(NH3predict$Ref)

#
# Model 2 (A2)
NH3predict$Application.rate <- as.numeric(NH3predict$Application.rate)

#
# Model 3 (A3)
NH3predict$SoilpH <- as.numeric(NH3predict$SoilpH)
summary(NH3predict$SoilpH)
NH3predict[which(NH3predict$SoilpH>10),]

# pH categorisation
NH3predict$pH_cat <- factor(ifelse(NH3predict$SoilpH<=7, "normal", "chalky"))
ddply(NH3predict, c("pH_cat"), summarise, N=length(NH3loss))

ggplot(NH3predict, aes(x=SoilpH, fill=pH_cat)) + geom_histogram(binwidth=0.1, alpha=.5, position="identity")
ggplot(NH3predict, aes(x=SoilpH, fill=pH_cat)) + geom_density(alpha=.3)

#
# Model 4 (A4)
NH3predict$Clay <- as.numeric(as.character(NH3predict$Clay))

#
# Model 5 (A5)
NH3predict$SOC <- as.numeric(as.character(NH3predict$SOC))
# SOC categorisation 
# Organic soils: >30% organic matter
# >12% SOC
NH3predict$SOC_cat <- factor(ifelse(NH3predict$SOC>=12, "organic", "mineral"))
ggplot(NH3predict, aes(x=SOC, fill=SOC_cat)) + geom_histogram(binwidth=0.1, alpha=.5, position="identity")
ggplot(NH3predict, aes(x=SOC, fill=SOC_cat)) + geom_density(alpha=.3)

#
# Model 6
NH3predict$Crop <- factor(NH3predict$Crop)
# Cover crop categorisation
NH3predict$Cover <- NH3predict$Crop
levels(NH3predict$Cover)[levels(NH3predict$Cover)=="alfalfa"] <- "Arable" 
levels(NH3predict$Cover)[levels(NH3predict$Cover)=="cereal"] <- "Arable"
levels(NH3predict$Cover)[levels(NH3predict$Cover)=="grass"] <- "Grass"
levels(NH3predict$Cover)[levels(NH3predict$Cover)=="beans"] <- "Arable"
levels(NH3predict$Cover)[levels(NH3predict$Cover)=="maize"] <- "Arable"
levels(NH3predict$Cover)[levels(NH3predict$Cover)=="none"] <- "None"
levels(NH3predict$Cover)[levels(NH3predict$Cover)=="other"] <- "Other"
levels(NH3predict$Cover)[levels(NH3predict$Cover)=="trees"] <- "Other"
summary(NH3predict$Cover)

#
# Model 7-8
NH3predict$Rainfall <- as.numeric(as.character(NH3predict$Rainfall))
NH3predict$AirTemperature <- as.numeric(as.character(NH3predict$AirTemperature))
NH3predict$Duration <- as.numeric(as.character(NH3predict$Duration))
NH3predict$Rain.mean <- NH3predict$Rainfall/NH3predict$Duration

# Model 9
# Binarising irrigation
NH3predict$irrigation_cat <- factor(ifelse(NH3predict$irrigation=="0", 0, 1))

# Model 10
NH3predict$WaterPerDay <- as.numeric(as.character(NH3predict$WaterPerDay))

#
#
# Reordering factor levels

NH3predict$Fertiliser <- factor(NH3predict$Fertiliser, levels = c("urea+", "UAN", "ammonium+1", "ammonium+2"))
NH3predict$Method <- relevel(NH3predict$Method, "Micromet")
NH3predict$Application <- relevel(NH3predict$Application, "Broadcast") 
NH3predict$Place <- relevel(NH3predict$Place, "Outdoor") 
NH3predict$pH_cat <- relevel(NH3predict$pH_cat, "normal") 

# Model response variable
NH3predict$NH3loss <- as.numeric(as.character(NH3predict$NH3loss))


# Full model ----

# Model for calculating weights
model <- lmer(sqrt(NH3loss +4) ~ Fertiliser + Application + Method + Place + (1|Ref),
              control=lmerControl(optimizer="bobyqa"),
              data=NH3predict)

# Calculating weights
NH3predict$Weight <- 1/lm(abs(resid(model)) ~ fitted(model))$fitted.values^2

# Model A_Base 
m1_lmer_p1 <- lmer(sqrt(NH3loss+4) ~ Fertiliser + Application + Method + Place + (1|Ref),
                   weights=Weight,
                   control=lmerControl(optimizer="bobyqa"),
                   data=NH3predict)

# Combined model
NH3predict_mALL <- NH3predict[which(is.na(NH3predict$pH_cat)==FALSE & is.na(NH3predict$Clay)==FALSE & is.na(NH3predict$AirTemperature)==FALSE),]

mALL_lmer <- lmer(sqrt(NH3loss+4) ~ Fertiliser + Application + Method + Place + 
                    pH_cat + Fertiliser:pH_cat + 
                    Clay + AirTemperature +
                    (1|Ref),
                  weights=Weight,
                  control=lmerControl(optimizer="bobyqa"),
                  data=NH3predict_mALL)
summary(mALL_lmer)


#
#
# SENSITIVITY ----

# Sensitivity to outliers ----

NH3predict_urea <- dplyr::filter(NH3predict, Fertiliser=="urea+")

table(#NH3predict_urea$Fertiliser,
      NH3predict_urea$Method,
      NH3predict_urea$pH_cat)
      #NH3predict_mALL$Application,
      #NH3predict_mALL$Place)

NH3predict_urea <- dplyr::filter(NH3predict_urea, Method!="15N")
NH3predict_urea$Method <- factor(NH3predict_urea$Method)

table(#NH3predict_urea$Fertiliser,
  NH3predict_urea$Method,
  #NH3predict_urea$pH_cat,
  NH3predict_urea$Application)
  #NH3predict_mALL$Place)

NH3predict_urea <- dplyr::filter(NH3predict_urea, Application!="Injected")
NH3predict_urea$Application <- factor(NH3predict_urea$Application)

table(#NH3predict_urea$Fertiliser,
  NH3predict_urea$Method,
  #NH3predict_urea$pH_cat,
  NH3predict_urea$Application)
#NH3predict_mALL$Place)

NH3predict_urea <- dplyr::filter(NH3predict_urea, Method!="Semi-open")
NH3predict_urea$Method <- factor(NH3predict_urea$Method)

table(#NH3predict_urea$Fertiliser,
  NH3predict_urea$Method,
  #NH3predict_urea$pH_cat,
  #NH3predict_urea$Application)
  NH3predict_urea$Place)

NH3predict_urea <- dplyr::filter(NH3predict_urea, Place=="Outdoor")
NH3predict_urea$Place <- factor(NH3predict_urea$Place)


mALL_urea <- lmer(sqrt(NH3loss+4) ~ Application + Method +  
                    pH_cat +  
                    Clay + AirTemperature +
                    (1|Ref),
                  weights=Weight,
                  control=lmerControl(optimizer="bobyqa"),
                  data=NH3predict_urea)

#
#
# Compare models ----

anova(mALL_lmer)
anova(mALL_urea)

# Get effect sizes

#
# Fertiliser

# Effect sizes
at <- list("Application"="Broadcast", "Method"="Micromet", "Place"="Outdoor")
mALL_est_F <- emmeans(mALL_lmer, ~ Fertiliser, type="response", at=at)
mALL_est_F  <- data.frame(mALL_est_F)

#
# Application

# Effect sizes
at <- list("Fertiliser"="urea+", "Method"="Micromet", "Place"="Outdoor")
mALL_est_A <- emmeans(mALL_lmer, ~ Application, type="response", at=at)
mALL_est_A <- data.frame(mALL_est_A)

mALL_est_A_urea <- emmeans(mALL_urea, ~ Application, type="response", at=at)
mALL_est_A_urea  <- data.frame(mALL_est_A_urea)


#
# Method

# Effect sizes
at <- list("Fertiliser"="urea+", "Application"="Broadcast", "Place"="Outdoor")
mALL_est_M <- emmeans(mALL_lmer, ~ Method, type="response", at=at)
mALL_est_M <- data.frame(mALL_est_M)

mALL_est_M_urea <- emmeans(mALL_urea, ~ Method, type="response", at=at)
mALL_est_M_urea <- data.frame(mALL_est_M_urea)

#
# Place

# Effect size
at <- list("Fertiliser"="urea+", "Application"="Broadcast", "Method"="Micromet")
mALL_est_P <- emmeans(mALL_lmer, ~ Place, type="response", at=at)
mALL_est_P <- data.frame(mALL_est_P)

#

## 
# Collating a dataframe with all effect sizes

# Vectors for emmeans-objects and fectors
Model_vec <- c("mALL_est_A", "mALL_est_M")
Model_vec_urea <- c("mALL_est_A_urea", "mALL_est_M_urea")
Factor_vec <- c("Application", "Method")

# Creating an empty dataframe
mALL_estimates <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())
mALL_estimates_urea <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())

# Loop for extracting all estimates
for(i in c(1:2)){ # Cycling through models
  m_est <- data.frame(get(Model_vec[i]))
  mALL_estimates <- rbind(mALL_estimates, data.frame(Factor=Factor_vec[i],
                                                 Level=m_est[,1],
                                                 Effect=m_est[,2],
                                                 SE=m_est[,3],
                                                 df=m_est[,4],
                                                 lower.CL=m_est[,5],
                                                 upper.CL=m_est[,6]))
}

for(i in c(1:2)){ # Cycling through models
  m_est <- data.frame(get(Model_vec_urea[i]))
  mALL_estimates_urea <- rbind(mALL_estimates_urea, data.frame(Factor=Factor_vec[i],
                                                       Level=m_est[,1],
                                                       Effect=m_est[,2],
                                                       SE=m_est[,3],
                                                       df=m_est[,4],
                                                       lower.CL=m_est[,5],
                                                       upper.CL=m_est[,6]))
}


# Compare effect sizes ----

# Chech the effects are different
mALL_estimates$Effect == mALL_estimates_urea$Effect

# Effects in tables
mALL_estimates
mALL_estimates_urea

# In a figure ----
mALL_estimates$Model <- rep("full", length(mALL_estimates$Factor))
mALL_estimates_urea$Model <- rep("Urea subset", length(mALL_estimates_urea$Factor))
mALL_estimate_comp <- rbind(mALL_estimates, mALL_estimates_urea)

# Plot
ggplot(mALL_estimate_comp, 
       aes(y=Level, x=Effect, xmin=lower.CL, xmax=upper.CL, 
           colour=Model)) +
  theme_bw(base_size = 11) +
  facet_grid(Factor ~ ., scales = "free_y", space = "free_y") +
  geom_point(size=3, position=position_dodge(.5)) + 
  geom_errorbarh(height=.1, position=position_dodge(.5))  +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  scale_colour_manual(values=c("steelblue", "firebrick")) +
  labs(x="Effect Size (NH3 loss (%))", y="") 

# Save the figure
ggsave(here::here("output/sensitivity_plus4.png"), last_plot())  


#
#
# SENSITIVITY ----

# Sensitivity to sqrt transformation ----

# Exclude negative effects from the data
NH3_noneg <- filter(NH3predict, NH3loss>0)

# Rerun model (no negatives)
# Model for calculating weights (no negatives)
model <- lmer(sqrt(NH3loss) ~ Fertiliser + Application + Method + Place + (1|Ref),
              control=lmerControl(optimizer="bobyqa"),
              data=NH3_noneg)

# Calculating weights (no negatives)
NH3_noneg$Weight <- 1/lm(abs(resid(model)) ~ fitted(model))$fitted.values^2

# Model A_Base (no negatives)
m1_lmer_p1_nn <- lmer(sqrt(NH3loss) ~ Fertiliser + Application + Method + Place + (1|Ref),
                   weights=Weight,
                   control=lmerControl(optimizer="bobyqa"),
                   data=NH3_noneg)

#
# Model assumptions (no negatives)
par(mfrow=c(1,3))
plot(fitted(m1_lmer_p1_nn), NH3_noneg$Weight*resid(m1_lmer_p1_nn)) 
abline(h=0) # OK
qqnorm(resid(m1_lmer_p1_nn))
qqline(resid(m1_lmer_p1_nn)) 
hist(resid(m1_lmer_p1_nn)) 
par(mfrow=c(1,1))

#
#
# Compare models ----

anova(m1_lmer_p1)
anova(m1_lmer_p1_nn)

# Get effect sizes

#
# Fertiliser

# Effect sizes
at <- list("Application"="Broadcast", "Method"="Micromet", "Place"="Outdoor")
m1_est_F <- emmeans(m1_lmer_p1, ~ Fertiliser, type="response", at=at)
m1_est_F  <- data.frame(m1_est_F)

m1_est_F_nn <- emmeans(m1_lmer_p1_nn, ~ Fertiliser, type="response", at=at)
m1_est_F_nn  <- data.frame(m1_est_F_nn)

#
# Application

# Effect sizes
at <- list("Fertiliser"="urea+", "Method"="Micromet", "Place"="Outdoor")
m1_est_A <- emmeans(m1_lmer_p1, ~ Application, type="response", at=at)
m1_est_A <- data.frame(m1_est_A)

m1_est_A_nn <- emmeans(m1_lmer_p1_nn, ~ Application, type="response", at=at)
m1_est_A_nn <- data.frame(m1_est_A_nn)

#
# Method

# Effect sizes
at <- list("Fertiliser"="urea+", "Application"="Broadcast", "Place"="Outdoor")
m1_est_M <- emmeans(m1_lmer_p1, ~ Method, type="response", at=at)
m1_est_M <- data.frame(m1_est_M)

m1_est_M_nn <- emmeans(m1_lmer_p1_nn, ~ Method, type="response", at=at)
m1_est_M_nn <- data.frame(m1_est_M_nn)

#
# Place

# Effect size
at <- list("Fertiliser"="urea+", "Application"="Broadcast", "Method"="Micromet")
m1_est_P <- emmeans(m1_lmer_p1, ~ Place, type="response", at=at)
m1_est_P <- data.frame(m1_est_P)

m1_est_P_nn <- emmeans(m1_lmer_p1_nn, ~ Place, type="response", at=at)
m1_est_P_nn <- data.frame(m1_est_P_nn)

#

## 
# Collating a dataframe with all effect sizes

# Vectors for emmeans-objects and fectors
Model_vec <- c("m1_est_F", "m1_est_A", "m1_est_M", "m1_est_P")
Model_vec_nn <- c("m1_est_F_nn", "m1_est_A_nn", "m1_est_M_nn", "m1_est_P_nn")
Factor_vec <- c("Fertiliser","Application", "Method", "Place")

# Creating an empty dataframe
m1_estimates <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())
m1_estimates_nn <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())

# Loop for extracting all estimates
for(i in c(1:4)){ # Cycling through models
  m_est <- data.frame(get(Model_vec[i]))
  m1_estimates <- rbind(m1_estimates, data.frame(Factor=Factor_vec[i],
                                                 Level=m_est[,1],
                                                 Effect=m_est[,2],
                                                 SE=m_est[,3],
                                                 df=m_est[,4],
                                                 lower.CL=m_est[,5],
                                                 upper.CL=m_est[,6]))
}

for(i in c(1:4)){ # Cycling through models
  m_est <- data.frame(get(Model_vec_nn[i]))
  m1_estimates_nn <- rbind(m1_estimates_nn, data.frame(Factor=Factor_vec[i],
                                                 Level=m_est[,1],
                                                 Effect=m_est[,2],
                                                 SE=m_est[,3],
                                                 df=m_est[,4],
                                                 lower.CL=m_est[,5],
                                                 upper.CL=m_est[,6]))
}


# Compare effect sizes ----

# Chech the effects are different
m1_estimates$Effect == m1_estimates_nn$Effect
m1_estimates$Effect - m1_estimates_nn$Effect
m1_estimates$SE - m1_estimates_nn$SE

# Effects in tables
m1_estimates
m1_estimates_nn

# In a figure ----
m1_estimates$Data <- rep("full", length(m1_estimates$Factor))
m1_estimates_nn$Data <- rep("no negatives", length(m1_estimates_nn$Factor))
m1_estimate_comp <- rbind(m1_estimates, m1_estimates_nn)

# Plot
ggplot(m1_estimate_comp, 
       aes(y=Level, x=Effect, xmin=lower.CL, xmax=upper.CL, 
           colour=Data)) +
  theme_bw(base_size = 11) +
  facet_grid(Factor ~ ., scales = "free_y", space = "free_y") +
  geom_point(size=3, position=position_dodge(1)) + 
  geom_errorbarh(height=.1, position=position_dodge(1))  +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  scale_colour_manual(values=c("steelblue", "firebrick")) +
  labs(x="Effect Size (NH3 loss (%))", y="") 

# Save the figure
ggsave(here::here("output/sensitivity_plus4.png"), last_plot())  

compare_performance(m1_lmer_p1, m1_lmer_p1_nn)

#A1,2,3,4
compare_performance(mALL5_lmer, #A1
                    mALLnoTemp_lmer, #A2
                    mALL4_lmer, #A3
                    mALL_lmer) #A4

# END ----
