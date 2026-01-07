library(mgcv)
library(gratia)
library(lme4)
library(ggplot2)
library(glmmTMB)
library(DHARMa)

NH3predict_mALL5 <- readRDS("data/Wrangled/NH3predict_mALL5.rds")

#
#
#

# LO5 model (most complex, least data)
mALL5_lmer <- lmer(sqrt(NH3loss+4) ~ Fertiliser + Application + Method  +  pH_cat + 
                     Fertiliser:pH_cat + 
                     Application.rate + 
                     AirTemperature + 
                     Rain.mean + 
                     AirTemperature:Rain.mean + 
                     (1|Ref),
                   weights=Weight,
                   control=lmerControl(optimizer="bobyqa"),
                   data=NH3predict_mALL5)
par(mfrow=c(1,3))
qqnorm(resid(mALL5_lmer));qqline(resid(mALL5_lmer))
hist(resid(mALL5_lmer))
plot(fitted(mALL5_lmer),resid(mALL5_lmer),);abline(h=0)
par(mfrow=c(1,1))

summary(mALL5_lmer)


mALL5_lmer_nn <- lmer(sqrt(NH3loss) ~ Fertiliser + Application + Method  +  pH_cat + 
                        Fertiliser:pH_cat + 
                        Application.rate + 
                        AirTemperature + 
                        Rain.mean + 
                        AirTemperature:Rain.mean + 
                        (1|Ref),
                   weights=Weight,
                   control=lmerControl(optimizer="bobyqa"),
                   data=filter(NH3predict_mALL5, NH3loss>0))

mALL5_lmer_c <- lmer(sqrt(NH3loss+4) ~ Fertiliser + Application + Method  +  pH_cat + 
                        Fertiliser:pH_cat + 
                        Application.rate + 
                        AirTemperature + 
                        Rain.mean + 
                        AirTemperature:Rain.mean + 
                        (1|Ref),
                      weights=Weight,
                      control=lmerControl(optimizer="bobyqa"),
                      data=filter(NH3predict_mALL5, NH3loss>0))

mALL5_glm <- glmmTMB(NH3loss ~ Fertiliser + Application + Method  +  pH_cat + 
                       Fertiliser:pH_cat + 
                       Application.rate + 
                       AirTemperature + 
                       Rain.mean + 
                       AirTemperature:Rain.mean + 
                       (1|Ref),
                     family=Gamma(link="log"), 
                     data=filter(NH3predict_mALL5, NH3loss>0))
plot(simulateResiduals(mALL5_glm))
summary(mALL5_glm)

mALL5_gam <- gam(NH3loss ~ Fertiliser + Application + Method  +  pH_cat +
                       Fertiliser:pH_cat +
                       s(Application.rate) +
                       s(AirTemperature) +
                       s(Rain.mean) +
                       s(AirTemperature,Rain.mean) +
                       s(Ref, bs="re"),
                     family=scat(link="inverse"),
                     data=NH3predict_mALL5)
appraise(mALL5_gam)
draw(mALL5_gam, residuals=T)
summary(mALL5_gam)


#
#
#

mALL_lmer <- mALL5_lmer
mALL_lmer_nn <- mALL5_lmer_nn
mALL_lmer_c <- mALL5_lmer_c
mALL_GLM <- mALL5_glm
mALL_GAM <- mALL5_gam


#
# Fertiliser

# Effect sizes
at <- list("Application"="Broadcast", "Method"="Micromet", "Place"="Outdoor")
m1_est_F <- emmeans(mALL_lmer, ~ Fertiliser, type="response", at=at)
m1_est_F  <- data.frame(m1_est_F)

m1_est_F_nn <- emmeans(mALL_lmer_nn, ~ Fertiliser, type="response", at=at)
m1_est_F_nn  <- data.frame(m1_est_F_nn)

m1_est_F_c <- emmeans(mALL_lmer_c, ~ Fertiliser, type="response", at=at)
m1_est_F_c  <- data.frame(m1_est_F_c)

m1_est_F_gam <- emmeans(mALL_GAM, ~ Fertiliser, type="response", at=at)
m1_est_F_gam  <- data.frame(m1_est_F_gam)

m1_est_F_glm <- emmeans(mALL_GLM, ~ Fertiliser, type="response", at=at)
m1_est_F_glm  <- data.frame(m1_est_F_glm)

#
# Application

# Effect sizes
at <- list("Fertiliser"="urea+", "Method"="Micromet", "Place"="Outdoor")
m1_est_A <- emmeans(mALL_lmer, ~ Application, type="response", at=at)
m1_est_A <- data.frame(m1_est_A)

m1_est_A_nn <- emmeans(mALL_lmer_nn, ~ Application, type="response", at=at)
m1_est_A_nn <- data.frame(m1_est_A_nn)

m1_est_A_c <- emmeans(mALL_lmer_c, ~ Application, type="response", at=at)
m1_est_A_c <- data.frame(m1_est_A_c)

m1_est_A_gam <- emmeans(mALL_GAM, ~ Application, type="response", at=at)
m1_est_A_gam <- data.frame(m1_est_A_gam)

m1_est_A_glm <- emmeans(mALL_GLM, ~ Application, type="response", at=at)
m1_est_A_glm  <- data.frame(m1_est_A_glm)


#
# Method

# Effect sizes
at <- list("Fertiliser"="urea+", "Application"="Broadcast", "Place"="Outdoor")
m1_est_M <- emmeans(mALL_lmer, ~ Method, type="response", at=at)
m1_est_M <- data.frame(m1_est_M)

m1_est_M_nn <- emmeans(mALL_lmer_nn, ~ Method, type="response", at=at)
m1_est_M_nn <- data.frame(m1_est_M_nn)

m1_est_M_c <- emmeans(mALL_lmer_c, ~ Method, type="response", at=at)
m1_est_M_c <- data.frame(m1_est_M_c)

m1_est_M_gam <- emmeans(mALL_GAM, ~ Method, type="response", at=at)
m1_est_M_gam <- data.frame(m1_est_M_gam)

m1_est_M_glm <- emmeans(mALL_GLM, ~ Method, type="response", at=at)
m1_est_M_glm  <- data.frame(m1_est_M_glm)

#

## 
# Collating a dataframe with all effect sizes

# Vectors for emmeans-objects and fectors
Model_vec <- c("m1_est_F", "m1_est_A", "m1_est_M")
Model_vec_nn <- c("m1_est_F_nn", "m1_est_A_nn", "m1_est_M_nn")
Model_vec_c <- c("m1_est_F_c", "m1_est_A_c", "m1_est_M_c")
Model_vec_gam <- c("m1_est_F_gam", "m1_est_A_gam", "m1_est_M_gam")
Model_vec_glm <- c("m1_est_F_glm", "m1_est_A_glm", "m1_est_M_glm")
Factor_vec <- c("Fertiliser","Application", "Method")

# Creating an empty dataframe
m1_estimates <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())
m1_estimates_nn <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())
m1_estimates_c <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())
m1_estimates_gam <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())
m1_estimates_glm <- data.frame(Factor=factor(), Level=factor(), Effect=numeric(), SE=numeric(), df=numeric(), lower.CL=numeric(), upper.CL=numeric())

# Loop for extracting all estimates
for(i in c(1:3)){ # Cycling through models
  m_est <- data.frame(get(Model_vec[i]))
  m1_estimates <- rbind(m1_estimates, data.frame(Factor=Factor_vec[i],
                                                 Level=m_est[,1],
                                                 Effect=m_est[,2],
                                                 SE=m_est[,3],
                                                 df=m_est[,4],
                                                 lower.CL=m_est[,5],
                                                 upper.CL=m_est[,6]))
}

for(i in c(1:3)){ # Cycling through models
  m_est <- data.frame(get(Model_vec_nn[i]))
  m1_estimates_nn <- rbind(m1_estimates_nn, data.frame(Factor=Factor_vec[i],
                                                       Level=m_est[,1],
                                                       Effect=m_est[,2],
                                                       SE=m_est[,3],
                                                       df=m_est[,4],
                                                       lower.CL=m_est[,5],
                                                       upper.CL=m_est[,6]))
}

for(i in c(1:3)){ # Cycling through models
  m_est <- data.frame(get(Model_vec_c[i]))
  m1_estimates_c <- rbind(m1_estimates_c, data.frame(Factor=Factor_vec[i],
                                                       Level=m_est[,1],
                                                       Effect=m_est[,2],
                                                       SE=m_est[,3],
                                                       df=m_est[,4],
                                                       lower.CL=m_est[,5],
                                                       upper.CL=m_est[,6]))
}

for(i in c(1:3)){ # Cycling through models
  m_est <- data.frame(get(Model_vec_gam[i]))
  m1_estimates_gam <- rbind(m1_estimates_gam, data.frame(Factor=Factor_vec[i],
                                                       Level=m_est[,1],
                                                       Effect=m_est[,2],
                                                       SE=m_est[,3],
                                                       df=m_est[,4],
                                                       lower.CL=m_est[,5],
                                                       upper.CL=m_est[,6]))
}

for(i in c(1:4)){ # Cycling through models
  m_est <- data.frame(get(Model_vec_glm[i]))
  m1_estimates_glm <- rbind(m1_estimates_glm, data.frame(Factor=Factor_vec[i],
                                                        Level=m_est[,1],
                                                        Effect=m_est[,2],
                                                        SE=m_est[,3],
                                                        df=m_est[,4],
                                                        lower.CL=m_est[,5],
                                                        upper.CL=m_est[,6]))
}


# Compare effect sizes ----

# Effects in tables
m1_estimates
m1_estimates_nn
m1_estimates_c
m1_estimates_gam
m1_estimates_glm

# In a figure ----
m1_estimates$Model <- rep("LMM (all data)", length(m1_estimates$Factor))
m1_estimates_nn$Model <- rep("LMM (subset, no +4)", length(m1_estimates_nn$Factor))
m1_estimates_c$Model <- rep("LMM (subset with +4)", length(m1_estimates_c$Factor))
m1_estimates_gam$Model <- rep("GAM (all data)", length(m1_estimates_gam$Factor))
m1_estimates_glm$Model <- rep("GLM (subset)", length(m1_estimates_glm$Factor))
m1_estimate_comp <- rbind(m1_estimates, m1_estimates_nn, m1_estimates_c, m1_estimates_gam, m1_estimates_glm)

# Plot
ggplot(m1_estimate_comp, 
       aes(y=Level, x=Effect, xmin=lower.CL, xmax=upper.CL, 
           colour=Model)) +
  theme_bw(base_size = 11) +
  facet_grid(Factor ~ ., scales = "free_y", space = "free_y") +
  geom_point(size=3, position=position_dodge(.5)) + 
  geom_errorbarh(height=.1, position=position_dodge(.5))  +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  scale_colour_manual(values=c("steelblue", "firebrick", "darkgreen", "black", "orange")) +
  labs(x="Effect Size (NH3 loss (%))", y="") 

# Save the figure
ggsave(here::here("output/sensitivity_plus4.png"), last_plot())  

# Are GAMs, allowing non-linear associations, better? 
library(performance)
compare_performance(mALL_lmer,
                    #mALL_lmer_nn,
                    #mALL_lmer_c,
                    #mALL_GLM,
                    mALL_GAM,
                    rank = TRUE)

summary(NH3predict$Method)
NH3predict[which(NH3predict$Method=="15N"),]

# Performance values 
# LO1, LO2, LO3, LO4, LO5
compare_performance(mALL1_lmer,#LO1
                    mALL2_lmer,#LO2
                    mALL3_lmer,#LO3
                    mALL4_lmer,#LO4
                    mALL5_lmer)#LO5


levels(mALL1_estimates$Level) <- c("Urea+",
                                   "UAN",
                                   "Ammonium+1",
                                   "Ammonium+2",
                                   "Broadcast",
                                   "Incorporated",
                                   "Injected",
                                   "Liquid",
                                   "Micromet",
                                   "15N",
                                   "Closed chamber",
                                   "Drager-Tube",
                                   "Semi-open",
                                   "Ventilated chambers",
                                   "Windtunnel",
                                   "Normal",
                                   "Chalky",
                                   "Urea+ : Normal",
                                   "Urea+ : Chalky",
                                   "UAN : Normal",
                                   "UAN : Chalky",  
                                   "Ammonium+1 : Normal",
                                   "Ammonium+1 : Chalky",
                                   "Ammonium+2 : Normal",
                                   "Ammonium+2 : Chalky")
levels(mALL5_estimates$Level) <- c("Urea+",
                                   "UAN",
                                   "Ammonium+1",
                                   "Ammonium+2",
                                   "Broadcast",
                                   "Incorporated",
                                   "Injected",
                                   "Liquid",
                                   "Micromet",
                                   "15N",
                                   "Closed chamber",
                                   "Drager-Tube",
                                   "Semi-open",
                                   "Ventilated chambers",
                                   "Windtunnel",
                                   "Outdoor",
                                   "Laboratory",
                                   "Normal",
                                   "Chalky",
                                   "Urea+ : Normal",
                                   "Urea+ : Chalky",
                                   "UAN : Normal",
                                   "UAN : Chalky",  
                                   "Ammonium+1 : Normal",
                                   "Ammonium+1 : Chalky",
                                   "Ammonium+2 : Normal",
                                   "Ammonium+2 : Chalky")

# Full figure 
LO1 <-
ggplot(mALL1_estimates, aes(y=Level, x=Effect, xmin=lower.CL, xmax=upper.CL, colour=Factor)) +
  theme_bw(base_size = 11) +
  facet_grid(Factor ~ ., scales = "free_y", space = "free_y") +
  geom_point(size=3) + 
  geom_errorbarh(height=.1)  +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  scale_colour_manual(values=c("darkgreen", "firebrick3", "purple4", "mediumblue", "palegreen4")) +
  #scale_color_brewer(palette = "Dark2") +
  labs(x="Effect Size (NH3 loss (%))", y="") +
  theme(legend.position = "none") +
  xlim(-0.5, 35)+
  ggtitle("b)")

# Full figure 
A1 <-
ggplot(mALL5_estimates, aes(y=Level, x=Effect, xmin=lower.CL, xmax=upper.CL, colour=Factor)) +
  theme_bw(base_size = 11) +
  facet_grid(Factor ~ ., scales = "free_y", space = "free_y") +
  geom_point(size=3) + 
  geom_errorbarh(height=.1)  +
  geom_vline(xintercept=0, colour="black", linetype="dashed") +
  scale_colour_manual(values=c("darkgreen", "firebrick3", "purple4", "mediumblue","orange", "palegreen4")) +
  #scale_color_brewer(palette = "Dark2") +
  labs(x="Effect Size (NH3 loss (%))", y="") +
  theme(legend.position = "none")+
  xlim(-0.5, 35)+
  ggtitle("a)")

ggsave("Model_LO1.png", LO1, width=15, height=15, units="cm")
ggsave("Model_A1.png", A1, width=15, height=15, units="cm")

library(patchwork)
combined_plot <- A1 + LO1 + plot_layout(ncol = 2, 
                                        guides = "collect")

combined_plot
ggsave("Model_A1_LO1.png", combined_plot, width=30, height=15, units="cm")
