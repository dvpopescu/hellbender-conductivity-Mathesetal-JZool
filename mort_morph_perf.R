library(survival)
library(lme4)
library(nlme)
library(mclogit)
library(dplyr)
library(FSA)
library(PairedData)
library(emmeans)
library(ggplot2)
library(tidyverse)
library(plotrix)
library(ggfortify)

# This code contains 3 sections for three separate analyses to assess the 
# effects of conductivity treatments on eastern hellbender larval 
#    1) mortality (using counts of surviving animals - Poisson models)
#    2) morphology (PCA analysis, mass, Snout-Vent Length, total length,
#        tail height, head width, mid width)
#    3) performance (number of bursts, speed, total distance)


# MORTALITY - COUNT DATA --------------------------------------------------


mort <- read.csv("mort_Oct23.csv", header = T)
summary(mort)

mort$tank <- as.factor(mort$tank)
mort$clutch <- as.factor(mort$clutch)
mort$rack <- as.factor(mort$rack)
mort$ba <- as.factor(mort$ba)
mort$cond <- as.factor(mort$cond)
mort$impact <- as.factor(mort$impact)
mort$ud <- as.factor(mort$ud)
mort$switch <- as.factor(mort$switch)
mort$alternate <- as.factor(mort$alternate)

# examine data
boxplot(mort$pre_prop~mort$cond * mort$ba)
boxplot(mort$pre_prop~mort$alternate*mort$switch)

group_by(mort, c(alternate)) %>%
  summarise(
    count = n(),
    median = median(pre, na.rm = TRUE),
    IQR = IQR(pre, na.rm = TRUE)
  )

mort.a <- subset(mort, ba == "after")
mort.a
mort.a$alternate <- factor(mort.a$cond,
                           levels = c('low','high'), ordered = TRUE)
mort.a
boxplot(mort.a$pre_prop~mort.a$switch*mort.a$cond)

# BA in unswitched treatments AND switched treatments before (or all Controls), 
# to look at interaction between BA and Cond in Control state treatments
mort.c <- subset(mort, impact == "c")
summary(mort.c)
mort.c$alternate <- factor(mort.c$alternate,
                             levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
m.c1 <- glmer(pre ~ ba * cond + (1|clutch), family = 'poisson', data=mort.c) #lowest AIC
m.c2 <- glmer(pre ~ ba * cond + (1|tank), family = 'poisson', data=mort.c)
m.c3 <- glmer(pre ~ ba * cond + (1|tank)+ (1|clutch), family = 'poisson', data=mort.c)
AIC(m.c1, m.c2, m.c3)

summary(m.c1)
car::Anova(m.c1, type = "II") #time*cond and cond S
emmeans(m.c1, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
                                    # (before low and after low) 
                                    # (before high and after high) * <
                                    # (before low and before high) <
                                    # (after low and after high) * <

boxplot(pre_prop ~ alternate, data=mort.c, ylim = c(0,1), xlab ="Treatment Condition", ylab = "Survival by Tank", names = c("before low", "after low", "before high", "after high"))
stripchart(pre_prop ~ alternate, data = mort.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    


### compare switched to their non-switched counterparts

hi.to.lo <- subset(mort, (switch == "n" & alternate == "afterhigh") | 
                     (switch == "y" & alternate =="afterlow") | alternate == "beforehigh")
hi.to.lo$alternate <- factor(hi.to.lo$alternate,
                             levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(pre_prop ~ alternate, data = hi.to.lo) # no catching up when exposed to good conditions
boxplot(pre_prop*100 ~ alternate, data = hi.to.lo, ylim = c(0,100), xlab ="Treatment Condition", ylab = " Survival by Tank (%)", names = c("before high", "after high", "after low"))
stripchart(pre_prop*100 ~ alternate, data = hi.to.lo,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)       

m.hi.to.lo1 <- glmer(pre ~ alternate + (1|clutch), family='poisson', data=hi.to.lo) #lowest AIC
m.hi.to.lo2 <- glmer(pre ~ alternate + (1|tank), family='poisson', data=hi.to.lo)
m.hi.to.lo3 <- glmer(pre ~ alternate + (1|clutch)+ (1|tank), family='poisson', data=hi.to.lo)
AIC(m.hi.to.lo1, m.hi.to.lo2, m.hi.to.lo3)

summary(m.hi.to.lo1)
car::Anova(m.hi.to.lo1, type = "II") #swtich S
emmeans(m.hi.to.lo1, pairwise ~ alternate)

lo.to.hi <- subset(mort, switch == "n" & alternate == "afterlow" | 
                     switch == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi$alternate <- factor(lo.to.hi$alternate,
                             levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(pre ~ alternate, data = lo.to.hi) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(pre_prop*100 ~ alternate, data = lo.to.hi, ylim = c(0,100), xlab ="Treatment Condition", ylab = " Survival by Tank (%)", names = c("before low", "after low", "after high"))
stripchart(pre_prop*100 ~ alternate, data = lo.to.hi,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)       

m.lo.to.hi1 <- glmer(pre ~ alternate + (1|clutch), family='poisson', data=lo.to.hi) #lowest AIC
m.lo.to.hi2 <- glmer(pre ~ alternate + (1|tank), family='poisson', data=lo.to.hi)
m.lo.to.hi3 <- glmer(pre ~ alternate + (1|clutch)+ (1|tank), family='poisson', data=lo.to.hi)
AIC(m.lo.to.hi1, m.lo.to.hi2, m.lo.to.hi3)

summary(m.lo.to.hi1)
car::Anova(m.lo.to.hi1, type = "II") #switch NS
emmeans(m.lo.to.hi1, pairwise ~ alternate)

### END MORTALITY ####


# MORPHOLOGY --------------------------------------------------------------



morph <- read.csv("morph_Oct23.csv", header = T)
summary(morph)

morph$tank <- as.factor(morph$tank)
morph$clutch <- as.factor(morph$clutch)
morph$rack <- as.factor(morph$rack)
morph$ba <- as.factor(morph$ba)
morph$cond <- as.factor(morph$cond)
morph$impact <- as.factor(morph$impact)
morph$ud <- as.factor(morph$ud)
morph$switched <- as.factor(morph$switched)
morph$taillength <- morph$Total.Length-morph$SVL

hist(morph$body.mass)
hist(morph$Total.Length)
hist(morph$SVL)
hist(morph$Head.Width)
hist(morph$Mid.Width)
hist(morph$Tail.Width)
hist(morph$Tail.Height)

### Test for correlation between morph metrics ###
df<-cor(morph[,c(2:8,17)], use="complete.obs")
write.csv(df,file="morph.cor.csv")

### PCA ###
morphPCA <- prcomp(na.omit(morph[,c(2:8,17)]), center = TRUE, scale. = TRUE)
summary(morphPCA)
str(morphPCA)
autoplot(morphPCA, x = 1, y = 2, loadings.label = TRUE) 
autoplot(morphPCA, x = 2, y = 3, loadings.label = TRUE) # PC2 separates SVL and total length 
autoplot(morphPCA, x = 2, y = 4, loadings.label = TRUE)
autoplot(morphPCA, x = 2, y = 5, loadings.label = TRUE) # PC5 separates tail metrics 
autoplot(morphPCA, x = 2, y = 7, loadings.label = TRUE) # PC7 separates total length from SVL
morphPCA$rotation

morphPC<-morphPCA$x
morph$PC1 <- morphPC[1]
morph$PC2 <- morphPC[2]
summary(morph)

cor(morphPCA$x, na.omit(morph[,c(2:8,17)],))

### PC2 models ###

morph.c <- subset(morph, impact == "c")
PC.c <- lmer(PC2 ~ ba * cond + (1|tank), data=morph.c) #lowest AIC

summary(PC.c)
car::Anova(PC.c, type = "II")
emmeans(PC.c, pairwise ~ ba*cond)

### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(PC2 ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions

PC.hi.to.lo1 <- lmer(PC2 ~ alternate + (1|tank), data=hi.to.lo1) # lowest AIC
PC.hi.to.lo2 <- lmer(PC2 ~ alternate + (1|clutch), data=hi.to.lo1)
PC.hi.to.lo3 <- lmer(PC2 ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(PC.hi.to.lo1, PC.hi.to.lo2, PC.hi.to.lo3)

summary(PC.hi.to.lo1)
car::Anova(PC.hi.to.lo1, type = "II") #switch NS
emmeans(PC.hi.to.lo1, pairwise ~ alternate)


lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                      switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(PC2 ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)

PC.lo.to.hi1 <- lmer(PC2 ~ alternate + (1|tank), data=lo.to.hi1)
PC.lo.to.hi2 <- lmer(PC2 ~ alternate + (1|clutch), data=lo.to.hi1)
PC.lo.to.hi3 <- lmer(PC2 ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1) #lowest AIC
AIC(PC.lo.to.hi1, PC.lo.to.hi2, PC.lo.to.hi3)

summary(PC.lo.to.hi3)
car::Anova(PC.lo.to.hi3, type = "II")
emmeans(PC.lo.to.hi3, pairwise ~ alternate)

### body mass models ###
morph.c <- subset(morph, impact == "c")
mass1.c <- lmer(body.mass ~ ba * cond + (1|tank), data=morph.c) #lowest AIC
mass2.c <- lmer(body.mass ~ ba * cond + (1|clutch), data=morph.c)
mass3.c <- lmer(body.mass ~ ba * cond + (1|tank) + (1|clutch), data=morph.c)
AIC(mass1.c, mass2.c, mass3.c)

summary(mass1.c)
car::Anova(mass1.c, type = "II") #cond S
emmeans(mass1.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

morph.c$alternate <- factor(morph.c$alternate,
                           levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(body.mass ~ ba*cond, data=morph.c)
boxplot(body.mass ~ alternate, data=morph.c, xlab ="Treatment Condition", ylim = c(0.2,1), ylab = " Body Mass (g)", names = c("before low", "after low", "before high", "after high"))
stripchart(body.mass ~ alternate, data = morph.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   

### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                     switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                             levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(body.mass ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(body.mass ~ alternate, data=hi.to.lo1, ylim=c(0.2,1), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Body Mass (g)")
stripchart(body.mass ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)       
axis.break(axis=2,breakpos=0.3,pos=NULL,bgcol="white",breakcol="black", style="slash",brw=0.05)
mass1.hi.to.lo1 <- lmer(body.mass ~ alternate + (1|tank), data=hi.to.lo1) # lowest AIC
mass2.hi.to.lo1 <- lmer(body.mass ~ alternate + (1|clutch), data=hi.to.lo1)
mass3.hi.to.lo1 <- lmer(body.mass ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(mass1.hi.to.lo1, mass2.hi.to.lo1, mass3.hi.to.lo1)

summary(mass1.hi.to.lo1)
car::Anova(mass1.hi.to.lo1, type = "II") #switch NS
emmeans(mass1.hi.to.lo1, pairwise ~ alternate)


lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                     switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(body.mass ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(body.mass ~ alternate, data = lo.to.hi1, ylim=c(0.2,1), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Body Mass (g)")
stripchart(body.mass ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

mass1.lo.to.hi1 <- lmer(body.mass ~ alternate + (1|tank), data=lo.to.hi1) #lowest AIC
mass2.lo.to.hi1 <- lmer(body.mass ~ alternate + (1|clutch), data=lo.to.hi1)
mass3.lo.to.hi1 <- lmer(body.mass ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(mass1.lo.to.hi1, mass2.lo.to.hi1, mass3.lo.to.hi1)

summary(mass1.lo.to.hi1)
car::Anova(mass1.lo.to.hi1, type = "II") #switch S
emmeans(mass1.lo.to.hi1, pairwise ~ alternate)

### total length models ###
morph <- read.csv("morph_Oct23.csv", header = T)

morph.c <- subset(morph, impact == "c")
length1.c <- lmer(Total.Length ~ ba * cond + (1|tank), data=morph.c)
summary(length1.c)
car::Anova(length1.c, type = "II")
emmeans(length1.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(Total.Length ~ ba*cond, data=morph.c)
morph.c$alternate <- factor(morph.c$alternate,
                            levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(Total.Length*2.54 ~ ba*cond, data=morph.c)
boxplot(Total.Length*2.54 ~ alternate, data=morph.c, ylim = c(3.5,5.5), xlab ="Treatment Condition", ylab = " Total Length (cm)", names = c("before low", "after low", "before high", "after high"))
stripchart(Total.Length*2.54 ~ alternate, data = morph.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   

### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                     switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(Total.Length*2.54 ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(Total.Length*2.54 ~ alternate, data=hi.to.lo1, ylim = c(3.5,5.5), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Total Length (cm)")
stripchart(Total.Length*2.54 ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)      

length.hi.to.lo1 <- lmer(Total.Length*2.54 ~ alternate + (1|tank), data=hi.to.lo1)
summary(length.hi.to.lo1)
emmeans(length.hi.to.lo1, pairwise ~ alternate)


lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                     switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(Total.Length*2.54 ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(Total.Length*2.54 ~ alternate, data = lo.to.hi1, ylim = c(3.5,5.5), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Total Length (cm)")
stripchart(Total.Length*2.54 ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

length.lo.to.hi1 <- lmer(Total.Length*2.54 ~ alternate + (1|tank), data=lo.to.hi1)
summary(length.lo.to.hi1)
emmeans(length.lo.to.hi1, pairwise ~ alternate)

### SVL models ###
morph <- read.csv("morph_Oct23.csv", header = T)

morph.c <- subset(morph, impact == "c")
SVL1.c <- lmer(SVL*2.54 ~ ba * cond + (1|tank), data=morph.c)
SVL2.c <- lmer(SVL*2.54 ~ ba * cond + (1|clutch), data=morph.c) #lowest AIC
SVL3.c <- lmer(SVL*2.54 ~ ba * cond + (1|tank) + (1|clutch), data=morph.c)
AIC(SVL1.c, SVL2.c, SVL3.c)

summary(SVL2.c)
car::Anova(SVL2.c, type = "II") #ba (time) and cond S
emmeans(SVL2.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(SVL ~ ba*cond, data=morph.c)
morph.c$alternate <- factor(morph.c$alternate,
                            levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(SVL*2.54 ~ ba*cond, data=morph.c)
boxplot(SVL*2.54 ~ alternate, data=morph.c, ylim = c(2,3.5), xlab ="Treatment Condition", ylab = " SVL (cm)", names = c("before low", "after low", "before high", "after high"))
stripchart(SVL*2.54 ~ alternate, data = morph.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   

### compare switched to their non-switched conterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(SVL ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(SVL*2.54 ~ alternate, data=hi.to.lo1, ylim = c(2,3.5), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="SVL (cm)")
stripchart(SVL*2.54 ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)      

SVL.hi.to.lo1 <- lmer(SVL*2.54 ~ alternate + (1|tank), data=hi.to.lo1)
SVL.hi.to.lo2 <- lmer(SVL*2.54 ~ alternate + (1|clutch), data=hi.to.lo1) #lowest AIC
SVL.hi.to.lo3 <- lmer(SVL*2.54 ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1) 
AIC(SVL.hi.to.lo1, SVL.hi.to.lo2, SVL.hi.to.lo3)

summary(SVL.hi.to.lo2)
car::Anova(SVL.hi.to.lo2, type = "II") #switch NS
emmeans(SVL.hi.to.lo2, pairwise ~ alternate)

lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(SVL ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(SVL*2.54 ~ alternate, data = lo.to.hi1, ylim = c(2,3.5), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="SVL (cm)")
stripchart(SVL*2.54 ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

SVL.lo.to.hi1 <- lmer(SVL*2.54 ~ alternate + (1|tank), data=lo.to.hi1)
SVL.lo.to.hi2 <- lmer(SVL*2.54 ~ alternate + (1|clutch), data=lo.to.hi1) #lowest AIC
SVL.lo.to.hi3 <- lmer(SVL*2.54 ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(SVL.lo.to.hi1, SVL.lo.to.hi2, SVL.lo.to.hi3)

summary(SVL.lo.to.hi2)
car::Anova(SVL.lo.to.hi2, type = "II") #switch S
emmeans(SVL.lo.to.hi2, pairwise ~ alternate)

### tail length models ###
morph <- read.csv("morph_Oct23.csv", header = T)
morph$taillength <- (morph$Total.Length-morph$SVL)*2.54

morph.c <- subset(morph, impact == "c")
taillength.c1 <- lmer(taillength  ~ ba * cond + (1|tank), data=morph.c)
taillength.c2 <- lmer(taillength  ~ ba * cond + (1|clutch), data=morph.c) #lowest AIC
taillength.c3 <- lmer(taillength  ~ ba * cond + (1|tank) + (1|clutch), data=morph.c)
AIC(taillength.c1, taillength.c2, taillength.c3)

summary(taillength.c2)
car::Anova(taillength.c2, type = "II") #cond S
emmeans(taillength.c2, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(taillength ~ ba*cond, data=morph.c)
morph.c$alternate <- factor(morph.c$alternate,
                            levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(taillength ~ ba*cond, data=morph.c)
boxplot(taillength ~ alternate, data=morph.c, ylim = c(1,3), xlab ="Treatment Condition", ylab = " Tail Length (cm)", names = c("before low", "after low", "before high", "after high"))
stripchart(taillength ~ alternate, data = morph.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   

### compare switched to their non-switched conterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$taillength <- (hi.to.lo1$Total.Length-hi.to.lo1$SVL)*2.54
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(taillength ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(taillength ~ alternate, data=hi.to.lo1, ylim = c(1,3), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Tail Length (cm)")
stripchart(taillength ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)      

taillength.hi.to.lo1 <- lmer(taillength ~ alternate + (1|tank), data=hi.to.lo1) #lowest AIC
taillength.hi.to.lo2 <- lmer(taillength ~ alternate + (1|clutch), data=hi.to.lo1)
taillength.hi.to.lo3 <- lmer(taillength ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(taillength.hi.to.lo1, taillength.hi.to.lo2, taillength.hi.to.lo3)

summary(taillength.hi.to.lo1)
car::Anova(taillength.hi.to.lo1, type = "II")
emmeans(taillength.hi.to.lo1, pairwise ~ alternate)


lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                      switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$taillength <- (lo.to.hi1$Total.Length-lo.to.hi1$SVL)*2.54
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(taillength ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(taillength ~ alternate, data = lo.to.hi1, ylim = c(1,3), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Tail Length (cm)")
stripchart(taillength ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

taillength.lo.to.hi1 <- lmer(taillength ~ alternate + (1|tank), data=lo.to.hi1) #lowest AIC
taillength.lo.to.hi2 <- lmer(taillength ~ alternate + (1|clutch), data=lo.to.hi1)
taillength.lo.to.hi3 <- lmer(taillength ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(taillength.lo.to.hi1, taillength.lo.to.hi2, taillength.lo.to.hi3)

summary(taillength.lo.to.hi1)
car::Anova(taillength.lo.to.hi1, type = "II") #switch S
emmeans(taillength.lo.to.hi1, pairwise ~ alternate)

### Head Width models ###
morph <- read.csv("morph_Oct23.csv", header = T)

morph.c <- subset(morph, impact == "c")
Head1.c <- lmer(Head.Width*2.54 ~ ba * cond + (1|tank), data=morph.c) #lowest AIC
Head2.c <- lmer(Head.Width*2.54 ~ ba * cond + (1|clutch), data=morph.c)
Head3.c <- lmer(Head.Width*2.54 ~ ba * cond + (1|tank) + (1|clutch), data=morph.c)
AIC(Head1.c, Head2.c, Head3.c)

summary(Head1.c)
car::Anova(Head1.c, type = "II") #cond S
emmeans(Head1.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(Head.Width*2.54 ~ ba*cond, data=morph.c)
morph.c$alternate <- factor(morph.c$alternate,
                            levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(Head.Width*2.54 ~ ba*cond, data=morph.c)
boxplot(Head.Width*2.54 ~ alternate, data=morph.c, ylim = c(0.4,1), xlab ="Treatment Condition", ylab = " Head Width (cm)", names = c("before low", "after low", "before high", "after high"))
stripchart(Head.Width*2.54 ~ alternate, data = morph.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   

### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(Head.Width ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(Head.Width*2.54 ~ alternate, data=hi.to.lo1, ylim = c(0.4,1), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Head Width (cm)")
stripchart(Head.Width*2.54 ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE) 


Head.hi.to.lo1 <- lmer(Head.Width*2.54 ~ alternate + (1|tank), data=hi.to.lo1) #lowest AIC
Head.hi.to.lo2 <- lmer(Head.Width*2.54 ~ alternate + (1|clutch), data=hi.to.lo1)
Head.hi.to.lo3 <- lmer(Head.Width*2.54 ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(Head.hi.to.lo1, Head.hi.to.lo2, Head.hi.to.lo3)

summary(Head.hi.to.lo1)
car::Anova(Head.hi.to.lo1, type = "II") #switch S
emmeans(Head.hi.to.lo1, pairwise ~ alternate)


l.lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(Head.Width ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(Head.Width*2.54 ~ alternate, data = lo.to.hi1, ylim = c(0.4,1), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Head Width (cm)")
stripchart(Head.Width*2.54 ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

Head.lo.to.hi1 <- lmer(Head.Width*2.54 ~ alternate + (1|tank), data=lo.to.hi1)
Head.lo.to.hi2 <- lmer(Head.Width*2.54 ~ alternate + (1|clutch), data=lo.to.hi1) #lowest AIC
Head.lo.to.hi3 <- lmer(Head.Width*2.54 ~ alternate + (1|tank) +  (1|clutch), data=lo.to.hi1)
AIC(Head.lo.to.hi1, Head.lo.to.hi2, Head.lo.to.hi3)

summary(Head.lo.to.hi2)
car::Anova(Head.lo.to.hi2, type = "II") #switch S
emmeans(Head.lo.to.hi2, pairwise ~ alternate)

### Mid Width models ###
morph <- read.csv("morph_Oct23.csv", header = T)

morph.c <- subset(morph, impact == "c")
mid1.c <- lmer(Mid.Width*2.54 ~ ba * cond + (1|tank), data=morph.c) #lowest AIC
mid2.c <- lmer(Mid.Width*2.54 ~ ba * cond + (1|clutch), data=morph.c)
mid3.c <- lmer(Mid.Width*2.54 ~ ba * cond + (1|tank) + (1|clutch), data=morph.c)
AIC(mid1.c, mid2.c, mid3.c)

summary(mid1.c)
car::Anova(mid1.c, type = "II")
emmeans(mid1.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(Mid.Width ~ ba*cond, data=morph.c)
boxplot(Mid.Width*2.54 ~ ba*cond, data=morph.c)
morph.c$alternate <- factor(morph.c$alternate,
                            levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(Mid.Width*2.54 ~ ba*cond, data=morph.c)
boxplot(Mid.Width*2.54 ~ alternate, data=morph.c, ylim = c(0.3,0.75), xlab ="Treatment Condition", ylab = " Mid Width (cm)", names = c("before low", "after low", "before high", "after high"))
stripchart(Mid.Width*2.54 ~ alternate, data = morph.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   

### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(Mid.Width ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(Mid.Width*2.54 ~ alternate, data=hi.to.lo1, ylim = c(0.3,0.75), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Mid Width (cm)")
stripchart(Mid.Width*2.54 ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE) 
mid.hi.to.lo1 <- lmer(Mid.Width*2.54 ~ alternate + (1|tank), data=hi.to.lo1) #lowest AIC
mid.hi.to.lo2 <- lmer(Mid.Width*2.54 ~ alternate + (1|clutch), data=hi.to.lo1)
mid.hi.to.lo3 <- lmer(Mid.Width*2.54 ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC (mid.hi.to.lo1, mid.hi.to.lo2, mid.hi.to.lo3)

summary(mid.hi.to.lo1)
emmeans(mid.hi.to.lo1, pairwise ~ alternate)


l.lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(Mid.Width ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(Mid.Width*2.54 ~ alternate, data = lo.to.hi1, ylim = c(0.3,0.75), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Mid Width (cm)")
stripchart(Mid.Width*2.54 ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

mid.lo.to.hi1 <- lmer(Mid.Width*2.54 ~ alternate + (1|tank), data=lo.to.hi1) #lowest AIC
mid.lo.to.hi2 <- lmer(Mid.Width*2.54 ~ alternate + (1|clutch), data=lo.to.hi1)
mid.lo.to.hi3 <- lmer(Mid.Width*2.54 ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(mid.lo.to.hi1, mid.lo.to.hi2, mid.lo.to.hi3)

summary(mid.lo.to.hi1)
emmeans(mid.lo.to.hi1, pairwise ~ alternate)

### Tail Width models ###
morph <- read.csv("morph_Oct23.csv", header = T)

morph.c <- subset(morph, impact == "c")
tail1.c <- lmer(Tail.Width*2.54 ~ ba * cond + (1|tank), data=morph.c) #lowest AIC
tail2.c <- lmer(Tail.Width*2.54 ~ ba * cond + (1|clutch), data=morph.c)
tail3.c <- lmer(Tail.Width*2.54 ~ ba * cond + (1|tank) + (1|clutch), data=morph.c)
AIC(tail1.c, tail2.c, tail3.c)

summary(tail1.c)
car::Anova(tail1.c, type = "II") #cond*switch and cond S
emmeans(tail1.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(Tail.Width ~ ba*cond, data=morph.c)
morph.c$alternate <- factor(morph.c$alternate,
                            levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(Tail.Width*2.54 ~ ba*cond, data=morph.c)
boxplot(Tail.Width*2.54 ~ alternate, data=morph.c, ylim = c(0.15,0.55), xlab ="Treatment Condition", ylab = " Base of Tail Width (cm)", names = c("before low", "after low", "before high", "after high"))
stripchart(Tail.Width*2.54 ~ alternate, data = morph.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   

### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(Tail.Width ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(Tail.Width*2.54 ~ alternate, data=hi.to.lo1, ylim = c(0.15,0.55), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Base of Tail Width (cm)")
stripchart(Tail.Width*2.54 ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE) 

tail.hi.to.lo1 <- lmer(Tail.Width*2.54 ~ alternate + (1|tank), data=hi.to.lo1) 
tail.hi.to.lo2 <- lmer(Tail.Width*2.54 ~ alternate + (1|clutch), data=hi.to.lo1)
tail.hi.to.lo3 <- lmer(Tail.Width*2.54 ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(tail.hi.to.lo1, tail.hi.to.lo2, tail.hi.to.lo3) #all the same

summary(tail.hi.to.lo1)
car::Anova(tail.hi.to.lo1, type = "II") #switch S
emmeans(tail.hi.to.lo1, pairwise ~ alternate)


l.lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(Tail.Width ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(Tail.Width*2.54 ~ alternate, data = lo.to.hi1, ylim = c(0.15,0.55), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Base of Tail Width (cm)")
stripchart(Tail.Width*2.54 ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

tail.lo.to.hi1 <- lmer(Tail.Width*2.54 ~ alternate + (1|tank), data=lo.to.hi1)
tail.lo.to.hi2 <- lmer(Tail.Width*2.54 ~ alternate + (1|clutch), data=lo.to.hi1) #lowest AIC
tail.lo.to.hi3 <- lmer(Tail.Width*2.54 ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(tail.lo.to.hi1, tail.lo.to.hi2, tail.lo.to.hi3)

summary(tail.lo.to.hi2)
car::Anova(tail.lo.to.hi2, type = "II") #switch NS
emmeans(tail.lo.to.hi2, pairwise ~ alternate)

### Tail Height models ###
morph <- read.csv("morph_Oct23.csv", header = T)

morph.c <- subset(morph, impact == "c")
tailheight1.c <- lmer(Tail.Height*2.54 ~ ba * cond + (1|tank), data=morph.c) #lowest AIC
tailheight2.c <- lmer(Tail.Height*2.54 ~ ba * cond + (1|clutch), data=morph.c)
tailheight3.c <- lmer(Tail.Height*2.54 ~ ba * cond + (1|tank) + (1|clutch), data=morph.c)
AIC(tailheight1.c, tailheight2.c, tailheight3.c)

summary(tailheight1.c)
car::Anova(tailheight1.c, type = "II") # ba, cond, and ba*cond S
emmeans(tailheight1.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(Tail.Height ~ ba*cond, data=morph.c)
morph.c$alternate <- factor(morph.c$alternate,
                            levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
boxplot(Tail.Height*2.54 ~ ba*cond, data=morph.c)
boxplot(Tail.Height*2.54 ~ alternate, data=morph.c, ylim = c(0.3,0.9), xlab ="Treatment Condition", ylab = "Tail Height (cm)", names = c("before low", "after low", "before high", "after high"))
stripchart(Tail.Height*2.54 ~ alternate, data = morph.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   
### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(morph, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(Tail.Height ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(Tail.Height*2.54 ~ alternate, data=hi.to.lo1, ylim = c(0.3,0.9), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Tail Height (cm)")
stripchart(Tail.Height*2.54 ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE) 

tailheight.hi.to.lo1 <- lmer(Tail.Height*2.54 ~ alternate + (1|tank), data=hi.to.lo1) #lowest AIC
tailheight.hi.to.lo2 <- lmer(Tail.Height*2.54 ~ alternate + (1|clutch), data=hi.to.lo1)
tailheight.hi.to.lo3 <- lmer(Tail.Height*2.54 ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(tailheight.hi.to.lo1, tailheight.hi.to.lo2, tailheight.hi.to.lo3)

summary(tailheight.hi.to.lo1)
car::Anova(tailheight.hi.to.lo1, type = "II") # switch S
emmeans(tailheight.hi.to.lo1, pairwise ~ alternate)


l.lo.to.hi1 <- subset(morph, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(Tail.Height ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(Tail.Height*2.54 ~ alternate, data = lo.to.hi1)
boxplot(Tail.Height*2.54 ~ alternate, data = lo.to.hi1, ylim = c(0.3,0.9), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Tail Height (cm)")
stripchart(Tail.Height*2.54 ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    


tailheight.lo.to.hi1 <- lmer(Tail.Height*2.54 ~ alternate + (1|tank), data=lo.to.hi1)
tailheight.lo.to.hi2 <- lmer(Tail.Height*2.54 ~ alternate + (1|clutch), data=lo.to.hi1) #lowest AIC
tailheight.lo.to.hi3 <- lmer(Tail.Height*2.54 ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(tailheight.lo.to.hi1, tailheight.lo.to.hi2, tailheight.lo.to.hi3)

summary(tailheight.lo.to.hi2)
car::Anova(tailheight.lo.to.hi2, type = "II") # switch S
emmeans(tailheight.lo.to.hi2, pairwise ~ alternate)



# PERFORMANCE -------------------------------------------------------------


### Bursts - Poisson ###
perf <- read.csv("PerfNov2.csv", header = T) 
summary(perf)
cor(perf[,c(7,11,13,14)], use="complete.obs")

perf.c <- subset(perf, impact == "c")
bursts.c1 <- glmer(bursts ~ ba * cond + (1|tank), data=perf.c, family="poisson") #lowest AIC
bursts.c2 <- glmer(bursts ~ ba * cond + (1|clutch), data=perf.c, family="poisson")
bursts.c3 <- glmer(bursts ~ ba * cond + (1|tank) + (1|clutch), data=perf.c, family="poisson")
AIC(bursts.c1, bursts.c2, bursts.c3)

bursts.nb1 <- glmer.nb(bursts ~ ba * cond + (1|tank), data=perf.c)
summary(bursts.nb1)
AIC(bursts.c1,bursts.nb1) # poisson lower AIC

summary(bursts.c1)
car::Anova(bursts.c1, type = "II") #none S
emmeans(bursts.c1, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(bursts ~ ba*cond, data=perf.c)
perf.c$alternate <- factor(perf.c$alternate,
                            levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(bursts ~ ba*cond, data=perf.c)
boxplot(bursts ~ alternate, data=perf.c, ylim = c(0,12), xlab ="Treatment Condition", ylab = "Bursts", names = c("before low", "after low", "before high", "after high"))
stripchart(bursts ~ alternate, data = perf.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   
### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(perf, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(bursts ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(bursts ~ alternate, data = hi.to.lo1)
boxplot(bursts ~ alternate, data=hi.to.lo1, ylim = c(0,12), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Bursts")
stripchart(bursts ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
bursts.hi.to.lo1 <- lmer(bursts ~ alternate + (1|tank), data=hi.to.lo1)
bursts.hi.to.lo2 <- lmer(bursts ~ alternate + (1|clutch), data=hi.to.lo1) #lowest AIC
bursts.hi.to.lo3 <- lmer(bursts ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(bursts.hi.to.lo1, bursts.hi.to.lo2, bursts.hi.to.lo3)

summary(bursts.hi.to.lo2)
car::Anova(bursts.hi.to.lo2, type = "II") #switch NS
emmeans(bursts.hi.to.lo2, pairwise ~ alternate)


lo.to.hi1 <- subset(perf, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(bursts ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(bursts ~ alternate, data = lo.to.hi1)
boxplot(bursts ~ alternate, data = lo.to.hi1, ylim = c(0,12), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Bursts")
stripchart(bursts ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    


bursts.lo.to.hi1 <- lmer(bursts ~ alternate + (1|tank), data=lo.to.hi1) # lowest AIC
bursts.lo.to.hi2 <- lmer(bursts ~ alternate + (1|clutch), data=lo.to.hi1) #lowest AIC
bursts.lo.to.hi3 <- lmer(bursts ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(bursts.lo.to.hi1, bursts.lo.to.hi2, bursts.lo.to.hi3)

summary(bursts.lo.to.hi1)
car::Anova(bursts.lo.to.hi1, type = "II") #switch NS
emmeans(bursts.lo.to.hi1, pairwise ~ alternate)

### Speed 1s ###
perf <- read.csv("PerfNov2.csv", header = T)

perf.c <- subset(perf, impact == "c")
speed1.c <- lmer(speed.1/100 ~ ba * cond + (1|tank), data=perf.c) #lowest AIC
speed2.c <- lmer(speed.1/100 ~ ba * cond + (1|clutch), data=perf.c)
speed3.c <- lmer(speed.1/100 ~ ba * cond + (1|tank) + (1|clutch), data=perf.c)
AIC(speed1.c, speed2.c, speed3.c)

summary(speed1.c)
car::Anova(speed1.c, type = "II") #ba S
emmeans(speed1.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(speed.1/100 ~ ba*cond, data=perf.c)
perf.c$alternate <- factor(perf.c$alternate,
                           levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(speed.1/100 ~ alternate, data=perf.c, ylim = c(0,0.25), xlab ="Treatment Condition", ylab = "Burst Speed (m/s)", names = c("before low", "after low", "before high", "after high"))
stripchart(speed.1/100 ~ alternate, data = perf.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)   
### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(perf, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(speed.1 ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(speed.1/100 ~ alternate, data=hi.to.lo1, ylim = c(0,0.25), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Burst Speed (m/s)")
stripchart(speed.1/100 ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

speed1.hi.to.lo1 <- lmer(speed.1/100 ~ alternate + (1|tank), data=hi.to.lo1) #lowest AIC
speed1.hi.to.lo2 <- lmer(speed.1/100 ~ alternate + (1|clutch), data=hi.to.lo1)
speed1.hi.to.lo3 <- lmer(speed.1/100 ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(speed1.hi.to.lo1, speed1.hi.to.lo2, speed1.hi.to.lo3)

summary(speed1.hi.to.lo1)
car::Anova(speed1.hi.to.lo1, type = "II") #switch S
emmeans(speed1.hi.to.lo1, pairwise ~ alternate)


l.lo.to.hi1 <- subset(perf, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(speed.1 ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(speed.1/100 ~ alternate, data = lo.to.hi1, ylim = c(0,0.25), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Burst Speed (m/s)")
stripchart(speed.1/100 ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

speed1.lo.to.hi1 <- lmer(speed.1/100 ~ alternate + (1|tank), data=lo.to.hi1) #lowest AIC
speed1.lo.to.hi2 <- lmer(speed.1/100 ~ alternate + (1|clutch), data=lo.to.hi1) #lowest AIC
speed1.lo.to.hi3 <- lmer(speed.1/100 ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(speed1.lo.to.hi1, speed1.lo.to.hi2, speed1.lo.to.hi3)

summary(speed1.lo.to.hi1)
car::Anova(speed1.lo.to.hi1, type = "II") #switch NS
emmeans(speed1.lo.to.hi1, pairwise ~ alternate)

### Speed 0.5s ###
perf <- read.csv("PerfNov2.csv", header = T)

perf.c <- subset(perf, impact == "c")
speed.5.c <- lmer(speed.5 ~ ba * cond + (1|tank), data=perf.c)
summary(speed.5.c)
car::Anova(speed.5.c, type = "II") #ba S
emmeans(speed.5.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(speed.5 ~ ba*cond, data=perf.c)

### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(perf, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
kruskal.test(speed.5 ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(speed.5 ~ alternate, data = hi.to.lo1)

speed5.hi.to.lo1 <- lmer(speed.5 ~ alternate + (1|tank), data=hi.to.lo1)
summary(speed5.hi.to.lo1)
emmeans(speed5.hi.to.lo1, pairwise ~ alternate)


l.lo.to.hi1 <- subset(perf, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
kruskal.test(speed.5 ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(speed.5 ~ alternate, data = lo.to.hi1)

speed5.lo.to.hi1 <- lmer(speed.5 ~ alternate + (1|tank), data=lo.to.hi1)
summary(speed5.lo.to.hi1)
emmeans(speed5.lo.to.hi1, pairwise ~ alternate)

### distance ###
perf <- read.csv("PerfNov2.csv", header = T)

perf.c <- subset(perf, impact == "c")
dist.c <- lmer(first.dist/100 ~ ba * cond + (1|tank), data=perf.c) #lowest AIC
dist.c2 <- lmer(first.dist/100 ~ ba * cond + (1|clutch), data=perf.c)
dist.c3 <- lmer(first.dist/100 ~ ba * cond + (1|tank) + (1|clutch), data=perf.c)
AIC(dist.c, dist.c2, dist.c3)

summary(dist.c)
car::Anova(dist.c, type = "II") #ba S
emmeans(dist.c, pairwise ~ ba*cond) # this provides all the relevant comparisons (+ a couple of irrelevant comparisons)
# (before and after low)
# (before and after high)
# (before low and high)
# (after low and high)

boxplot(first.dist/100 ~ ba*cond, data=perf.c)
perf.c$alternate <- factor(perf.c$alternate,
                           levels = c('beforelow','afterlow','beforehigh','afterhigh'),ordered = TRUE)
boxplot(first.dist/100 ~ alternate, data=perf.c, ylim = c(0,1.25), xlab ="Treatment Condition", ylab = "Burst Distance (m)", names = c("before low", "after low", "before high", "after high"))
stripchart(first.dist/100 ~ alternate, data = perf.c,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)  

### compare switched to their non-switched counterparts

hi.to.lo1 <- subset(perf, switched == "n" & alternate == "afterhigh" | 
                      switched == "y" & alternate =="afterlow" | alternate == "beforehigh")
hi.to.lo1$alternate <- factor(hi.to.lo1$alternate,
                              levels = c('beforehigh','afterhigh','afterlow'),ordered = TRUE)
kruskal.test(first.dist ~ alternate, data = hi.to.lo1) # no catching up when exposed to good conditions
boxplot(first.dist ~ alternate, data = hi.to.lo1)
boxplot(first.dist/100 ~ alternate, data=hi.to.lo1, ylim = c(0,1.25), names = c("before high", "after high", "after low"), xlab="Treatment Condition", ylab="Burst Distance (m)")
stripchart(first.dist/100 ~ alternate, data = hi.to.lo1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

dist.hi.to.lo1 <- lmer(first.dist/100 ~ alternate + (1|tank), data=hi.to.lo1) #lowest AIC
dist.hi.to.lo2 <- lmer(first.dist/100 ~ alternate + (1|clutch), data=hi.to.lo1)
dist.hi.to.lo3 <- lmer(first.dist/100 ~ alternate + (1|tank) + (1|clutch), data=hi.to.lo1)
AIC(dist.hi.to.lo1, dist.hi.to.lo2, dist.hi.to.lo3)

summary(dist.hi.to.lo1)
car::Anova(dist.hi.to.lo1, type = "II") #switch S
emmeans(dist.hi.to.lo1, pairwise ~ alternate)


lo.to.hi1 <- subset(perf, switched == "n" & alternate == "afterlow" | 
                        switched == "y" & alternate =="afterhigh" | alternate == "beforelow")
lo.to.hi1$alternate <- factor(lo.to.hi1$alternate,
                              levels = c('beforelow','afterlow','afterhigh'),ordered = TRUE)
kruskal.test(first.dist ~ alternate, data = lo.to.hi1) # good conditions gives them a leg up to withstand high cond (maybe at a cost to size)
boxplot(first.dist/100 ~ alternate, data = lo.to.hi1, ylim = c(0,1.25), names = c("before low", "after low", "after high"), xlab="Treatment Condition", ylab="Burst Distance (m)")
stripchart(first.dist/100 ~ alternate, data = lo.to.hi1,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = 1,    # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)    

dist.lo.to.hi1 <- lmer(first.dist/100 ~ alternate + (1|tank), data=lo.to.hi1) #lowest AIC
dist.lo.to.hi2 <- lmer(first.dist/100 ~ alternate + (1|clutch), data=lo.to.hi1) #lowest AIC
dist.lo.to.hi3 <- lmer(first.dist/100 ~ alternate + (1|tank) + (1|clutch), data=lo.to.hi1)
AIC(dist.lo.to.hi1, dist.lo.to.hi2, dist.lo.to.hi3)

summary(dist.lo.to.hi1)
car::Anova(dist.lo.to.hi1, type = "II") #switch S
emmeans(dist.lo.to.hi1, pairwise ~ alternate)

