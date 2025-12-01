
surv <- read.csv("surv.csv", header=T)
summary(surv)
str(surv)
head(surv)
attach(surv)

library(survival)
library(survminer)

# set up the data using week as the time period (right-censored)

so <- Surv(time = days, event = FateCode, type='right')
so

# effect of site on survival

fit1 <- survfit(so ~ cond, data = surv)
summary(fit1)
surv_summary(fit1)

fit1

ggsurvplot(fit1, data = surv, pval = TRUE)
ggsurvplot(fit1, data = surv, fun = "event")
ggsurvplot(fit1, data = surv, fun = "cumhaz")

ggsurvplot(fit1, data = surv,
           conf.int = TRUE,
           xlab= "Days",
           ylab= "Survival Probability",
           pval = FALSE,
           fun = NULL,
           break.x.by = 50,
           risk.table = FALSE,
           size = 0.7,
           palette = c("darkgray","darkgray","black","black"),
           linetype = c("solid","dashed","solid","dashed"),
           legend = c(0.15,0.2),
           font.x = 14,
           font.y =14,
           font.legend = 14,
           legend.title = "",
           legend.labs = c("high", "high-to-low", "low",
                           "low-to-high"))
help(ggsurvplot)

# Pairwise survdiff
res <- pairwise_survdiff(Surv(days, FateCode) ~ cond,
                         data = surv)
res
summary(res)
