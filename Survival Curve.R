library(survival)
library(maxstat)
library(survminer)
df <- read.csv("Hybrid Survival Curve.csv", sep=",", header=T)
df.cut <- surv_cutpoint(df, time = "Disease.Free.Survival.Days..Overall.", event = "Survival.Code", variables = "Hybrid",
              minprop = 0.1, progressbar = TRUE)
plot(df.cut, "Hybrid", palette = "npg")
df.cat <- surv_categorize(df.cut)
df.cat$Disease.Free.Survival.Days..Overall. <- as.numeric(df.cat$Disease.Free.Survival.Days..Overall.)
fit <- survfit(Surv(Disease.Free.Survival.Days..Overall.,Survival.Code) ~Hybrid, data = df.cat)
ggsurvplot(fit, data = df.cat, risk.table = TRUE, conf.int = TRUE,pval = T)

