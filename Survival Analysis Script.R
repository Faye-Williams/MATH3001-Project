library(survival)
library(survminer)

install.packages("ggplot2")
install.packages("ggpubr")
install.packages("dplyr")
install.packages("cmprsk")

library(ggplot2)
library(ggpubr)
library(cmprsk)
library(dplyr)
attach(ovarian)
data(cancer, package="survival")

#Simple model without plots

library(survival)
library(tidyverse)
library(survminer)

# Run Kaplan-Meier on the data
mod.ovarian <- survfit(Surv(futime, fustat) ~ 1, data = ovarian)

# Kaplan-Meier Survival Curve
ggsurvplot(mod.ovarian)

# Cumulative Hazard
ggsurvplot(mod.ovarian, fun = function(y) -log(y))

install.packages(rms)
library(rms)
mod.lung <- psm(Surv(futime, fustat) ~ 1, data = ovarian)
survplot(mod.lung, what="hazard")

s <- Surv(futime)

model <- survfit(s ~ 1)
summary(model)
plot(model)



#Survival curve
plot(survfit(Surv(futime, fustat) ~ 1, data = ovarian),
           xlab = "Days",
           ylab = "Overall survival probability")

#Summary
summary(survfit(Surv(futime, fustat) ~ 1, data = ovarian), type="fh")

#NA cum hazard plot
vfit <- survfit(Surv(futime, fustat) ~ 1, data = ovarian, type="fh")
plot(vfit, cumhaz=TRUE, xlab="Days", ylab="Cumulative hazard")

#PLOT KM WITH BOTH TREATMENT GROUPS

rx_fit <- survfit(Surv(futime, fustat) ~ rx, data=ovarian)
plot(rx_fit)


#KM AND NA ALTERNATE METHOD
with(ovarian, plot(Surv(futime, fustat), fun = "cumhaz", 
                main = "Cumulativa hazards function",
                xlab = "Duration"))
with(ovarian, plot(Surv(futime, fustat),
                main = "Survival function",
                xlab = "Duration"))


plot(survfit(Surv(futime, fustat) ~ 1, data = ovarian),
     xlab = "Days",
     ylab = "Overall survival probability")

######VETERANS#####

plot(survfit(Surv(time, status) ~ 1, data = veteran),
     xlab = "Days",
     ylab = "Overall survival probability")

#Summary
summary(survfit(Surv(time, status) ~ 1, data = veteran), type="fh")

#NA cum hazard plot
vfit <- survfit(Surv(time, status) ~ 1, data = veteran, type="fh")
plot(vfit, cumhaz=TRUE, xlab="Days", ylab="Cumulative hazard")

#PLOT KM WITH BOTH TREATMENT GROUPS

trt_fit <- survfit(Surv(time, status) ~ trt, data=veteran)
plot(trt_fit)


ggsurvplot(rx_fit)

ggsurvplot(trt_fit, conf.int=TRUE, pval=TRUE, 
          legend.labs=c("1", "2"), legend.title="Treatment Group",  
          palette=c("red", "blue"), 
          title="Kaplan-Meier Curve for Lung Cancer Survival in Veteran data set",
          surv.median.line = "hv")

ggsurvplot(
  trt_fit,
  fun = "event",
  linetype = "strata", # Change line type by groups
  pval = TRUE,
  conf.int = TRUE,
  ggtheme = theme_bw(),
  palette = c("red", "blue"),
  title = "Kaplan-Meier Cumulative Risk Function Estimate
  of lung cancer survival in veteran data")

ggsurvplot(
  trt_fit,
  fun = "cumhaz",
  linetype = "strata", # Change line type by groups
  pval = TRUE,
  conf.int = TRUE,
  ggtheme = theme_bw(),
  palette = c("red", "blue"),
  title = "Kaplan-Meier Cumulative Hazard Function Estimate
  of lung cancer survival in veteran data")

km_fit_vet <- survfit(Surv(time, status) ~ trt + celltype +
                          karno + diagtime + age + prior, data = veteran)


survdiff(formula = Surv(time, status) ~ trt, data = veteran)

#Thus p value is 92.2%, as our p value is 0.93

#The Chi-Squared test statistic is 0 with 1 degree of freedom and the 
#corresponding p-value is 0.9. Since this p-value is greater than 0.05, 
#we accept the null hypothesis.

#In other words, we do not have sufficient evidence to say that there is a 
#statistically significant difference in survival between the two groups.

survfit(formula = Surv(time, status) ~ trt, data = veteran)

install.packages("eha")
library(eha)
logrank(Surv(time, status), group = trt, data = veteran)


#ACCEPT THE SURVIVAL IS SAME

#LOGRANK FOR OVARIAN

survdiff(formula = Surv(futime, fustat) ~ rx, data = ovarian)

logrank(Surv(futime, fustat), group = rx, data = ovarian)

#The Chi-Squared test statistic is 1.1 with 1 degree of freedom and the 
#corresponding p-value is 0.3. Since this p-value is greater than 0.05, 
#we accept the null hypothesis.

#In other words, we do not have sufficient evidence to say that there is a 
#statistically significant difference in survival between the two groups.

#ACCEPT THE SURVIVAL IS SAME

#COX PH FOR VET

cox_veteran <- coxph(Surv(time, status) ~ trt + celltype + karno                   
             + diagtime + age + prior , data = veteran)
summary(cox_veteran)

#COX PH FOR OVARIAN

cox_ovarian <- coxph(Surv(futime, fustat) ~ age + resid.ds +
                       rx + ecog.ps, data = ovarian)
summary(cox_ovarian)

#
library(ggfortify)
install.packages(ggfortify)
library(ggplot2)

aa_fit <-aareg(formula = Surv(time, status) ~ trt + celltype +
                 karno + diagtime + age + prior , 
               data = veteran)
autoplot(aa_fit)

#CODE FOR WEIBULL NEED TO SWITCH UP
#https://bookdown.org/mpfoley1973/data-sci/survival-curve-estimation.html
data.frame(t = rep(1:80, 3),
           alpha = c(rep(1.5, 80), rep(1, 80), rep(0.75, 80)),
           lambda = rep(0.03, 240)) %>%
  mutate(
    f = dweibull(x = t, shape = alpha, scale = 1 / 0.03),
    S = pweibull(q = t, shape = alpha, scale = 1 / 0.03, lower.tail = FALSE),
    h = f / S  # same as alpha * lambda^alpha * t^(alpha-1)
  ) %>%
  ggplot(aes(x = t, y = h, color = as.factor(alpha))) +
  geom_line() +
  theme(legend.position = "top") +
  labs(y = "hazard", x = "time", color = "alpha",
       title = "Weibul hazard function at varying levels of alpha",
       subtitle = "Lambda = 0.03",
       caption = "alpha = 1 is special case of exponential function.")

#https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
# Veteran deets - https://stat.ethz.ch/R-manual/R-devel/library/survival/html/veteran.html

#bkrn infoo - https://si.biostat.washington.edu/sites/default/files/modules//SISCR_2016_4_all.pdf

#https://bioconnector.github.io/workshops/r-survival.html#kaplan-meier_plots
#https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
#http://ehar.se/r/ehar2/proportional-hazards-and-cox-regression.html#the-log-rank-test
#https://www.r-bloggers.com/2021/08/log-rank-test-in-r-survival-curve-comparison/

