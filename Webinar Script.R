#Contemporary Methods in Causal Inference for Program Evaluation Webinar
# by Noah Greifer | noah@unc.edu
# University of North Carolina at Chapel Hill
# 11/1/2018

#This example estimates the effect of mothers having a job
#that provides or subsidizes child care
#on the length that they breastfeed their children
#National Longitudinal Survey of Youth 1979 (NLSY79)
#and the NLSY79 Children and Youth

#install the relevant packages
install.packages(c("cobalt", "MatchIt", "Matching",
                 "rgenoud", "optmatch", "designmatch",
                 "WeightIt", "optweight", "jtools"))

#load dataset
load("NLSY79.rData")

#define propensity score formula
cov.names <- c(
    "childRace", #RACE OF CHILD (MOTHER'S RACIAL/ETHNIC COHORT FROM SCREENER)                CRACE
    "hrsWorked4thQtr", #USUAL HOURS WORKED BY MOTHER AT ALL JOBS IN 4TH QTR BEFORE BIRTH OF CHILD  NHRJBS01
    "earnings4thQtr", #USUAL EARNINGS OF MOTHER AT ALL JOBS IN 4TH QTR BEFORE BIRTH OF CHILD      EARNJB01
    "postBirthLeaveWks", #weeks after birth that mother returned to work (restricted to 12 weeks or less)
    "highestGrade" , #HIGHEST GRADE COMPLETED AS OF MAY 1 SURVEY YEAR (REVISED)
    "residenceRegion", #REGION OF CURRENT RESIDENCE
    "ruralUrban"  , #IS R'S CURRENT RESIDENCE URBAN/RURAL?
    "totalWelfare",  # TOTAL AMOUNT AFDC FOOD STAMPS OR OTH WELFARE/SSI RECEIVED DURING CAL YR
    "weeksWorked", #NUMBER OF WEEKS WORKED IN PAST CALENDAR YEAR
    "hoursPerWeek", #hours per week worked in the last calendar year
    "maternityLeave", #FRINGE BENEFITS CURRENT JOB/MOST RECENT JOB - MATERNITY/PATERNITY LEAVE
    "flexibleSchedule", #job allows flexible schedule
    "argueChores", #FREQUENCY R & HUSBAND/PARTNER ARGUE ABOUT-CHORES & RESPONSIBILITIES
    "dentalInsurance", #company provides dental insurance
    "lifeInsurance", #company provides life insurance
    "profitSharing", #company provides profit sharing
    "retirement", #company provided retirement plan
    "training") #company provided training opportunities

#obtain the propensity score formula
library(cobalt)

cov.form <- f.build("childCare", cov.names)

bal.tab(cov.form, data = NLSY79.data, estimand = "ATT",
        m.threshold = .1)

#Matching

#Historical methods
library(MatchIt)

#1:1 nearest neighbor without replacement, logistic regression PS
m.out1 <- matchit(cov.form, data = NLSY79.data,
                  method = "nearest",
                  distance = "logit",
                  replace = FALSE,
                  ratio = 1)
(b1 <- bal.tab(m.out1, m.threshold = .1))

#1:3 nearest neighbor without replacement, probit PS
m.out2 <- matchit(cov.form, data = NLSY79.data,
                  method = "nearest",
                  distance = "probit",
                  replace = FALSE,
                  ratio = 3)
(b2 <- bal.tab(m.out2, m.threshold = .1))

#Machine Learning

#Genetic matching
m.out3 <- matchit(cov.form, data = NLSY79.data,
                  method = "genetic",
                  distance = "logit",
                  replace = FALSE,
                  ratio = 1, 
                  pop.size = 150,
                  ks = FALSE,
                  print.level = 0)
(b3 <- bal.tab(m.out3, m.threshold = .1))

#Optimization hybrid methods

#1:1 optimal without replacement, logistic regression PS
m.out4 <- matchit(cov.form, data = NLSY79.data,
                  method = "optimal",
                  distance = "logit",
                  replace = FALSE,
                  ratio = 1)
(b4 <- bal.tab(m.out4, m.threshold = .1))

#Optimal full matching, logistic regression PS
m.out5 <- matchit(cov.form, data = NLSY79.data,
                  method = "full",
                  distance = "logit")
(b5 <- bal.tab(m.out5, m.threshold = .1))

#Optimization

#designmatch with balance constraints: 0 mean differences, 0 proportion differences
library(designmatch)
NLSY79.data.sorted <- NLSY79.data[order(NLSY79.data$childCare, decreasing = TRUE),]

cat.vars <- cov.names[sapply(NLSY79.data[,cov.names], function(x)
    is.factor(x) || length(unique(x)) <= 2)]
cont.vars <- setdiff(cov.names, cat.vars)

d.out1 <- bmatch(NLSY79.data.sorted$childCare,
                 n_controls = 1,
                 total_groups = sum(NLSY79.data.sorted$childCare),
                 mom = list(covs = NLSY79.data.sorted[,cont.vars],
                            tols = rep(0, length(cont.vars))),
                 fine = list(covs = NLSY79.data.sorted[,cat.vars]))
(b6 <- bal.tab(d.out1, formula = cov.form, data = NLSY79.data.sorted, 
               m.threshold = .1))

#Balance comparison
b <- bal.tab(cov.form, data = NLSY79.data,
             weights = list(NN1to1 = get.w(m.out1),
                            NN1to3 = get.w(m.out2),
                            Genetic = get.w(m.out3),
                            Optimal = get.w(m.out4),
                            Full = get.w(m.out5)),
             estimand = "ATT", 
             method = c(rep("matching", 4), "weighting"),
             m.threshold = .1)

print(b, disp.bal.tab = FALSE)

#Effect estimation
library(jtools)

#MatchIt output
summ(lm(wksBreastfed ~ childCare, data = NLSY79.data,
        weights = m.out4$weights))

#designmatch output
with(NLSY79.data.sorted,
     t.test(wksBreastfed[d.out1$t_id],
            wksBreastfed[d.out1$c_id]))

#Weighting

#Historical Methods

#PS Matching with logistic PS
library(WeightIt)
w.out1 <- weightit(cov.form, data = NLSY79.data,
                   estimand = "ATT", method = "ps")
(b6 <- bal.tab(w.out1, m.threshold = .1))

#Machine Learning

#GBM

w.out2 <- weightit(cov.form, data = NLSY79.data,
                   estimand = "ATT", method = "gbm",
                   stop.method = "es.max",
                   interaction.depth = 4,
                   verbose = FALSE)
(b7 <- bal.tab(w.out2, m.threshold = .1))

#Optimization Hybrids

#CBPS, over-identified

w.out3 <- weightit(cov.form, data = NLSY79.data,
                   estimand = "ATT", method = "cbps")
(b8 <- bal.tab(w.out3, m.threshold = .1))

#Optimization

#Entropy balancing

w.out4 <- weightit(cov.form, data = NLSY79.data,
                   estimand = "ATT", method = "ebal")
(b9 <- bal.tab(w.out4))

#optweights, perfect balance
library(optweight)
ow.out1 <- optweight(cov.form, data = NLSY79.data,
                     estimand = "ATT", tols = 0,
                     min.w = 1e-8)
(b10 <- bal.tab(ow.out1, estimand = "ATT"))

#optweights, approximate balance
ow.out2 <- optweight(cov.form, data = NLSY79.data,
                     estimand = "ATT", tols = .01,
                     min.w = 1e-8)
(b11 <- bal.tab(ow.out2, estimand = "ATT",
                m.threshold = .1))

#Effect estimation
library(jtools)

#WeightIt output
summ(lm(wksBreastfed ~ childCare, data = NLSY79.data,
        weights = w.out1$weights), robust = TRUE)

#optweight output
summ(lm(wksBreastfed ~ childCare, data = NLSY79.data,
        weights = ow.out1$weights), robust = TRUE)
