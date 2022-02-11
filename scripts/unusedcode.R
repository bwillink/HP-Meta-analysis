#Would this be the case without accounting for covariance in sampling errors?
  
rma.mv(g ~ 1 + sei,
       V = var.g,
       random = list( ~ 1 | Experiment, ~ 1 | ID),
       data = dat_1,
       method = "REML")

#The results are qualitatively similar and both suggest bias. Here, I show for models with and without sampling error covariance within experiments

rma.mv(g ~ 1 + var.g,
       V = var.g,
       random = list( ~ 1 | Experiment, ~ 1 | ID),
       data = dat_1,
       method = "REML")

#Using the robust variance estimator instead of the variance-covariance error matrix within studies

Q1.B.R <-
  rma.mv(
    g ~ Trait.type:Gradient.category + var.g -1,
    V = var.g,
    random = list( ~ 1 | ID, ~ 1 | Experiment),
    data = dat_1,
    method = "REML"
  )

Q1.B.R_est <- clubSandwich::coef_test(Q1.B.R, vcov = "CR2")

forest.default(x= Q1.B.R_est$beta, sei =  Q1.B.R_est$SE, ci.lb = (Q1.B.R_est$beta - 1.96 * Q1.B.R_est$SE) , ci.ub =(Q1.B.R_est$beta + 1.96 * Q1.B.R_est$SE), 
               annotate=TRUE, showweights=T, header=F, 
               slab = c("variance","Environment fecundity", "Environment survivorship", "Pollution fecundity", "Pollution survivorship", "Resource fecundity", "Resource survivorship"))

#DATA set 2

#Without the var-covar matrix

# meta regression with SE
rma.mv(g ~ 1 + sei,
       V = var.g,
       random = list( ~ 1 | Experiment, ~ 1 | ID),
       data = dat_2,
       method = "REML")

#Using the robust variance estimator instead of the variance-covariance error matrix within studies

Q2.B.R <-
  rma.mv(
    g ~ Trait.type:Gradient.category + var.g -1,
    V = var.g,
    random = list( ~ 1 | ID, ~ 1 | Experiment),
    data = dat_2,
    method = "REML"
  )

Q2.B.R_est <- clubSandwich::coef_test(Q2.B.R, vcov = "CR2")

forest.default(x= Q2.B.R_est$beta, sei =  Q2.B.R_est$SE, ci.lb = (Q2.B.R_est$beta - 1.96 * Q2.B.R_est$SE) , ci.ub =(Q2.B.R_est$beta + 1.96 * Q2.B.R_est$SE), 
               annotate=TRUE, showweights=T, header=F, 
               slab = c("variance","Environment fecundity", "Environment intensity", "Environment prevalence", "Environment survivorship", "Pollution fecundity", "Pollution intensity", "Pollution prevalence", "Pollution survivorship", "Resource fecundity", "Resource intensity", "Resource prevalence", "Resource survivorship"))


### Weighting the effects of (small) studies
#A related question that one could want to check is how different studies contribute to the overall effects to see. In other words is our main result being driven by a few studies? Here I compute and plot a few measures of outlier/influence diagnostics.

#I'm also calculating the influence statistics cooks distance and dfbetas at the level of studies. In the metafor man pages: 

#"Cook's distance can be interpreted as the Mahalanobis distance between the entire set of predicted values once with the 𝑖th case included and once with the 𝑖th case excluded from the model fitting." 

#"the DFBETAS value(s) essentially indicate(s) how many standard deviations the estimated coefficient(s) change(s) after excluding the 𝑖th case from the model fitting."

#A value is influential if "The hat value is larger than 3×(𝑝/𝑘)." where p is the number of coefficients and k the number of cases.

#### Data set 1: fitness effects of environmental stress

# Cook's distance
Cook_dat1 <- cooks.distance(Q1, progbar=TRUE, dat_1$Study, reestimate=TRUE, parallel="multicore", ncpus = 2)

# how are these distances distributes?
hist(Cook_dat1)

# Which are the most influential studies?
Cook_dat1[which(Cook_dat1 > 1)]

# dfbetas
dfbetas_dat1 <- dfbetas.rma.mv(Q1, progbar=TRUE, dat_1$Study, reestimate=TRUE, parallel="multicore")

# Which are the most influential studies?
dfbetas_dat1[abs(dfbetas_dat1[,1])>1 | abs(dfbetas_dat1[,2])>1
             | abs(dfbetas_dat1[,3])>1 | abs(dfbetas_dat1[,4])>1
             | abs(dfbetas_dat1[,5])>1 | abs(dfbetas_dat1[,6])>1,]

# hat values
hat_dat1 <- hatvalues(Q1) 

# Which are the most influential effects?
inf_cutoff <- 3 * length(Q1$beta)/nrow(dat_1)

dat_1$Study[hat_dat1 > inf_cutoff]

#### Data set 2: resistance vs fitness effects of environmental stress

# Cook's distance
Cook_dat2 <- cooks.distance(Q2, progbar=TRUE, dat_2$Study, reestimate=TRUE, parallel="multicore", ncpus = 4)

# how are these distances distributes?
hist(Cook_dat2)

# Which are the most influential studies?
Cook_dat2[which(Cook_dat2 > 3)]

# dfbetas
dfbetas_dat2 <- dfbetas.rma.mv(Q2, progbar=TRUE, dat_2$Study, reestimate=TRUE, parallel="multicore")

# Which are the most influential studies?
dfbetas_dat2[abs(dfbetas_dat2[,1])>1 | abs(dfbetas_dat2[,2])>1
             | abs(dfbetas_dat2[,3])>1 | abs(dfbetas_dat2[,4])>1
             | abs(dfbetas_dat2[,5])>1 | abs(dfbetas_dat2[,6])>1,]

# hat values
hat_dat2 <- hatvalues(Q2) 

# Which are the most influential effects?
inf_cutoff <- 3 * length(Q2$beta)/nrow(dat_2)

dat_2$Study[hat_dat2 > inf_cutoff]

#### Heterogeneity

#If the variance-covariance matrix is included following the metafor example, total I2 becomes nearly 100%

# we account for covariance in sampling errors by using V, but assuming equal heterogeneity between fecundity vs survival and between environment, pollution and resources.

# matrix of sampling variances and covariances W 
W <- solve(varcovmat_1_PD$mat)

# model matrix
X <- model.matrix(Q1)

# estimate heterogeneity
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

# overall I2: how much of total variance is attributed to heterogeneity
100 * sum(Q1$sigma2) / (sum(Q1$sigma2) + (Q1$k-Q1$p)/sum(diag(P)))

# Heterogeneity within experiments, between experiments, between parasites
100 * Q1$sigma2 / (sum(Q1$sigma2) + (Q1$k-Q1$p)/sum(diag(P)))

sav <- confint(Q1)

# confidence intervals for I2
100 * sav[[1]]$random[1,2:3] /(sum(Q1$sigma2) + (Q1$k-Q1$p)/sum(diag(P))) ### CI for ID-level I^2
100 * sav[[2]]$random[1,2:3] / (sum(Q1$sigma2)+ (Q1$k-Q1$p)/sum(diag(P))) ### CI for the experiment-level I^2
#100 * sav[[3]]$random[1,2:3] / (sum(Q1$sigma2) + (Q1$k-Q1$p)/sum(diag(P)))  ### CI for the parasite-level I^2

#If the variance-covariance matrix is included following the metafor example, total I2 becomes nearly 100%

# we account for covariance in sampling errors by using V, but assuming equal heterogeneity between fecundity vs survival and between environment, pollution and resources.

# matrix of sampling variances and covariances W 
W <- solve(varcovmat_2_PD$mat)

# model matrix
X <- model.matrix(Q2.B)

# estimate heterogeneity
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

# overall I2: how much of total variance is attributed to heterogeneity
100 * sum(Q2.B$sigma2) / (sum(Q2.B$sigma2) + (Q2.B$k-Q2.B$p)/sum(diag(P)))

# Heterogeneity between experiments, within experiments, between parasites
100 * Q2.B$sigma2 / (sum(Q2.B$sigma2) + (Q2.B$k-Q2.B$p)/sum(diag(P)))

sav2 <- confint(Q2.B)

# confidence intervals for I2
100 * sav2[[1]]$random[1,2:3] /(sum(Q2.B$sigma2) + (Q2.B$k-Q2.B$p)/sum(diag(P))) ### CI for ID-level I^2
100 * sav2[[2]]$random[1,2:3] / (sum(Q2.B$sigma2)+ (Q2.B$k-Q2.B$p)/sum(diag(P))) ### CI for the experiment-level I^2
#100 * sav2[[3]]$random[1,2:3] / (sum(Q2.B$sigma2) + (Q2.B$k-Q2.B$p)/sum(diag(P)))  ### CI for the parasite-level I^2

# Descriptive statistics
d_stats <-
  data.frame(
    stat = c(
      "Number of papers",
      "Number of effects",
      "Number of experiments",
      "Number of arthropod species",
      "Number of mollusc species",
      "Number of fish species",
      "Number of amphibian species",
      "Number of reptile species",
      "Number of bird species",
      "Number of mammal species",
      "Number of viral taxa",
      "Number of bacterial taxa",
      "Number of Multiple infection",
      "Number of fungal taxa",
      "Number of protozan taxa",
      "Number of helminth taxa",
      "Number of environment experiments",
      "Number of pollution experiments",
      "Number of resource experiments",
      "Number of resistance effects",
      "Number of prevalence effects",
      "Number of intensity effects",
      "Number of fitness effects",
      "Number of survival effects",
      "Number of fecundity effects"
    ),
    Question_1 = c(
      length(unique(dat_1$Study)),
      length(dat_1$g),
      length(unique(dat_1$Experiment)),
      length(unique(dat_1$Host[which(dat_1$Host.type == "Arthropod")])),
      length(unique(dat_1$Host[which(dat_1$Host.type == "Mollusc")])),
      length(unique(dat_1$Host[which(dat_1$Host.type == "Fish")])),
      length(unique(dat_1$Host[which(dat_1$Host.type == "Amphibian")])),
      length(unique(dat_1$Host[which(dat_1$Host.type == "Reptile")])),
      length(unique(dat_1$Host[which(dat_1$Host.type == "Bird")])),
      length(unique(dat_1$Host[which(dat_1$Host.type == "Mammal")])),
      length(unique(dat_1$Parasite[which(dat_1$Parasite.type ==
                                           "Virus")])),
      length(unique(dat_1$Parasite[which(dat_1$Parasite.type ==
                                           "Bacteria")])),
      length(unique(dat_1$Parasite[which(dat_1$Parasite.type ==
                                           "Multiple")])),
      length(unique(dat_1$Parasite[which(dat_1$Parasite.type ==
                                           "Fungus")])),
      length(unique(dat_1$Parasite[which(dat_1$Parasite.type ==
                                           "Protozoan")])),
      length(unique(dat_1$Parasite[which(dat_1$Parasite.type ==
                                           "Helminth")])),
      length(dat_1$Experiment[which(dat_1$Gradient.category ==
                                      "Environment")]),
      length(dat_1$Experiment[which(dat_1$Gradient.category ==
                                      "Pollution")]),
      length(dat_1$Experiment[which(dat_1$Gradient.category ==
                                      "Resource")]),
      length(dat_1$ID[which(dat_1$Trait.category == "Epidemiological")]),
      length(dat_1$ID[which(dat_1$Trait.type == "Prevalence")]),
      length(dat_1$ID[which(dat_1$Trait.type == "Intensity")]),
      length(dat_1$ID[grep("demographic", dat_1$Trait.category)]),
      length(dat_1$ID[which(dat_1$Trait.type == "Survivorship")]),
      length(dat_1$ID[which(dat_1$Trait.type == "Fecundity")])
    ),
    Question_2 = c(
      length(unique(dat_2$Study)),
      length(dat_2$g),
      length(unique(dat_2$Experiment)),
      length(unique(dat_2$Host[which(dat_2$Host.type == "Arthropod")])),
      length(unique(dat_2$Host[which(dat_2$Host.type == "Mollusc")])),
      length(unique(dat_2$Host[which(dat_2$Host.type == "Fish")])),
      length(unique(dat_2$Host[which(dat_2$Host.type == "Amphibian")])),
      length(unique(dat_2$Host[which(dat_2$Host.type == "Reptile")])),
      length(unique(dat_2$Host[which(dat_2$Host.type == "Bird")])),
      length(unique(dat_2$Host[which(dat_2$Host.type == "Mammal")])),
      length(unique(dat_2$Parasite[which(dat_2$Parasite.type ==
                                           "Virus")])),
      length(unique(dat_2$Parasite[which(dat_2$Parasite.type ==
                                           "Bacteria")])),
      length(unique(dat_2$Parasite[which(dat_2$Parasite.type ==
                                           "Multiple")])),
      length(unique(dat_2$Parasite[which(dat_2$Parasite.type ==
                                           "Fungus")])),
      length(unique(dat_2$Parasite[which(dat_2$Parasite.type ==
                                           "Protozoan")])),
      length(unique(dat_2$Parasite[which(dat_2$Parasite.type ==
                                           "Helminth")])),
      length(dat_2$Experiment[which(dat_2$Gradient.category ==
                                      "Environment")]),
      length(dat_2$Experiment[which(dat_2$Gradient.category ==
                                      "Pollution")]),
      length(dat_2$Experiment[which(dat_2$Gradient.category ==
                                      "Resource")]),
      length(dat_2$ID[which(dat_2$Trait.category == "Epidemiological")]),
      length(dat_2$ID[which(dat_2$Trait.type == "Prevalence")]),
      length(dat_2$ID[which(dat_2$Trait.type == "Intensity")]),
      length(dat_2$ID[grep("demographic", dat_2$Trait.category)]),
      length(dat_2$ID[which(dat_2$Trait.type == "Survivorship")]),
      length(dat_2$ID[which(dat_2$Trait.type == "Fecundity")])
    )
  )

d_stats %>%
  kbl() %>%
  kable_material(c("striped", "hover"), full_width = F)


