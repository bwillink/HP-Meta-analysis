
### Note that this works correctly but produces nonsensical output right now because it uses a dummy column for the
# Hedge's g and var(g) columns
setwd("C:/RData")
dat = read.csv("Multcomp_test.csv")
head(dat)

# Dummy labeling for the Hedge's g and var(g) columns - ***Current output is meaningless***
# Chose "Experiment" column for this because it is numerical
gs = dat[,"Experiment"]
var_gs = dat[,"Experiment"]

# This code produces a variance-covariance of sampling errors matrix, n x n, with n = number of effect sizes
# requires the dataset (dat), Hedge's gs (gs), and var(g) (var_gs) to be defined above
# Uses "Experiment", "Level_C", "N_C", and "N_X" columns of AV and BW's dataset

varcovmat = matrix(0, nrow = dim(dat)[1], ncol = dim(dat)[1])
for (i in 1:dim(dat)[1]) {
  for (j in 1:dim(dat)[1]) {
    if (i == j) {varcovmat[i,j] = var_gs[i]}else {
      if (dat[i, "Experiment"] == dat[j, "Experiment"] & dat[i, "Level_C"] == dat[j, "Level_C"]) {
        varcovmat[i,j] = 1/dat[i,"N_C"] + gs[i]*gs[j]/(dat[i,"N_C"] + dat[j,"N_C"] + dat[j,"N_X"]) 
      }
    }
  }
}

                   