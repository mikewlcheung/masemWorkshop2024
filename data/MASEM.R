## -----------------------------------------------------------------------------
## Required packages in this workshop
lib2install <- c("metaSEM", "lavaan", "semPlot", "readxl")

## Install them automatically if they are not available on your computer
for (i in lib2install) {
  if (!(i %in% rownames(installed.packages()))) install.packages(i)
}


















## ----message=FALSE------------------------------------------------------------
## Load the library for MASEM
library(metaSEM)

## Load the library to read XLSX file
library(readxl)

## Read the study characteristics
study <- read_xlsx("Digman97.xlsx", sheet="Info")

## Display a few studies
head(study)
 
## Create an empty list to store the correlation matrices
Digman97.data <- list()
  
## Read 1 to 14 correlation matrices
for (i in 1:14) {
  ## Read each sheet and convert it into a matrix
  mat <- as.matrix(read_xlsx("Digman97.xlsx", sheet=paste0("Study ", i)))
  ## Add the row names
  rownames(mat) <- colnames(mat)
  ## Save it into a list
  Digman97.data[[i]] <- mat
}

## Add the names of the studies
names(Digman97.data) <- study$Study

## Show the first few studies
head(Digman97.data)

## Extract the sample sizes
Digman97.n <- study$n
Digman97.n

## Extract the cluster
Digman97.cluster <- study$Cluster
Digman97.cluster

## Display the no. of studies
pattern.na(Digman97.data, show.na=FALSE)

## Display the cumulative sample sizes
pattern.n(Digman97.data, Digman97.n)


## -----------------------------------------------------------------------------
model <- "## Factor loadings
          ## Alpha is measured by A, C, and ES
          Alpha =~ A + C + ES
          ## Beta is measured by E and I
          Beta =~ E + I
          ## Factor correlation between Alpha and Beta
          Alpha ~~ Beta"

## Display the model
plot(model, color="yellow")

## Convert the lavaan syntax into a RAM model as the metaSEM only knows the RAM model
## It is important to ensure that the variables are arranged in A, C, ES, E, and I.
RAM <- lavaan2RAM(model, obs.variables=c("A","C","ES","E","I"))
RAM




## -----------------------------------------------------------------------------
## method="FEM": fixed-effects TSSEM
fixed1 <- tssem1(Digman97.data, Digman97.n, method="FEM")

## summary of the findings
summary(fixed1)

## extract coefficients
coef(fixed1)


## -----------------------------------------------------------------------------
fixed2 <- tssem2(fixed1, RAM=RAM)
summary(fixed2)

plot(fixed2, color="green")


## -----------------------------------------------------------------------------
# Display the original study characteristic
table(Digman97.cluster)     

## Younger participants: "Children" and "Adolescents"
## Older participants: "Mature adults"
sample <- ifelse(Digman97.cluster %in% c("Children", "Adolescents"), 
                 yes="Younger participants", no="Older participants")

table(sample)

## cluster: variable for the analysis with cluster
fixed1.cluster <- tssem1(Digman97.data, Digman97.n, method="FEM", cluster=sample)

summary(fixed1.cluster)


## -----------------------------------------------------------------------------
fixed2.cluster <- tssem2(fixed1.cluster, RAM=RAM)

summary(fixed2.cluster)

## Setup two plots side-by-side
layout(t(1:2))

## Plot the first group
plot(fixed2.cluster[[1]], col="green")
title("Younger participants")

## Plot the second group
plot(fixed2.cluster[[2]], col="yellow")
title("Older participants")




## -----------------------------------------------------------------------------
## method="REM": Random-effects model
random1 <- tssem1(Digman97.data, Digman97.n, method="REM", RE.type="Diag")

summary(random1)

## Extract the fixed-effects estimates
(est_fixed <- coef(random1, select="fixed"))

## Average correlation matrix
## Convert the estimated vector to a symmetrical matrix
## where the diagonals are fixed at 1 (for a correlation matrix)
averageR <- vec2symMat(est_fixed, diag=FALSE)
dimnames(averageR) <- list(c("A", "C", "ES", "E", "I"),
                           c("A", "C", "ES", "E", "I"))
averageR




## ----warning=FALSE------------------------------------------------------------
random2 <- tssem2(random1, RAM=RAM)
summary(random2)

## Plot the parameter estimates
plot(random2, color="green")


## ----message=FALSE------------------------------------------------------------
library(readxl)
library(metaSEM)

## Read the Excel file
df <- read_xlsx("Nohe15.xlsx", sheet="Table A1")

head(df)

## Variable names
my.var <- c("W1", "S1", "W2", "S2")

## Split each row as a list
my.list <- split(df, 1:nrow(df))

## Select the correlation coefficients and convert it into a correlation matrix
my.cor <- lapply(my.list, 
                 function(x) {## Convert the correlations into a correlation matrix
                              mat <- vec2symMat(unlist(x[c("W1-S1","W1-W2","W1-S2","S1-W2","S1-S2","W2-S2")]), diag = FALSE)
                              ## Add the dimensions for ease of reference
                              dimnames(mat) <- list(my.var, my.var)
                              ## Return the correlation matrix
                              mat})

## Add the study names for ease of reference
names(my.cor) <- df$Study

head(my.cor)

## Sample sizes
df$N

## Display the no. of studies
pattern.na(my.cor, show.na=FALSE)

## Display the cumulative sample sizes
pattern.n(my.cor, df$N)


## ----message=FALSE------------------------------------------------------------
## Fit the stage one model (fixed-effects model)
fix1 <- tssem1(my.cor, df$N, method="FEM")

summary(fix1)

## Fit the stage one model (random-effects model)
rand1 <- tssem1(my.cor, df$N, method="REM")

summary(rand1)

## Extract the fixed-effects estimates
R1 <- coef(rand1, select = "fixed")
R1

## Convert it to a correlation matrix
R1 <- vec2symMat(R1, diag = FALSE)
R1

## Add the dimension names for ease of reference
dimnames(R1) <- list(my.var, my.var)
R1


## -----------------------------------------------------------------------------
## Model in lavaan syntax
## The parameter labels are not necessary. But they make the output easier to follow.
model1 <- "W2 ~ w2w*W1 + s2w*S1    # path coefficeients   
           S2 ~ s2s*S1 + w2s*W1
           S1 ~~ w1withs1*W1       # correlation
           S2 ~~ w2withs2*W2
           W1 ~~ 1*W1              # variances of the independent variables are fixed at 1
           S1 ~~ 1*S1"

## Display the model
## See the layout argument in ?semPaths
## layout can be tree, circle, spring, tree2, circle2
plot(model1, color="yellow", layout="spring")

## Convert the above model to RAM specification used by metaSEM.
## We also have to specify how the variables are arranged in the data.
RAM1 <- lavaan2RAM(model1, obs.variables = my.var)
RAM1

## Fit the stage two model
rand2a <- tssem2(rand1, RAM=RAM1)

summary(rand2a)

plot(rand2a, color="green", layout="spring")


## -----------------------------------------------------------------------------
model2 <- "W2 ~ same*W1 + diff*S1     # regression coefficeients   
           S2 ~ same*S1 + diff*W1
           S1 ~~ w1cs1*W1             # correlation
           S2 ~~ w2cs2*W2
           W1 ~~ 1*W1                 # variances of the independent variables are fixed at 1
           S1 ~~ 1*S1"

## Display the model
plot(model2, color="yellow")

## Convert the above model to RAM specification used by metaSEM.
## We also have to specify how the variables are arranged in the data.
RAM2 <- lavaan2RAM(model2, obs.variables = my.var)

## Fit the stage two model
rand2b <- tssem2(rand1, RAM=RAM2)

summary(rand2b)

## Comparing nested models with anova()
anova(rand2a, rand2b)

plot(rand2b, color="green")


## -----------------------------------------------------------------------------
## Convert the correlation matrices into a dataframe, which is required in osmasem()
Nohe.df <- Cor2DataFrame(x=my.cor, n=df$N)

## Standardize the moderator "Lag" to improve numerical stability and add it into the dataframe
Nohe.df$data$Lag <- scale(df$Lag)
    
head(Nohe.df$data)

## Fit the model without any moderator
## The results are similat to that in TSSEM.
fit0 <- osmasem(model.name="No moderator", RAM=RAM1, data=Nohe.df)
summary(fit0)

plot(fit0)

## Create a matrix to present "Lag" moderator
## We use Lag to moderate the following paths:
## W1 -> W2
## W1 -> S2
## S1 -> W2
## S1 -> S2
A1 <- create.modMatrix(RAM=RAM1, output="A", mod="Lag")
A1

## Fit the model with Lag as a moderator
fit1 <- osmasem(model.name="Lag as a moderator", RAM=RAM1, Ax=A1, data=Nohe.df)
summary(fit1)

## Comparing the models with and without the covariate
anova(fit1, fit0)

## Effect of w2w when Lag is at the mean value
mxEval(w2w, fit1$mx.fit)

## Effect of w2w when Lag is at +1 SD as Lag is standardized
mxEval(w2w + w2w_1, fit1$mx.fit)

## Effect of w2w when Lag is at -1 SD as Lag is standardized
mxEval(w2w - w2w_1, fit1$mx.fit)




## -----------------------------------------------------------------------------
library(metaSEM)

## Read the SPSS dataset
my.df <- foreign::read.spss("Schutte21.sav", use.value.labels = TRUE, to.data.frame=TRUE)

## A function to convert rows into a 3x3 correlation matrix
create.matrix <- function(x, type=c(1, 2, 3)) {
  mat <- matrix(NA, ncol=3, nrow=3)
  diag(mat) <- 1
  type <- as.character(type)
  ## Mindfulness, EI, Gratitude
  ## 1: Mindfulness and EI
  ## 2: Mindfulness and Gratitude
  ## 3: EI and Gratitude
  switch(type,
         "1" = mat[1, 2] <- mat[2, 1] <- unlist(x),
         "2" = mat[1, 3] <- mat[3, 1] <- unlist(x),
         "3" = mat[2, 3] <- mat[3, 2] <- unlist(x))
  mat
}

## Convert the data to correlation matrices using the create.matrix().
my.cor <- lapply(split(my.df, seq(nrow(my.df))),
                 function(x, y) create.matrix(x["Effect_size"], x["Type_of_Association"]))

## Variable names
varlist <- c("Mindfulness", "EI", "Gratitude")

## Add the variable names in the dimnames.
my.cor <- lapply(my.cor, function(x) {dimnames(x) <- list(varlist, varlist); x}  )

## Add the study names
names(my.cor) <- my.df$Study_name

## Correlation matrices in the analysis
head(my.cor)

## Sample sizes
my.n <- my.df$N
my.n

## Number of studies in each cell
pattern.na(my.cor, show.na = FALSE)

## Total sample sizes in each cell
pattern.n(my.cor, my.n)


## -----------------------------------------------------------------------------
## Stage 1 analysis: find an average correlation matrix
stage1 <- tssem1(my.cor, my.n)
summary(stage1)

## Average correlation matrix
meanR <- vec2symMat(coef(stage1, select = "fixed"), diag = FALSE)
dimnames(meanR) <- list(varlist, varlist)
meanR

## Absolute heterogeneity variance: tau^2
tau2 <- vec2symMat(coef(stage1, select = "random"), diag = FALSE)
dimnames(tau2) <- list(varlist, varlist)
tau2

## Relative heterogeneity index: I^2
I2 <- vec2symMat(summary(stage1)$I2.values[, "Estimate"], diag = FALSE)
dimnames(I2) <- list(varlist, varlist)
I2


## ----message=FALSE------------------------------------------------------------
## Proposed model
model <- "## cp: c prime (c'), the common notation for the direct effect.
          Gratitude ~ cp*Mindfulness + b*EI
          EI ~ a*Mindfulness
          Mindfulness ~~ 1*Mindfulness
          ## Define direct, indirect, and total effects
          Direct := cp
          Indirect := a*b
          Total := a*b + cp"

plot(model, color="yellow")

RAM1 <- lavaan2RAM(model, obs.variables = varlist)
RAM1

## Stage 2 analysis: fit the path model
## Likelihood-based CI: intervals.type = "LB"
stage2 <- tssem2(stage1, RAM=RAM1, intervals.type = "LB")
summary(stage2)

plot(stage2, color="green")


## ----message=FALSE------------------------------------------------------------
## Proposed model
model <- "Gratitude ~ cp*Mindfulness + b*EI
          EI ~ a*Mindfulness
          Mindfulness ~~ 1*Mindfulness
          ## Define direct, indirect, and total effects
          Direct := cp
          Indirect := a*b
          Total := a*b + cp
          ## Apply a constraint on the direct and indirect effects
          cp == a*b"

RAM2 <- lavaan2RAM(model, obs.variables = varlist)
RAM2

stage2b <- tssem2(stage1, RAM=RAM2, intervals.type = "LB")
summary(stage2b)

## Compare the models with and without c=a*b
anova(stage2, stage2b)


## ----eval=FALSE---------------------------------------------------------------
## library(metaSEM)
## 
## ## Use all cores-2, i.e., keep 2 cores for other things.
## mxOption(key='Number of Threads', value=parallel::detectCores()-2)


## -----------------------------------------------------------------------------
## X1 is okay
X1 <- matrix(c(1, 0.8, 0.8, 1), nrow=2)
X1
is.pd(X1)

## X2 is not okay as the correlation is 1.2.
X2 <- matrix(c(1, 1.2, 1.2, 1), nrow=2)
X2
is.pd(X2)

## X3 is not okay as the correlation matrix is not symmetric.
X3 <- matrix(c(1, 0.6, 0.8, 1), nrow=2)
X3
isSymmetric(X3)






## -----------------------------------------------------------------------------
pattern.na(Hunter83$data, show.na=FALSE)

