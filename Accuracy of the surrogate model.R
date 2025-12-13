# Preamble:

library(plotly)
library(rgl)
library(evd)
library(randtoolbox)
library(lhs)
library(DiceDesign)
library(writexl)
library(readxl)
library(rkriging)
setwd("C:/Users/ferry/OneDrive/Dokumente/__UNI/9. Semester/Stallfudien 2/Projekt 2")
################################################################################


##  Fragestellung: accuraly of the surrogate model
# test different kernel and basis functions

# Kernels
kernels <- c("Gaussian", "Matern12", "Matern32", "Matern52", "RQ")

# Basisfunktionen PREP 
bf1 <- function(x) x
bf2 <- function(x) x^2
bf3 <- function(x) c(x, x^2)
bf4 <- function(x){
  out <- x
  n <- length(x)
  for(i in 1:n){
    for(j in 1:n){
      out <- c(out, x[i]*x[j])
    }
  }
  out
}
bf5 <- function(x) c(1,x)
bf6 <- function(x) c(1,x^2)
bf7 <- function(x) c(1,x,x^2)

bfs <- c(bf1, bf2, bf3, bf4, bf5, bf6, bf7)

MSPEs <- Bias <- data.frame(matrix(NA, 7, 5))
colnames(MSPEs) <-  colnames(Bias) <-  kernels
rownames(MSPEs) <- rownames(Bias) <- c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6", "bf7")


## Data Prep
data <- readxl::read_excel("Design of Experiments_out.xlsx")
data <- as.matrix(data)

dat <- data[, 2:7]
dat <- scale(dat)

X <- dat 
y <- data[,8]
p <- 6

## (Fruit) Loop
for(i in 1:length(bfs)){
  for(j in kernels){
    kernel.parameters <- list(
      type = j,
      lengthscale = rep(1, 6), # initial value
      lengthscale.lower.bound = rep(1e-3, p), # lower bound
      lengthscale.upper.bound = rep(1e3, p), # upper bound
      alpha = 2 # für den RQ Kernel
    )
    
    #### eine Naeherungsschaetzung fuer die lengthscale parameter
    mod <- Universal.Kriging(X, y, basis.function = bfs[[i]],
                             kernel.parameters = kernel.parameters,
                             nlopt.parameters = list(maxeval = 1000)
    )
    
    ## Update kernel.parameters
    kernel.parameters$lengthscale <- mod$get_lengthscale()
    
    #### Kreuvalidierung
    Pred <- numeric(200)
    
    for(k in 1:200) {
      mod <- Universal.Kriging(X[-k, ], y[-k], basis.function = bfs[[i]],
                               kernel.parameters = kernel.parameters,
                               nlopt.parameters = list(maxeval = 1000))
      
      Pred[k] <- Predict.Kriging(mod, t(X[k, ]))$mean
      
      if(k%%20 == 0) print(paste0("at ", Sys.time(), " i = ", i, ", j = ", j, ", k = ", k))
    }
    
    ## MSPE
    MSPEs[i,j] <- mean((y - Pred)^2)
    
    ## Bias
    Bias[i,j] <- mean(y - Pred)
  }
}

save(MSPEs, Bias, file = "Accuracy Results.RData")
