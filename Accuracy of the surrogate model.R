# Preamble:

library(readxl)
library(rkriging)

################################################################################


##  Fragestellung: accuraly of the surrogate model
# test different kernel and basis functions

##  erstelle Ergebnis DF
kernels <- c("Gaussian", "Matern12", "Matern32", "Matern52", "RQ", "RQ", "RQ")

# Basisfunktionen PREP 
bf1 <- function(x) x
bf2 <- function(x) c(1,x)
bf3 <- function(x) x^2
bf4 <- function(x) c(1,x^2)
bf5 <- function(x) c(x,x^2)
bf6 <- function(x) c(1,x,x^2)
bf7 <- function(x){
  out <- x
  n <- length(x)
  for(i in 1:n){
    for(j in i:n){
      out <- c(out, x[i]*x[j])
    }
  }
  out
}
bf8 <- function(x){
  out <- c(1,x)
  n <- length(x)
  for(i in 1:n){
    for(j in i:n){
      out <- c(out, x[i]*x[j])
    }
  }
  out
}
bf9 <- function(x) log(x)
bf10 <- function(x) c(1, log(x))
bf11 <- function(x) sqrt(x)
bf12 <- function(x) c(1, sqrt(x))

bfs <- c(bf1, bf2, bf3, bf4, bf5, bf6, bf7, bf8, bf9, bf10, bf11, bf12)

MSPEs <- Bias <- data.frame(matrix(NA, 12, 7))
colnames(MSPEs) <-  colnames(Bias) <-  kernels
rownames(MSPEs) <- rownames(Bias) <- c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6", 
                                       "bf7", "bf8", "bf9", "bf10", "bf11", "bf12")


## Data Prep
data <- readxl::read_excel("Design of Experiments_out.xlsx")
data <- as.matrix(data)

dat <- data[, 2:7]
# scaling of data
#f_mem kann so bleiben
dat[,"sigma_mem"] <- dat[,"sigma_mem"]/1000 #umwandeln in MPa
dat[,"E_mem"] <- dat[,"E_mem"] / 1000000 #umwandeln in GPa
#nu_mem kann so bleiben
dat[,"sigma_edg"] <- dat[,"sigma_edg"] / 1000000 #umwandeln in GPa (ist in der Tabelle mit MPa angegeben!)
dat[,"sigma_sup"] <- dat[,"sigma_sup"] / 1000000 #umwandeln in GPa (ist in der Tabelle mit MPa angegeben!)

X <- dat 
y <- data[,8]
p <- 6

## (Fruit) Loop
for(i in 1:length(bfs)){
  for(j in kernels){
    # setze Parameter fuer RQ Kernel
    if(j %in% 1:4) a <- 0
    if(j == 5) a <- 0.5 
    if(j == 6) a <- 1 
    if(j == 7) a <- 2 
    
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

