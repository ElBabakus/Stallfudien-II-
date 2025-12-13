library(rkriging)
library(ordinal)

#### lade Daten -------------------------------------------

#### params umrechnen
lognormal_params <- function(E,S){
  return(list(
    mu = log(E^2 / sqrt(E^2 + S^2)),
    sigma = sqrt(log(1 + S^2/E^2))
  ))
}

unif_min_max <- function(E,S){
  return(list(
    a = E - sqrt(3)*S,
    b = E + sqrt(3)*S
  ))
}

gumbel_params <- function(E,S){
  gamm <- -digamma(1) # Euler- Macheroni- Konstante
  
  return(list(
    mu = E - sqrt(6)*S/pi*gamm,
    beta = sqrt(6)*S/pi
  ))
}


#### Daten Vorbereiten -----------------------------------------

data <- readxl::read_excel("Design of Experiments_out.xlsx")
data <- as.matrix(data)

dat <- data[, 2:7]
dat <- scale(dat)

X <- dat 
y <- data[,8]
p <- 6

#### Modell Fitten ----------------------------------------

## mit dem besten Model und der besten Basisfunktion: Matern52 + function(x) x^2
kernel.parameters <- list(
  type = "Matern52",
  lengthscale = rep(1, 6), # initial value
  lengthscale.lower.bound = rep(1e-3, p), # lower bound
  lengthscale.upper.bound = rep(1e3, p) # upper bound
)
basis.function <- function(x) x^2

# eine Naeherungsschaetzung fuer die lengthscale parameter
mod <- Universal.Kriging(X, y, basis.function = basis.function,
                         kernel.parameters = kernel.parameters,
                         nlopt.parameters = list(maxeval = 1000)
)
l <- mod$get_lengthscale()
l

## Update kernel.parameters
kernel.parameters$lengthscale <- mod$get_lengthscale()

#### Kreuvalidierung --------------------------------------
Pred <- numeric(200)
for(k in 1:200) {
  mod <- Universal.Kriging(X[-k, ], y[-k], basis.function = basis.function,
                           kernel.parameters = kernel.parameters,
                           nlopt.parameters = list(maxeval = 1000))
  
  Pred[k] <- Predict.Kriging(mod, t(X[k, ]))$mean
  
  if(k%%20 == 0) print(paste0("at ", Sys.time(), ", k = ", k))
}

## MSPE
mean((y - Pred)^2)

## Bias
mean(y - Pred)


#### Fehlerwkt schätzen -----------------------------------
n <- 1e7

## ziehe fuer jede Variable n Zufallszahlen
X_star <- cbind(
  rgumbel(n, loc = gumbel_params(0.4, 0.12)$mu,
          scale = gumbel_params(0.4, 0.12)$beta),
  rlnorm(n, meanlog = lognormal_params(4000,800)$mu,
         sdlog = lognormal_params(4000,800)$sigma),
  rlnorm(n, meanlog = lognormal_params(600000,90000)$mu,
         sdlog = lognormal_params(600000, 90000)$sigma),
  runif(n, min = unif_min_max(0.4, 0.0115470053837925)$a,
        max = unif_min_max(0.4, 0.0115470053837925)$b),
  rlnorm(n, meanlog = lognormal_params(353677.651315323, 70735.5302630646)$mu,
         sdlog = lognormal_params(353677.651315323, 70735.5302630646)$sigma),
  rlnorm(n, meanlog = lognormal_params(400834.671490699, 80166.9342981399)$mu,
         sdlog = lognormal_params(400834.671490699, 80166.9342981399)$sigma)
)


## teile X_star spaltenweise durch das maximum der erhobenen Daten (NO)
# X_star <- sweep(X_star, 2, STATS = maximums, FUN = "/")
X_star <- scale(X_star)

## ziehe n mal den maximal erlaubten stress
sigma_mem_y <- rlnorm(n, meanlog = lognormal_params(11000,1650)$mu,
                      sdlog = lognormal_params(11000,1650)$sigma)

## berechne den output an den stellen X_star
f_star <- Predict.Kriging(mod, X_star)$mean
hist(f_star, breaks=100)

## wie oft ist der output größer als erlaubt?
mean(f_star > sigma_mem_y)


## TODO Montecarlo Fehler berechnen -> dabei kommt es mMn darauf an wie y
##      skaliert ist
