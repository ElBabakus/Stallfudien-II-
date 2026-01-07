library(rkriging)
library(evd)

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
# scaling of data
#f_mem kann so bleiben
dat[,"sigma_mem"] <- dat[,"sigma_mem"]/1000 #umwandeln in MPa
dat[,"E_mem"] <- dat[,"E_mem"] / 1000000 #umwandeln in GPa
#nu_mem kann so bleiben
dat[,"sigma_edg"] <- dat[,"sigma_edg"] / 1000000 #umwandeln in GPa (ist in der Tabelle mir MPa angegeben!)
dat[,"sigma_sup"] <- dat[,"sigma_sup"] / 1000000 #umwandeln in GPa (ist in der Tabelle mir MPa angegeben!)

X <- dat 
y <- data[,8]
p <- 6

#### Modell Fitten ----------------------------------------

## mit dem besten Model und der besten Basisfunktion bf3: Matern52 + function(x) x^2
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
mean((y - Pred)^2) # 77.35707

## Bias
mean(y - Pred) # -0.01867302


#### Fehlerwkt schätzen -----------------------------------
n <- 1e7

## ziehe fuer jede Variable n Zufallszahlen
set.seed(25122025)
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

# scaling of data
#f_mem kann so bleiben
X_star[,2] <- X_star[,2]/1000 #umwandeln in MPa
X_star[,3] <- X_star[,3] / 1000000 #umwandeln in GPa
#nu_mem kann so bleiben
X_star[,5] <- X_star[,5] / 1000000 #umwandeln in GPa (ist in der Tabelle mir MPa angegeben!)
X_star[,6] <- X_star[,6] / 1000000 #umwandeln in GPa (ist in der Tabelle mir MPa angegeben!)


## ziehe n mal den maximal erlaubten stress
set.seed(25122025)
sigma_mem_y <- rlnorm(n, meanlog = lognormal_params(11000,1650)$mu,
                      sdlog = lognormal_params(11000,1650)$sigma)

## berechne den output an den stellen X_star
f_star <- Predict.Kriging(mod, X_star)$mean

par(cex.axis = 1.5)
par(cex.lab = 1.5)
par(mar = c(5,7,2,2))
hist(f_star, breaks=100, main = "", xlab = "kPa", ylab =  expression(sigma[mem_max]),
     freq = FALSE, xlim = c(4000, 18000), ylim = c(0,7e-4), yaxs = "i", col = "grey")
  curve(dlnorm(x, meanlog = lognormal_params(11000,1650)$mu, 
             sdlog = lognormal_params(11000,1650)$sigma), from = 0, to = 20000
      , col = "red", lwd = 1.5, add = TRUE)
legend(
  x = 14539.43,
  y = 0.0006674010,
  legend = expression(sigma[mem_y]),
  col = "red",
  lwd = 1.5,
  bty = "n",
  cex = 1.5
)
box()

## wie oft ist der output größer als erlaubt?
mean(f_star > sigma_mem_y)


## TODO Montecarlo Fehler berechnen -> dabei kommt es mMn darauf an wie y
##      skaliert ist
