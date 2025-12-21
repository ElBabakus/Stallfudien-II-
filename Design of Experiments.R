library(lhs)
library(evd)

## Gain necessary parameters of distributions with given expecting value and 
## standard deviation.


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



## LHS Sampling Implementation:

# Create LHS Samping from 7-dim unit cube:
set.seed(201125)
U <- optimumLHS(n = 200, k = 7, eps = 0.1, verbose = FALSE)


# Now use corresponding quantile functions to 
# X1: sigma_mem_y: Lognormal(11000, 1650)
# X2: f_mem: Gumbel(0.4, 0.12)
# X3: sigma_mem: Lognormal(4000,800)
# X4: E_mem: Lognormal(600000,90000)
# X5: nu_mem: Uniform(0.4, 0.011)
# X6: sigma_edg: Lognormal(353678, 70736)
# X7: sigma_sup: Lognormal(400835, 80167)

X <- cbind(
  qlnorm(U[,1], meanlog = lognormal_params(11000,1650)$mu,
         sdlog = lognormal_params(11000,1650)$sigma),
  qgumbel(U[,2], loc = gumbel_params(0.4, 0.12)$mu,
          scale = gumbel_params(0.4, 0.12)$beta),
  qlnorm(U[,3], meanlog = lognormal_params(4000,800)$mu,
         sdlog = lognormal_params(4000,800)$sigma),
  qlnorm(U[,4], meanlog = lognormal_params(600000,90000)$mu,
         sdlog = lognormal_params(600000,90000)$sigma),
  qunif(U[,5], min = unif_min_max(0.4, 0.0115470053837925)$a,
        max = unif_min_max(0.4, 0.0115470053837925)$b),
  qlnorm(U[,6], meanlog = lognormal_params(353677.651315323,70735.5302630646)$mu,
         sdlog = lognormal_params(353677.651315323,70735.5302630646)$sigma),
  qlnorm(U[,7], meanlog = lognormal_params(400834.671490699,80166.9342981399)$mu,
         sdlog = lognormal_params(400834.671490699,80166.9342981399)$sigma)
)


# Visualisation

curve(dlnorm(x, meanlog = lognormal_params(11000,1650)$mu, 
             sdlog = lognormal_params(11000,1650)$sigma), from = 0, to = 20000
      , col = "darkgreen", lwd = 2, ylab = "Lognormal (11000,1650)", 
      main = expression(sigma[mem_y]), xlab = "kPa")
points(x = X[,1], y = dlnorm(X[,1], meanlog = lognormal_params(11000,1650)$mu, 
                             sdlog = lognormal_params(11000,1650)$sigma), cex = 1)



curve(dgumbel(x, loc = gumbel_params(0.4, 0.12)$mu, 
              scale = gumbel_params(0.4, 0.12)$beta), from = -2, to = 4,
      col = "darkgreen", lwd = 2, ylab = "Gumbel (0.4, 0.12)", 
      main = expression(f[mem]), xlab = "kPa")
points(x = X[,2], y = dgumbel(X[,2], loc = gumbel_params(0.4, 0.12)$mu, 
                              scale = gumbel_params(0.4, 0.12)$beta), cex = 1)


curve(dlnorm(x, meanlog = lognormal_params(4000,800)$mu, 
             sdlog = lognormal_params(4000,800)$sigma), from = 0, to = 10000
      , col = "darkgreen", lwd = 2, ylab = "Lognormal (4000,800)", 
      main = expression(sigma[mem]), xlab = "kPa")
points(x = X[,3], y = dlnorm(X[,3], meanlog = lognormal_params(4000,800)$mu, 
                             sdlog = lognormal_params(4000,800)$sigma), cex = 1)


curve(dlnorm(x, meanlog = lognormal_params(600000,90000)$mu, 
             sdlog = lognormal_params(600000,90000)$sigma), from = 0, to = 1000000
      , col = "darkgreen", lwd = 2, ylab = "Lognormal (60000,9000)", 
      main = expression(E[mem]), xlab = "kPa")
points(x = X[,4], y = dlnorm(X[,4], meanlog = lognormal_params(600000,90000)$mu, 
                             sdlog = lognormal_params(600000,90000)$sigma), cex = 1)


curve(dunif(x, min = unif_min_max(0.4, 0.0115470053837925)$a, 
            max = unif_min_max(0.4, 0.0115470053837925)$b),
      from = 0.35 , to = 0.45,
      col = "darkgreen", lwd = 2, ylab = "Uniform (0.4,0.011)",
      main = expression(nu[mem]), n = 10000) # n = 10k weil sonst Linien Schräg
points(x = X[,5], y = dunif(X[,5], min = unif_min_max(0.4, 0.0115470053837925)$a, 
                            max = unif_min_max(0.4, 0.0115470053837925)$b), cex = 1)


curve(dlnorm(x, meanlog = lognormal_params(353677.651315323,70735.5302630646)$mu, 
             sdlog = lognormal_params(353677.651315323,70735.5302630646)$sigma), 
      from = 0, to = 700000, 
      col = "darkgreen", lwd = 2, ylab = "Lognormal (353678,70736)", 
      main = expression(sigma[edg]), xlab = "kPa")
points(x = X[,6], y = dlnorm(X[,6], meanlog = lognormal_params(353677.651315323,70735.5302630646)$mu, 
                             sdlog = lognormal_params(353677.651315323,70735.5302630646)$sigma), cex = 1)



curve(dlnorm(x, meanlog = lognormal_params(400834.671490699,80166.9342981399)$mu, 
             sdlog = lognormal_params(400834.671490699,80166.9342981399)$sigma), 
      from = 0, to = 800000, 
      col = "darkgreen", lwd = 2, ylab = "Lognormal (400835,80167)", 
      main = expression(sigma[sup]), xlab = "kPa")
points(x = X[,7], y = dlnorm(X[,7], meanlog = lognormal_params(400834.671490699,80166.9342981399)$mu, 
                             sdlog = lognormal_params(400834.671490699,80166.9342981399)$sigma), cex = 1)



colnames(X) = c("sigma_mem_y", "f_mem", "sigma_mem", 
                "E_mem", "nu_mem", "sigma_edg", "sigma_sup")
rownames(X) = 1:200

# Exportiere den Datensatz in xlsx
#writexl::write_xlsx(as.data.frame(X), "Design of Experiments.xlsx")
