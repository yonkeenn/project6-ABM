require(deSolve)


diff_eqs2 <- function(t, x, params){
  
  with(as.list(c(x, params)), {
    # rates
    r1 <- (beta1*I1 + beta2*I2 + beta3*I3)*S
    r2 <- (a)*E
    r3 <- (gama1 + p1)*I1
    r4 <- (p1)*I1
    r5 <- (gama2 + p2)*I2
    r6 <- (p2)*I2
    r7 <- (gama3 + u)*I3
    r8 <- (gama1*I1 + gama1*I2 + gama1*I3)
    r9 <- (u)*I3
    
    # derivatives
    dS <- -r1
    dE <- r1 - r2
    dI1 <- r2 - r3
    dI2 <- r4 - r5
    dI3 <- r6 - r7
    dR <- r8
    dD <- r9
    
    
    #output
    return(list(c(dS, dE, dI1, dI2, dI3, dR, dD)))
    
  })
  
}

pop.size <- 1000
IO <- 9
xstart <- c(S=pop.size-IO, E=9, I1=0, I2=0, I3=0, R=0, D=0)
params <- c(beta1= 0.0002, beta2 = 0.0002, beta3 = 0.0002, a = 0.23, gama1 = 0.08, gama2 = 0.068, gama3 = 0.085, p1 = 0.1, p2 = 0.022, u = 0.057)
times <- seq(0, 52, by=0.01)

out <- as.data.frame(ode(xstart, times, diff_eqs2, params))


plot.diff_eqs2 <- function(out){
  par(mfrow=c(3,2))
  plot(out$time, out$S, xlab='Time',ylab='Susceptible', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$E, xlab='Time',ylab='Exposed', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$I1, xlab='Time',ylab='Infectious1', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$I2, xlab='Time',ylab='Infectious2', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$I3, xlab='Time',ylab='Infectious3', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$R, xlab='Time',ylab='Recovered', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$D, xlab='Time',ylab='Died', type='l', col='steelblue4', lwd=3)
}

## Analisis de Incertidumbre y Sensibilidad

parametros <- c("beta1", "beta2", "beta3", "a", "gama1", "gama2", "gama3", "p1", "p2", "u")
parametros_dist <- c("qnorm", "qnorm", "qnorm", "qnorm", "qnorm", "qnorm", "qnorm", "qnorm", "qnorm", "qnorm")
parametros_dist.arg <- list( list(mean=0.5, sd=0.15), 
               list(mean=0.1, sd=0.25), 
               list(mean=0.1, sd=0.13), 
               list(mean=0.2, sd=0.27),
               list(mean=0.133, sd=0.41), 
               list(mean=0.125, sd=0.25),
               list(mean=0.075, sd=0.39), 
               list(mean=0.033, sd=0.33),
               list(mean=0.042, sd=0.27), 
               list(mean=0.05, sd=0.11)
               )


SEIR_LHS_TEST <- LHS(model=NULL, parametros, 50, parametros_dist, parametros_dist.arg)

# Guardando los Datos:
write.csv(get.data(SEIR_LHS_TEST), file="SEIR_LHS_TEST.csv")

# Mis Resultado
#SEIR_MODEL_Result <- read.table("SEIR_MODEL_Result.csv",sep = ",",header = TRUE ) 
SEIR_MODEL_Result <- read.table("SEIR_MODEL_Result.csv") 


# Juntando los parametros
#acopladaLHS <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result)

# Ejecutando PRCC
#plotprcc(acopladaLHS)


#####

acopladaLHS_1 <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result[,2])
acopladaLHS_1 <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result[,1])
acopladaLHS_2 <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result[,2])
acopladaLHS_3 <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result[,3])
acopladaLHS_4 <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result[,4])
acopladaLHS_5 <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result[,5])
acopladaLHS_6 <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result[,6])
acopladaLHS_7 <- tell(SEIR_LHS_TEST, SEIR_MODEL_Result[,7])


plotprcc(acopladaLHS_1)
plotprcc(acopladaLHS_2)
plotprcc(acopladaLHS_3)
plotprcc(acopladaLHS_6)
plotprcc(acopladaLHS_7)
