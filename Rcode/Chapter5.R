# Applied Multiple Imputation
# Advantages, Pitfalls, New Developments and Applications in R
#
# Kristian Kleinke, Jost Reinecke, Daniel Salfran, Martin Spiess
################################################################

# Chapter 5
#######################

#################################################################
# Note: we use relative file paths                              #
# For the examples to work, please set your R working directory #
# to the parent folder: "AppliedMI"                             #
#################################################################

# Chapter 5.3.1 Reading in the CrimoC data

dat <- read.table("./data/crim4.dat", 
                  col.names=c("id", "FEMALE", "RE", "GY","HA", "ACRIM", "BCRIM", "CCRIM", "DCRIM"), 
                  na.strings="-999", 
                  colClasses = c("numeric", "factor", "factor", "factor", "factor", "numeric", "numeric","numeric", "numeric"))

head(dat)

# Chapter 5.3.2 Missing data inspection

colSums(is.na(dat))

round(colSums(is.na(dat))/nrow(dat)*100,2)

library(mice)
md.pattern(dat)

md.pattern(dat)[,c(1:5,9,7,6,8,10)]

md.pairs(dat)$mr

delinq = dat[,c("ACRIM", "BCRIM", "CCRIM", "DCRIM")]
mdp = md.pairs(delinq)$mr / colSums(is.na(delinq))
round(mdp,2)

round( cor( data.matrix(dat[,2:9]), use = "pairwise.complete.obs"), 2)

R <- is.na(dat)
colnames(R) <- paste("R.", colnames(dat), sep = "")
dat2 <-  data.frame(R,dat)
round( cor(dat2[,6:9], data.matrix( dat2[,11:18]), use = "pairwise.complete.obs"), 2 )

# Chapter 5.4.3 Multiple imputation by predictive mean matching

dat2 <- dat[,c("ACRIM", "BCRIM", "CCRIM", "DCRIM","FEMALE", "RE", "GY")] 
ini <- mice(dat2, maxit = 0)
ini$method
ini$predictorMatrix
test.mice.pmm <- mice(dat2, m = 5, maxit = 20, printFlag = FALSE, seed=123456 )
plot(test.mice.pmm)

imp.mice.pmm1 <- mice(dat2, m = 100, maxit = 20, donors = 1 , print = FALSE, seed=123456)
imp.mice.pmm5 <- mice(dat2, m = 100, maxit = 20, print = FALSE, seed=123456)
imp.mice.pmm20 <- mice(dat2, m = 100, maxit = 20, donors = 20, print = FALSE, seed=123456)

method <-ini$method
method[1:4] <- "midastouch"

test.mice.mt = mice(dat2, method=method, m = 5, maxit = 20 , seed = 123456)
plot(test.mice.mt)

imp.mice.mt <- mice(dat2, method = method, m = 100, maxit = 20, printFlag = FALSE, seed = 123456)

source("./functions/prepare.mplus.R")
dir.create("./imputations/pmm/midastouch", recursive = TRUE)
prepare.mplus(imp.mice.mt, dir = "./imputations/pmm/midastouch")

dir.create("./imputations/pmm/k1", recursive = TRUE)
dir.create("./imputations/pmm/k5", recursive = TRUE)
dir.create("./imputations/pmm/k20", recursive = TRUE)

prepare.mplus(imp.mice.pmm1, dir = "./imputations/pmm/k1")
prepare.mplus(imp.mice.pmm5, dir = "./imputations/pmm/k5")
prepare.mplus(imp.mice.pmm20, dir = "./imputations/pmm/k20")

# Chapter 5.5.1 norm2

library(norm2)
acrim.z <- ifelse(dat$ACRIM!=0,1,0)  
acrim.c <- ifelse(dat$ACRIM==0,NA,dat$ACRIM)

bcrim.z <- ifelse(dat$BCRIM!=0,1,0)  
bcrim.c <- ifelse(dat$BCRIM==0,NA,dat$BCRIM)

ccrim.z <- ifelse(dat$CCRIM!=0,1,0)  
ccrim.c <- ifelse(dat$CCRIM==0,NA,dat$CCRIM)

dcrim.z <- ifelse(dat$DCRIM!=0,1,0)  
dcrim.c <- ifelse(dat$DCRIM==0,NA,dat$DCRIM)

# dat3 data set to be imputed

dat3 <- cbind(acrim.z, acrim.c,
              bcrim.z, bcrim.c,
              ccrim.z, ccrim.c,
              dcrim.z, dcrim.c,
              dat[,c("FEMALE", "RE", "GY")]) # completely observed predictors

head(dat3[,1:8])

set.seed(123456)
emResult <- 
emNorm(cbind( acrim.z, acrim.c, bcrim.z, bcrim.c, ccrim.z, ccrim.c , dcrim.z , dcrim.c ) ~ FEMALE + GY + RE, data = dat3)

summary(emResult, show.variables=FALSE, show.patterns=FALSE, show.params=FALSE)

table(dat3$acrim.z)

table(dat3$acrim.z, acrim.c, useNA="ifany")


set.seed(123456)
mcmcResult <- mcmcNorm(emResult,iter = 1000)

par(mfrow=c(4,2))
for(i in 1:8){
  acf(mcmcResult$series.beta[,i],lag.max=1000, main=colnames(dat3)[i])}

par(mfrow=c(4,2))
for(i in 1:8)
  plot(mcmcResult$series.beta[,i],ylab="", main=colnames(dat3)[i])

emResult2 <- emNorm(cbind(acrim.z, acrim.c, bcrim.z, bcrim.c, ccrim.z, ccrim.c, dcrim.z, dcrim.c)~ FEMALE+GY+RE, data = dat3, prior = "ridge", prior.df = 200)

summary(emResult2, show.variables=FALSE, show.patterns=FALSE, show.params=FALSE)

set.seed(123456)
mcmcResult2 <- mcmcNorm(emResult2, iter = 1000)


par(mfrow=c(4,2))
for(i in 1:8){
  acf(mcmcResult2$series.beta[,i],lag.max=1000, main=colnames(dat3)[i])}

par(mfrow=c(4,2))
for(i in 1:8)
  plot(mcmcResult2$series.beta[,i],ylab="", main=colnames(dat3)[i])

set.seed(123456)
mcmcResult_m100 <- mcmcNorm(emResult2, iter = 100000, impute.every = 1000)

dir.create("./imputations/other/norm2/")

for (m in 1:100){
  
  com=mcmcResult_m100$imp.list[[m]]
  com=data.frame(com)
  
  # round zero indicators to 0 and 1
  
  com$acrim.z <- round(com$acrim.z) 
  com$acrim.z[com$acrim.z<0] <- 0 
  com$acrim.z[com$acrim.z>1] <- 1
  
  com$bcrim.z <- round(com$bcrim.z) 
  com$bcrim.z[com$bcrim.z<0] <- 0 
  com$bcrim.z[com$bcrim.z>1] <- 1
  
  com$ccrim.z <- round(com$ccrim.z)
  com$ccrim.z[com$ccrim.z<0] <- 0
  com$ccrim.z[com$ccrim.z>1] <- 1
  
  com$dcrim.z <- round(com$dcrim.z)
  com$dcrim.z[com$dcrim.z<0] <- 0
  com$dcrim.z[com$dcrim.z>1] <- 1
  
  # round count variables to integer values
  # make sure they range between 1 and 16
  
  com$acrim.c <- round(com$acrim.c)
  com$acrim.c[com$acrim.c<1] <- 1
  com$acrim.c[com$acrim.c>16] <- 16
  
  com$bcrim.c <- round(com$bcrim.c)
  com$bcrim.c[com$bcrim.c<1] <- 1
  com$bcrim.c[com$bcrim.c>16] <- 16
  
  com$ccrim.c <- round(com$ccrim.c)
  com$ccrim.c[com$ccrim.c<1] <- 1
  com$ccrim.c[com$ccrim.c>16] <- 16
  
  com$dcrim.c <- round(com$dcrim.c)
  com$dcrim.c[com$dcrim.c<1] <- 1
  com$dcrim.c[com$dcrim.c>16] <- 16
  
  # replace 1s by the respective count
  
  counts <- which(com$acrim.z==1)
  com$ACRIM <- com$acrim.z
  com$ACRIM[counts] <- com$acrim.c[counts]
  
  counts <- which(com$bcrim.z==1)
  com$BCRIM <- com$bcrim.z
  com$BCRIM[counts] <- com$bcrim.c[counts]
  
  counts <- which(com$ccrim.z==1) 
  com$CCRIM <- com$ccrim.z
  com$CCRIM[counts] <- com$ccrim.c[counts]
  
  counts <- which(com$dcrim.z==1) 
  com$DCRIM <- com$dcrim.z
  com$DCRIM[counts]<- com$dcrim.c[counts]
  
  out <- cbind(1:nrow(com), dat2[,c("FEMALE", "RE", "GY")],
               com[,c("ACRIM", "BCRIM", "CCRIM", "DCRIM")])
  
  write.table( out , paste( "./imputations/other/norm2/imp" , m , ".dat" , sep="") ,
               quote=F , row.names=F , col.names= F)
}

for (m in 1:100) cat(paste( "imp" , m , ".dat" , sep=""), "\n", file="./imputations/other/norm2/implist.dat", append=TRUE)

# Chapter 5.5.2 Amelia

library(Amelia)
set.seed(123456)
imp.amelia <- amelia(x=dat3, m=100, noms=c("FEMALE","RE","GY"))

set.seed(123456)
imp.amelia.ords <- amelia(x = dat2, m = 100, ords = c("ACRIM","BCRIM","CCRIM","DCRIM"),noms = c("FEMALE","RE","GY"))

dir.create("./imputations/other/Amelia/noms", recursive=TRUE)
dir.create("./imputations/other/Amelia/ords", recursive=TRUE)

write.amelia(imp.amelia.ords, file.stem = "./imputations/other/Amelia/ords/imp", extension = ".dat", format = "table",
             quote = FALSE, row.names = FALSE, col.names = FALSE)

for (m in 1:100) cat(paste( "imp" , m , ".dat" , sep=""), "\n", file=paste0("./imputations/other/Amelia/ords/implist.dat"),append=TRUE)

plot(imp.amelia.ords, lwd = 2, col = c("black","gray80"))

# obtain original variable from zero and count indicators
# write them to disk for further analysis in Mplus

for (m in 1:100){
  
  com=imp.amelia$imputations[[m]]
  
  # round zero indicators to zero and one
  
  com$acrim.z=round(com$acrim.z) 
  com$acrim.z[com$acrim.z<0]<-0 
  com$acrim.z[com$acrim.z>1]<-1
  
  com$bcrim.z=round(com$bcrim.z) 
  com$bcrim.z[com$bcrim.z<0]<-0 
  com$bcrim.z[com$bcrim.z>1]<-1
  
  com$ccrim.z=round(com$ccrim.z)
  com$ccrim.z[com$ccrim.z<0]<-0
  com$ccrim.z[com$ccrim.z>1]<-1
  
  com$dcrim.z=round(com$dcrim.z)
  com$dcrim.z[com$dcrim.z<0]<-0
  com$dcrim.z[com$dcrim.z>1]<-1
  
  # round count variables to integer values
  # make sure they range between 1 and 16
  
  com$acrim.c=round(com$acrim.c)
  com$acrim.c[com$acrim.c<1]<-1
  com$acrim.c[com$acrim.c>16]<-16
  
  com$bcrim.c=round(com$bcrim.c)
  com$bcrim.c[com$bcrim.c<1]<-1
  com$bcrim.c[com$bcrim.c>16]<-16
  
  com$ccrim.c=round(com$ccrim.c)
  com$ccrim.c[com$ccrim.c<1]<-1
  com$ccrim.c[com$ccrim.c>16]<-16
  com$dcrim.c=round(com$dcrim.c)
  com$dcrim.c[com$dcrim.c<1]<-1
  com$dcrim.c[com$dcrim.c>16]<-16
  
  # replace 1's by the respective count
  
  counts=which(com$acrim.z==1)
  com$ACRIM=com$acrim.z
  com$ACRIM[counts]=com$acrim.c[counts]
  
  counts=which(com$bcrim.z==1)
  com$BCRIM=com$bcrim.z
  com$BCRIM[counts]=com$bcrim.c[counts]
  
  counts=which(com$ccrim.z==1) 
  com$CCRIM=com$ccrim.z
  com$CCRIM[counts]=com$ccrim.c[counts]
  
  counts=which(com$dcrim.z==1) 
  com$DCRIM=com$dcrim.z
  com$DCRIM[counts]=com$dcrim.c[counts]
  
  out=cbind(1:nrow(com), 
            com[,c("FEMALE", "RE", "GY", "ACRIM", "BCRIM", "CCRIM", "DCRIM")])
  
  write.table( out , paste( "./imputations/other/Amelia/noms/imp" , m , ".dat" , sep="") ,
               quote=F , row.names=F , col.names= F)
}


for (m in 1:100)
{
  cat(paste( "imp" , m , ".dat" , sep=""), "\n",
      file="./imputations/other/Amelia/noms/implist.dat",append=TRUE)
}

# Chapter 5.5.3 mi

library(mi)
mdf <- missing_data.frame(dat2)
show(mdf)

mdf2 <- change(mdf, y = c("ACRIM", "BCRIM", "CCRIM", "DCRIM"), what = "type", to = "count")

mdf3 <- change(mdf, y = c("ACRIM", "BCRIM", "CCRIM", "DCRIM"), what = "method", to = "pmm")

set.seed(123456)
imp.mi <- mi(mdf3, n.chains = 100, parallel=FALSE)		

summary(imp.mi)$ACRIM

imputed.data <- mi::complete(imp.mi, m = 100)

# prepare data for analysis in Mplus

dir.create("./imputations/other/mi")
for (m in 1:100){
  com=imputed.data[[m]]
  out=cbind(1:nrow(com), 
            com[,c("FEMALE", "RE", "GY", "ACRIM", "BCRIM", "CCRIM", "DCRIM")])
  
  write.table( out , paste( "./imputations/other/mi/imp" , m , ".dat" , sep="") ,
               quote=F , row.names=F , col.names= F)
}

for (m in 1:100) cat(paste( "./imputations/other/mi/imp" , m , ".dat" , sep=""), "\n", file="implist.dat",append=TRUE)

# Chapter 5.5.4 aregImpute

dir.create("./imputations/other/aregImpute")
set.seed(123456)
imp.aregImpute <- aregImpute(~ ACRIM + BCRIM + CCRIM + DCRIM + FEMALE + RE + GY,
                             data = dat2,
                             type = "normpmm",
                             nk = 0,
                             n.impute = 100)	

for (m in 1:100){
  com <- cbind.data.frame(impute.transcan(imp.aregImpute, imputation = m, data = dat2,
                                          list.out = TRUE, pr = FALSE, check = FALSE))
  out <- cbind(1:nrow(com), com[,c("FEMALE", "RE", "GY", "ACRIM", "BCRIM", "CCRIM", "DCRIM")])
  write.table( out , paste( "./imputations/other/aregImpute/imp" , m , ".dat" , sep="") ,
               quote = FALSE , row.names = FALSE , 
               col.names = FALSE)
}

for (m in 1:100) cat(paste( "./imputations/other/aregImpute/imp" , m , ".dat" , sep=""), "\n", file="implist.dat",append=TRUE)

# Chapter 5.6.3.1 Multilevel multiple imputation in mice

dat.long <- reshape(data = dat2, varying = list(c("ACRIM","BCRIM","CCRIM","DCRIM")), direction = "long", v.names = "DELINQ", times = 0:3, timevar = "TIME", idvar = "ID")
head(dat.long)

library(mice)
ini <- mice(dat.long, maxit = 0)
method <- ini$method
method["DELINQ"] <- "2l.pan"
pred <- ini$pred
pred["DELINQ",] <- c(1,1,1,2,0,-2)

set.seed(123456)
test.imp.2l.pan <- mice(dat.long, method = method, pred = pred, m = 5, printFlag = FALSE)
plot(test.imp.2l.pan)

set.seed(123456)
imp.2l.pan <- mice(dat.long, method = method, pred = pred, m = 100, printFlag = FALSE)

dat.long$TIME2 <- dat.long$TIME^2

ini <- mice(dat.long, maxit = 0)
method <- ini$method
pred <- ini$pred
method["DELINQ"] <- "2l.pan"
pred["DELINQ",] <- c(1,1,1,2,0,-2,2)
print(pred)

set.seed(123456)
imp.2l.pan.quadratic <- mice(dat.long, method = method, pred = pred, m = 100, printFlag = FALSE)


# Chapter 5.6.3.2 pan – Multiple imputation under the multivariate linear mixed-effects model

library(pan)
dir.create("./plots")

set.seed(123456)
pan.convergence <- test.pan(group = dat.long$ID, y = matrix(dat.long$DELINQ, ncol = 1), x = dat.long[,c("TIME","FEMALE","RE","GY")], xcol = 1:5, zcol = 1:2, paniter = 5000)

iter = 5000
pdf(file = "./plots/panTS.pdf")
par(mfrow = c(4,2))
plot(1:iter, pan.convergence$sigma[1,1,],type="l", xlab="Iteration", ylab="Sigma",main="")
plot(1:iter, pan.convergence$psi[1,1,], type="l", xlab="Iteration", ylab="Psi (Intercept)")
plot(1:iter, pan.convergence$psi[2,1,], type="l", xlab="Iteration", ylab="Psi (Slope)")
for (i in 1:5){
  plot(1:iter, pan.convergence$beta[i,1,],type="l",
       xlab="Iteration",
       ylab=paste("Beta ",i-1))
}
dev.off()
pdf(file = "./plots/panACF.pdf")
par(mfrow = c(4,2))
acf(pan.convergence$sigma[1,1,], lag.max = iter,
    xlab = "lag", main = "Sigma")
acf(pan.convergence$psi[1,1,], lag.max = iter,
    xlab = "lag", main = "Psi (Intercept)")
acf(pan.convergence$psi[2,1,], lag.max = iter,
    xlab = "lag", main = "Psi (Slope)")
for (i in 1:5){
  acf(pan.convergence$beta[i,1,], lag.max = iter,
      xlab = "lag", main = paste("Beta",i-1))
}
dev.off()

source("./functions/mi.pan.R")
set.seed(123456)
imp.pan <- mi.pan(group = dat.long$ID,
                  y = matrix(dat.long$DELINQ, ncol = 1),
                  x = dat.long[,c("TIME","FEMALE","RE","GY")],
                  xcol = 1:5, zcol = 1:2, m = 100, paniter = 5000)

# Chapter 5.7 Repeated data analysis and pooling of results in R

# Chapter 5.7.1 t-test
load("./Rworkspace/pmm_mt.Rdata")
com <-  complete(imp.mice.mt, 1)
t.test(ACRIM~FEMALE,com)

ttest <- lm(ACRIM~FEMALE, com)
s <- summary(ttest)
s$coefficients

est <- se <- vector(length = imp.mice.mt$m, mode = "list")

for (i in 1:imp.mice.mt$m){
  com <- complete(imp.mice.mt, i)
  ttest <- lm(ACRIM~FEMALE, com)
  s <- summary(ttest)
  est[[i]] <- s$coefficients[2,1]
  se[[i]] <- s$coefficients[2,2]
}
library(norm2)
miinf <- miInference(est, se)
print(miinf)

miinf$est + c(-1,1) * qt(.975,miinf$df) * miinf$std.err

mira <- with(imp.mice.mt, lm(ACRIM~FEMALE))
result<-summary(pool(mira))
print(result[,1:5])
print(result[,6:10])

# Chapter 5.7.2 Linear model

est <- se <- vector(length=imp.mice.mt$m, mode = "list")

for (i in 1:imp.mice.mt$m){
  com <- complete(imp.mice.mt, i)
  model <- lm(ACRIM~FEMALE+RE+GY, com)
  s <- summary(model)
  est[[i]] <- s$coefficients[,1]
  se[[i]] <- s$coefficients[,2]
}

library(norm2)
miinf <- miInference(est, se)
print(miinf)

lower <- miinf$est -1 * qt(.975, miinf$df) * miinf$std.err
print(lower)
upper <- miinf$est +1 * qt(.975, miinf$df) * miinf$std.err
print(upper)

mira <- with(imp.mice.mt, lm(ACRIM~FEMALE+RE+GY))
result<-round(summary(pool(mira)),3)
print(result[,1:5])
print(result[,6:10])

# Chapter 5.7.3 Growth model

model <- '
  i =~ 1*ACRIM + 1*BCRIM + 1*CCRIM + 1*DCRIM
s =~ 0*ACRIM + 1*BCRIM + 2*CCRIM + 3*DCRIM
q =~ 0*ACRIM + 1*BCRIM + 4*CCRIM + 9*DCRIM
i ~~ 0*q
s ~~ 0*q
q ~~ 0*q
'

est <-  se <- vector(length=100, mode="list")
options(warn=2)
fail <- NULL
for (m in 1:100){
  fit <- try(growth(model, data = complete(imp.mice.mt,m)))
  if (any(class(fit)=="try-error")){fail=c(fail,m); next}
  s <- parameterEstimates(fit)
  est[[m]] <- s$est
  se[[m]] <- s$se
  nam <- s[,1:3]
}
print(fail)

miinf <- miInference(est, se)

output <- data.frame(nam, Est = miinf$est, SE = miinf$std.err,
                     df = miinf$df, p = miinf$p, Pct.mis = miinf$mis.inf)
output[,-3:-1] <- round(output[,-3:-1],3)
output[c(27:29,20:22),]

# Chapter 5.7.4 Generalized linear mixed-effects model

options(warn=1)

cor2fisher <- function(r) {(1/2*log( ( 1 + r) / ( 1 - r ) ))  }
fisher2cor <- function(z) {( exp(2*z) - 1 )/ ( exp(2*z) + 1 ) }

est <- se <- rvar <- rcor<- vector(length = 100,
                                   mode = "list")
for (i in 1:100){
  cat(i,"\n")
  com <- complete(imp.mice.mt, i)
  com.l <- reshape(data=com,
                   varying = list(c("ACRIM","BCRIM","CCRIM","DCRIM")),
                   direction = "long", v.names = "DELINQ",
                   times = 0:3, timevar = "TIME", idvar = "ID")
  com.l$ID <- as.factor(com.l$ID)
  poissonmodel <- glmmTMB::glmmTMB(
    DELINQ ~ (TIME+I(TIME^2)) + (TIME+I(TIME^2)|ID),
    data = com.l, family = poisson)
  s <- summary(poissonmodel)
  est[[i]] <- s$coefficients$cond[,1]
  se[[i]] <- s$coefficients$cond[,2]
  rvar[[i]] <- attr(
    glmmTMB::VarCorr(poissonmodel)[[1]]$"ID","stddev"
  )^2
  rcor[[i]] <- cor2fisher(attr(
    glmmTMB::VarCorr(poissonmodel)[[1]]$"ID",
    "correlation"
  ))
}

miInference(est,se)

exp(-1.57080)

library(abind)
apply(abind(rvar, along = 2), 1, mean)

as.dist(round(fisher2cor(apply(abind(rcor, along = 3), c(1, 2), mean)), 3))

# Chapter 5.7.5 Two-level negative binomial hurdle model

dir.create("./Routput/mira/")
library(glmmTMB)
options(warn=1)
for (i in 1:100)
{
  cat(i,"\n")
  com <- complete(imp.mice.mt, i)
  com.l <- reshape(data=com,
                   varying = list(
                     c("ACRIM","BCRIM","CCRIM","DCRIM")),
                   direction = "long",v.names = "DELINQ",
                   times = 0:3, timevar = "TIME", idvar = "ID")
  com.l$ID <- as.factor(com.l$ID)
  twolevel.hnbmodel <- try(
    glmmTMB(DELINQ ~ (TIME+I(TIME^2)+FEMALE+GY) + (1|ID),
            zi = ~(TIME+FEMALE+GY),
            data = com.l,
            family = list(family = "truncated_nbinom2",
                          link = "log")))
  save(twolevel.hnbmodel, file = paste0("./Routput/mira/twolevel.hnbmodel",i))
}

est.cond <- est.zi <- se.cond <- se.zi <- theta <- rvar <-
  vector(length = 100, mode = "list")
for (i in 1:100){
  load(paste0("./Routput/mira/twolevel.hnbmodel",i))
  s <- summary(twolevel.hnbmodel)
  
  est.cond[[i]] <- s$coefficients$cond[,1]
  est.zi[[i]] <- s$coefficients$zi[,1]
  
  se.cond[[i]] <- s$coefficients$cond[,2]
  se.zi[[i]] <- s$coefficients$zi[,2]
  
  theta[[i]] <- s$sigma
  
  rvar[[i]] <- attr(VarCorr(twolevel.hnbmodel)[[1]]$"ID","stddev")^2
}

miInference(est.cond,se.cond)

miInference(est.zi,se.zi)

mean(unlist(theta))

library(abind)
randomvar <- apply(abind(rvar, along = 2),1,mean)
randomvar

# Chapter 5.7.6 Correlation coefficients

library(miceadds)
correlations <- micombine.cor(imp.mice.mt,
                              variables = 1:4,
                              method = "pearson")
print(correlations[,1:6], digits=3)
print(correlations[,c(1,2,7:11)], digits=3)

print(attr(correlations, "r_matrix"), digits = 3)

# Chapter 5.7.7 Further useful tools for multiple imputation inference

# Chapter 5.7.7.1 Pooling explained variance in a linear model

pool.r.squared(with(imp.mice.mt,
                    lm(ACRIM~FEMALE+RE+GY)))

pool.r.squared(with(imp.mice.mt,
                    lm(ACRIM~FEMALE+RE+GY)),adjusted=TRUE)

# Chapter 5.7.7.2 Comparison of two nested models fitted to imputed data

fit1 <- with(imp.mice.mt, lm(ACRIM~FEMALE))
fit2 <- with(imp.mice.mt, lm(ACRIM~FEMALE+RE))
fit3 <- with(imp.mice.mt, lm(ACRIM~FEMALE+RE+GY))

comp1 <-  pool.compare(fit2, fit1, method= "Wald")
print(comp1$pvalue, digits=3)

comp2 <- pool.compare(fit3, fit2, method= "Wald")
print(comp2$pvalue, digits=2)

# Chapter 5.7.7.3 Combination of F-values

fvalues = vector(length = imp.mice.mt$m)

for (i in 1:imp.mice.mt$m){
  com <- complete(imp.mice.mt,i)
  model <- lm(ACRIM~FEMALE+RE+GY, com)
  s <- summary(model)
  fvalues[i] <- s$fstatistic["value"]
}
micombine.F(fvalues, 3)

# Chapter 5.7.7.4 Analysis of variance

crim4a <- read.csv2("./data/crim4a.csv",row.names = 1)
set.seed(123456)
imp.mice.pmm.crim4a <- mice(crim4a, m = 100, maxit = 20, printFlag = FALSE)
save(imp.mice.pmm.crim4a,file="./imputations/pmm/imp.mice.pmm.crim4a.Rdata")

result<-mi.anova(imp.mice.pmm.crim4a, "ACRIM~FEMALE*SCHOOL")

round(result$r.squared, 3)
colnames(result$anova.table)[7]<-"p.eta2"
round(result$anova.table, 3)

fvalues<- vector(length=imp.mice.pmm.crim4a$m, mode="list")

for (i in 1:imp.mice.pmm.crim4a$m){
  com <- complete(imp.mice.pmm.crim4a,i)
  model <- lm(ACRIM~FEMALE*SCHOOL, com)
  a <- anova(model)
  fvalues[[i]] <- a$`F value`
}

library(abind)
fvalues <- abind(fvalues, along = 2)
micombine.F(fvalues[1,],1)
micombine.F(fvalues[2,],2)
micombine.F(fvalues[3,],2)

# Chapter 5.7.7.5 Combination of χ2-values and other model fit statistics

model <- '
 i =~ 1*ACRIM + 1*BCRIM + 1*CCRIM + 1*DCRIM
s =~ 0*ACRIM + 1*BCRIM + 2*CCRIM + 3*DCRIM
q =~ 0*ACRIM + 1*BCRIM + 4*CCRIM + 9*DCRIM
i ~~ 0*q
s ~~ 0*q
q ~~ 0*q
'

library(semTools)
implist <- vector(length=100, mode="list")
for (i in 1:100){implist[[i]]<-complete(imp.mice.mt,i)}
growth.res <- growth.mi(model, implist)
summary(growth.res, fit.measures=TRUE, fmi=TRUE, ci=TRUE, asymptotic=TRUE)
fitmeasures(growth.res,"srmr_mplus")

# Chapter 5.7.7.6 datlist2mids

implist <- vector(length = 100, mode = "list")
for (i in 1:100){
  implist[[i]] <- read.table(paste0("./imputations/other/norm2/imp",i,".dat"),
                             col.names =
                               c("id", "FEMALE", "RE", "GY",
                                 "ACRIM","BCRIM","CCRIM","DCRIM"))
}
implist.mids <- miceadds::datlist2mids(implist)

mira <- with(implist.mids,
             glm(ACRIM~FEMALE+RE+GY, family=poisson))
result <- summary(pool(mira))
round(result[,1:5], 3)
round(result[,6:10], 3)

# Chapter 5.8 Diagnostic Tools

load("./data/dat3.Rdata")
ini <- mice(dat3, maxit = 0)
method <- ini$method
method[1:8] <- "norm"
set.seed(123456)
imp.mice.norm = mice(dat3, method = method, m = 100, maxit = 20, printFlag = FALSE)
save(imp.mice.norm,file="./imputations/other/micenorm/imp.mice.norm.Rdata")

imp.acrim.c <- data.frame(imp.mice.norm$imp$acrim.c)
head(imp.acrim.c[,1:5])

imp.acrim.c.long <- reshape(imp.acrim.c,
                            varying = list(
                              c(paste0("X", 1:100))),
                            direction = "long", v.names = "acrim.c",
                            times = 1:100, timevar = "IMP", idvar = "ID")

sub <- subset(imp.acrim.c.long,
              imp.acrim.c.long$IMP <= 20)
pdf("./plots/plausibility1.pdf")
boxplot(acrim.c~IMP,sub,
        xlab = "Imputation number", ylab = "acrim.c",
        pch=20, col="gray80")
dev.off()


load("./Rworkspace/pmm_mt.Rdata")
imp.acrim <- data.frame(imp.mice.mt$imp$ACRIM)
head(imp.acrim[,1:15])

imp.acrim.long <- reshape(imp.acrim,
                          varying=list(
                            c(paste0("X",1:100))),
                          direction="long",v.names="ACRIM",
                          times = 1:100, timevar = "IMP", idvar="ID")
sub.com = subset(dat, !is.na(dat$ACRIM))
pdf("./plots/plausibility2.pdf")
hist(sub.com$ACRIM,
     col= gray(.1, alpha=.5), freq=FALSE,
     xlab="Distribution of ACRIM", main="")
hist(imp.acrim.long$ACRIM, col= gray(.8, alpha=.5),
     freq=FALSE, add=TRUE, breaks=0:12)
legend("topright",
       legend = c("observed","imputed"),
       fill=c( gray(.1, alpha=.5), gray(.8, alpha=.5)))
dev.off()

imp.mice.mt.copy <- imp.mice.mt
imp.mice.mt.copy$m <- 5
library(lattice)
pdf("./plots/plausibility3.pdf")
stripplot(imp.mice.mt.copy,
          data = ACRIM~.imp,
          col = c("black", gray(.6, alpha=.5)),
          pch = 20)
dev.off() 

source("./functions/compare.percent.count.R")
result<-compare.percent.count(orig=dat,
                              imp = imp.mice.mt,
                              counts = c("ACRIM","BCRIM"))
result$ACRIM[,1:11]
result$BCRIM[,1:11]

source("./functions/compare.obs.imp.R")
result <-compare.obs.imp(dat,imp.mice.mt)
print(result,digits=3)

# Chapter 5.8.6 Assessing the suitability of the imputation method by means of Monte Carlo simulation

cc <- dat[complete.cases(dat),-1]
cdata <- vector(length = 1000, mode = "list")
set.seed(123456)
for (i in 1:1000){
  s <- sample(1:nrow(cc), nrow(dat), replace = TRUE)
  cdata[[i]] <- cc[s,]
}

set.seed(123456)
mdata <- cdata
for (i in 1:1000){
  s <- sample(1:nrow(dat), round(.33*nrow(dat)),
              replace = FALSE)
  mdata[[i]][s,"ACRIM"] <- NA
}

cd.est <- vector(length = 1000, mode = "list")
for (i in 1:1000){
  fit <- growth(model, data = cdata[[i]])
  s <- parameterEstimates(fit)
  cd.est[[i]] <- s$est
}

library(abind)
cd.est.matrix <- abind(cd.est, along = 2)
cd.est.matrix <- cd.est.matrix[c(20,21,22,27,28,29),]
Q = apply (cd.est.matrix, 1, mean)
names(Q) = c("Ivar","Svar","IS","I", "S","Q" )

imp <- vector(length = 1000, mode = "list")
library(mice)
set.seed(123456)
for (i in 1:1000){imp[[i]] <- mice(mdata[[i]], printFlag = FALSE)}

library(norm2)
mi.est <- vector(length = 1000, mode = "list")
for(i in 1:1000)
{
  est <- se <- vector(length = imp[[i]]$m, mode = "list")
  for (m in 1:imp[[i]]$m){
    fit <- growth(model, data = complete(imp[[i]],m))
    s <- parameterEstimates(fit)
    est[[m]] <- s$est
    se[[m]] <- s$se
  }
  miinf <- (miInference(est, se))
  mi.est[[i]] <- miinf$est
}

mi.est.matrix <- abind(mi.est, along=2)
mi.est.matrix <- mi.est.matrix[c(20,21,22,27,28,29),]
Qhat <- apply (mi.est.matrix, 1, mean)
names(Qhat) <- c("Ivar","Svar","IS","I", "S","Q" )

Bias <- Q-Qhat
round(Bias, 3)
PBias = abs(((Q-Qhat)/Q)*100)
round(PBias, 3)

dir.create("./montecarlo/chapter5",recursive=TRUE)
save(cdata, file="./montecarlo/chapter5/cdata.Rdata")
save(mdata, file="./montecarlo/chapter5/mdata.Rdata")
save(imp, file="./montecarlo/chapter5/imp.Rdata")
save(cd.est, file="./montecarlo/chapter5/cd.est.Rdata")
save(mi.est, file="./montecarlo/chapter5/mi.est.Rdata")

