# Applied Multiple Imputation
# Advantages, Pitfalls, New Developments and Applications in R
#
# Kristian Kleinke, Jost Reinecke, Daniel Salfran, Martin Spiess
################################################################

# Chapter 6
#######################

#################################################################
# Note: we use relative file paths                              #
# For the examples to work, please set your R working directory #
# to the parent folder: "AppliedMI"                             #
#################################################################

# Chapter 6.2.2 The countimp package

devtools::install_git(url =
        "https://github.com/kkleinke/countimp", branch = "master")

library(countimp)
help("mice.impute.poisson")

data(crim4w)
data(crim4l)

ini <- mice(crim4w, maxit=0) 
meth <- ini$method
pred <- ini$predictorMatrix

meth[6:9] <- "poisson"

pred["ACRIM",c("CCRIM","DCRIM")] <-0 
pred["BCRIM",c("ACRIM","DCRIM")] <-0 
pred["CCRIM",c("BCRIM","DCRIM")] <-0 
pred["DCRIM",c("BCRIM","CCRIM")] <-0

pred[,"id"] <- 0
pred[,c("HA")] <- 0

set.seed(123456)
imp.poisson <- countimp(data = crim4w, method = meth, m = 100,
               predictorMatrix = pred, maxit = 20, print = FALSE)

set.seed(123456)
meth[6:9] <- "quasipoisson"
imp.quasipoisson <- countimp(data = crim4w, method = meth, m =
                    100, predictorMatrix = pred, maxit = 20, 
                    print = FALSE)

set.seed(123456)
meth[6:9] <- "quasipoisson.boot"
imp.qp.boot <- countimp(data = crim4w, method = meth, m = 100,
               predictorMatrix = pred, maxit = 20, print = FALSE)

set.seed(123456)
meth[6:9] <- "hnb"
imp.hnb <- countimp(data = crim4w, method = meth, m = 100,
           predictorMatrix = pred, maxit = 20, print = FALSE)

pred[c("ACRIM", "BCRIM", "CCRIM", "DCRIM"), c("FEMALE", "RE","GY")] <- 2 
set.seed(123456)
imp.hnb2 <- countimp(data = crim4w, method = meth, m = 100, 
                     predictorMatrix = pred, maxit = 20, print = FALSE)

pred[c("ACRIM", "BCRIM", "CCRIM", "DCRIM"), c("FEMALE", "RE","GY")] <- 3 
set.seed(123456)
imp.hnb3 <- countimp(data = crim4w, method = meth, m = 100, 
                     predictorMatrix = pred, maxit = 20, print = FALSE)

data(crim4l)
crim4l$TIME2 <- crim4l$TIMEˆ2
ini <- mice(crim4l, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix
meth["DELINQ"] <- "2l.poisson"
pred["DELINQ", ] <- c(1, 1, 1, 2, 0, -2, 1)
options(warn = 1)

set.seed(123456)
imp.2l.poisson <- countimp(data = crim4l, m = 100, method =
                          meth, predictorMatrix = pred, maxit = 20)

meth["DELINQ"] <- "2l.nb"
set.seed(123456)
imp.2l.nb <- countimp(data = crim4l, method = meth, m = 100,
                      predictorMatrix = pred, maxit = 20)

meth["DELINQ"] <- "2l.hnb"
set.seed(123456)
imp.2l.hnb <- countimp(data = crim4l, method = meth, m = 100, 
                       predictorMatrix = pred, maxit = 20)

# Chapter 6.3 ImputeRobust – Imputation based on Generalized Additive Models for Location, Scale and Shape

devtools::install_git(url =
          "https://github.com/dsalfran/ImputeRobust", branch = "master")

ini <- mice(crim4w, maxit=0) 
meth <- ini$method
pred <- ini$predictorMatrix

meth[6:9] <- "gamlss"

pred["ACRIM",c("CCRIM","DCRIM")] <-0 
pred["BCRIM",c("ACRIM","DCRIM")] <-0 
pred["CCRIM",c("BCRIM","DCRIM")] <-0
pred["DCRIM",c("BCRIM","CCRIM")] <-0

pred[,"id"] <- 0
pred[,c("HA")] <- 0

set.seed(123456)
library(ImputeRobust)
imp.gamlss <- mice(data = crim4w, method = meth, m = 100, 
                     predictorMatrix = pred, maxit = 20, print = FALSE)

imp.gamlss.1par <- mice( data = crim4w, method = meth, m = 100, 
                         predictorMatrix = pred, maxit = 20, 
                         n.ind.par = 1, print = FALSE)

meth["DELINQ"] <- "gamlssPO"
imp.gamlssPO <- mice( data = crim4w, method = meth, m = 100,
                      predictorMatrix = pred, maxit = 20, print = FALSE)

meth["DELINQ"] <- "gamlssZIP"
imp.gamlssZIP <- mice( data = crim4w, method = meth, m = 100,
                       predictorMatrix = pred, maxit = 20, print = FALSE)


# Chapter 6.4 Quantile regression based multiple imputation

install.packages("Qtools")
library(Qtools)
meth["DELINQ"] <- "rq"
post=ini$post
post["DELINQ"] <- "imp[[j]][, i] <- round(squeeze(imp[[j]][,imp.rq i], c(0, 16)))"
imp.rq <- mice( data = crim4l, method = meth, m = 100, 
                predictorMatrix = pred, maxit = 20, post = post, print = FALSE)

# Chapter 6.5 Two-level predictive mean matching

library(miceadds)
search()

data(crim4l)
crim4l$TIME2 <- crim4l$TIMEˆ2
ini <- mice(crim4l, maxit=0)
meth <- ini$method
meth["DELINQ"] <- "2l.pmm"
pred <- ini$predictorMatrix
pred["DELINQ",] <- c(1,1,1,2,0,-2,1)
imp.2l.pmm <- mice( data = crim4l, method = meth, m = 100,
                  predictorMatrix = pred, maxit = 20, print = FALSE)


# additional R code

load("./imputations/chapter_6_comparison/imp.2l.poisson")
load("./imputations/chapter_6_comparison/imp.2l.nb")
load("./imputations/chapter_6_comparison/imp.2l.zip")
load("./imputations/chapter_6_comparison/imp.2l.zinb")
load("./imputations/chapter_6_comparison/imp.2l.hnb")
load("./imputations/chapter_6_comparison/imp.2l.pmm")

compare.percent.count<-
  function(orig,imp,counts){
    nam=names(orig[colSums(is.na(orig))>0])
    nam = nam[nam %in% counts]
    output <- vector(length(nam), mode="list")
    names(output)<-nam
    for (i in 1:length(nam)){
      data <- mice::complete(imp,"long")
      ry= is.na(orig[,nam[i]])
      ry[ry==TRUE]<-"imputed"
      ry[ry==FALSE]<-"observed"
      data <- cbind(data,ry)
      colnames(data)[length(colnames(data))]<-paste0("R.",nam[i])
      tab <- table(data[,nam[i]],data[,paste0("R.",nam[i])])
      
      
      output[[i]]<- rbind("imputed"=tab[,1]/sum(tab[,1])*100,
                          "observed"=tab[,2]/sum(tab[,2])*100)
      
    }
    return(output)
  }


result1 <-
  rbind(
    compare.obs.imp(crim4l,imp.2l.poisson)$DELINQ[2,],
    compare.obs.imp(crim4l,imp.2l.poisson)$DELINQ[1,],
    compare.obs.imp(crim4l,imp.2l.nb)$DELINQ[1,],
    compare.obs.imp(crim4l,imp.2l.zip)$DELINQ[1,],
    compare.obs.imp(crim4l,imp.2l.zinb)$DELINQ[1,],
    compare.obs.imp(crim4l,imp.2l.hnb)$DELINQ[1,],
    compare.obs.imp(crim4l,imp.2l.pmm)$DELINQ[1,]
  )

result2 <-
  rbind(
    compare.percent.count(crim4l,imp.2l.poisson,"DELINQ")$DELINQ[2,],
    compare.percent.count(crim4l,imp.2l.poisson,"DELINQ")$DELINQ[1,],
    compare.percent.count(crim4l,imp.2l.nb,"DELINQ")$DELINQ[1,],
    compare.percent.count(crim4l,imp.2l.zip,"DELINQ")$DELINQ[1,],
    compare.percent.count(crim4l,imp.2l.zinb,"DELINQ")$DELINQ[1,],
    compare.percent.count(crim4l,imp.2l.hnb,"DELINQ")$DELINQ[1,],
    compare.percent.count(crim4l,imp.2l.pmm,"DELINQ")$DELINQ[1,])

row.names(result2)<-row.names(result1)<-row.names(result1)<-c("observed", "2l.poisson","2l.nb","2l.zip","2l.zinb","2l.hnb","2l.pmm")

result<-cbind(result1[,-1],result2[,1:6])
library(Hmisc)
round(result,2) # Table 6.4

# Mplus Automation

impname=c("imp.2l.poisson","imp.2l.nb","imp.2l.zip","imp.2l.zinb","imp.2l.hnb","imp.2l.pmm")
mimeth=c("2l.poisson","2l.nb","2l.zip","2l.zinb","2l.hnb","2l.pmm")
for (m in 1:6){
  dir.create(paste0("./imputations/chapter_6_comparison/",mimeth[m]))}

for (m in 1:6){
  for (i in 1:100){
    com=complete(eval(parse(text = impname[m])),i)
    com=reshape(com[,-7],direction="wide", v.names="DELINQ", timevar="TIME", idvar="ID", varying = list(c("ACRIM","BCRIM","CCRIM","DCRIM")))
    out=cbind(1:nrow(com), 
              com[,c("FEMALE", "RE", "GY", "ACRIM", "BCRIM", "CCRIM", "DCRIM")])
    
    write.table( out , paste("./imputations/chapter_6_comparison/",mimeth[m], "/imp" , i , ".dat" , sep="") ,
                 quote=F , row.names=F , col.names= F)
  }
}


for (m in 1:6){
  for (i in 1:100){  
    cat(paste( "imp" , i , ".dat" , sep=""), "\n",
        file=paste0("./imputations/chapter_6_comparison/",mimeth[m],"/implist.dat"), append=TRUE)
  }}


library(MplusAutomation)
createModels("./imputations/chapter_6_comparison/template.txt")

oldwd=getwd()
setwd("./imputations/chapter_6_comparison/")
runModels(recursive=TRUE) # only works with an Mplus licence / when Mplus is installed!
setwd(oldwd)

allModelParameters <- readModels("./imputations/chapter_6_comparison", recursive = TRUE)

# model fit statistics

ModelFitSummary <- sapply(allModelParameters, "[", "summaries") 

modelfitresults<- do.call(rbind.data.frame,lapply(ModelFitSummary, subset, select=c(11,12,14,15,17,18,20,21)))
modelfitresults<-modelfitresults[c(4,2,6,5,1,3),] 
rownames(modelfitresults)<-mimeth

round(modelfitresults,2) # Table 6.5


### parameter estimates

unstandardizedOnly <- sapply(allModelParameters, "[", "parameters")
unstandardizedOnly <- sapply(unstandardizedOnly, "[", "unstandardized")

estimates <- data.frame(unstandardizedOnly[[1]][,1:2],
                        do.call(cbind.data.frame,lapply(unstandardizedOnly, subset, select=c(est,se)))
)

results=estimates
results <- results[25:48,]
results <- results[c(19:24,1:18),]
library(stringr)
results$paramHeader<-str_replace_all(results$paramHeader, ".ON", " ON ")
results$paramHeader<-str_replace_all(results$paramHeader, "Intercepts", "")
results$parameter <- paste(results$paramHeader,results$param,sep="")
row.names(results)<-results$parameter
ncol(results)
nrow(results)
head(results)
results=results[,c(-1,-2,-15)]
colnames(results)=c(rep(c("Est.","S.E."),6))

# order: 
# 1) hnb (est, se), 
# 2) nb (est, se), 
# 3) pmm (est, se), 
# 4) poisson (est, se)
# 5) zinb (est, se)
# 6 zip (est, se)

round(results,2) # Table 6.6
