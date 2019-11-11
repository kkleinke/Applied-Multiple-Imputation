compare.obs.imp <- function(orig,imp){
  nam <- names(imp$nmis[imp$nmis>0])
  output <- vector(length(nam), mode="list")
  names(output)<-nam
  for (i in 1:length(nam)){
    data <- mice::complete(imp,"long")
    ry <- is.na(orig[,nam[i]])
    ry[ry==TRUE] <- "imputed"
    ry[ry==FALSE] <- "observed"
    data <- cbind(data,ry)
    colnames(data)[length(colnames(data))] <- paste0("R.",nam[i])
    result<-aggregate(data[,nam[i]], list(data[,paste0("R.",nam[i])]),summary)
    result<-result[,2]
    rownames(result)=c("imputed","observed")
    output[[i]] <-result 
    }
return(output)
}