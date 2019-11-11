prepare.mplus <- function(imp.data, dir, order=c("FEMALE", "RE", "GY", "ACRIM", "BCRIM", "CCRIM", "DCRIM")){
  
  for (m in 1:imp.data$m){
    
    com=complete(imp.data,m) 
    
    out=cbind(1:nrow(com), 
              com[,order])
    
    write.table( out , paste(dir, "/imp" , m , ".dat" , sep="") ,
                 quote=F , row.names=F , col.names= F)
  }
  
  
  for (m in 1:imp.data$m)
  {
    cat(paste( "imp" , m , ".dat" , sep=""), "\n",
        file=paste0(dir,"/implist.dat"),append=TRUE)
  }
}