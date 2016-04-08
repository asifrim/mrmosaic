## function to calculate fisher sum
## p: vector of p-values
fisher.sum <- function(p, zero.sub=0.00001, na.rm=FALSE){
  if(any(p>1, na.rm=TRUE)||any(p<0, na.rm=TRUE))
    stop("You provided bad p-values")
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  p[p==0] <- zero.sub
  if(na.rm)
    p<- p[!is.na(p)]
  S= -2*sum(log(p))
  res <- data.frame(S=S, num.p=length(p))
  return(res)
}
