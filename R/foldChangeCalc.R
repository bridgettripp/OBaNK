linFold <- function(x) {
  library(gtools)
  x2 <- as.data.frame(t(x))
  x3 <- x2[order(x2[,1], decreasing = FALSE), ]
  x <- as.data.frame(t(x3))
  grps <- as.factor(x[1,])
  g1 <- sum(grps == levels(grps)[1])
  g2 <- sum(grps == levels(grps)[2])
  x <- as.data.frame(t(x[-1,]))#
  dcol <- colnames(x)
  fc <- data.frame()
  d1 <- type.convert(as.data.frame(t(x)))
  for (i in 1:dim(d1)[1]){ #features x foldchange
    z <- 1
    for (a in 1:g1){
      for (c in (g1+1):dim(d1)[2]){
        fc[i, z] <- foldchange(d1[i,a],d1[i,c])
        z <- z+1
      }}}
  fc = as.data.frame(t(fc)) #foldchange x feature
  colnames(fc) <- as.character(dcol)
  return(fc)}

logFold <- function(x){
  library(gtools)
  dcol <- colnames(x)
  fc <- data.frame()
  d1 <- type.convert(as.data.frame(t(x)))
  b <-  dim(d1)[2]/2
  for (i in 1:dim(d1)[1]){ #features x foldchange
  z <- 1
  for (a in 1:b){
    for (c in 1:b){
      fc[i,z] <- (d1[i, a+b]-d1[i,c])
      z <- z+1
    }}}
  fc = as.data.frame(t(fc)) #foldchange x feature
  colnames(fc) <- as.character(dcol)
  return(fc)}
