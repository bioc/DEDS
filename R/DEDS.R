###########################################################################
## Set Class
#require(methods)
#.initDEDS <- function() {
#setClass("DEDS", representation("list"), where=where)
#setClass("DEDS", representation("list"))
#}

############################
## using p values from different meausres
## X -- matrix of p values from different measures

deds.pval <- function(X, E=rep(0,ncol(X)), adj=c("fdr", "adjp"),
                          B=200) {
  ## draws unifom distribution along the columns of X
  genUniform <- function(X) { 
    X <- as.matrix(X)
    X2 <- apply(X, 2, function(z){
                runif(length(z), min=range(z)[1], max=range(z)[2])})
    X2
  }

  s <- svd(scale(X,scale=FALSE), nu=0)
  X.prime <- X %*% s$v
  bD <- euclidean(X, E)
  geneOrder <- order(bD)
  cat("Simulating null distribution ", B, " times, please wait...") 
  for (i in 1:B) {
    cat(i, "\t")
    if ((i%%5) == 0) cat("\n")
    Z.prime <- genUniform(X.prime)
    Z.prime <- Z.prime %*% t(s$v)
    bD <- cbind(bD, euclidean(Z.prime, E))
  }

  adj <- match.arg(adj)
  p <- switch(adj,
         fdr=deds.calcFDR(bD=bD[,-1], D=bD[,1], R=geneOrder-1),
         adjp=deds.calcAdjP(bD=bD[,-1], D=bD[,1], R=geneOrder-1))
  stats <- cbind(geneOrder=geneOrder, X[geneOrder,])
    
  res <- list(E=E, geneOrder=geneOrder, stats=stats, p=p, options=c("p","abs", "euclidean", adj))
  class(res) <- "DEDS"
  return(res)
  }          

deds.stat <- function(X, L, B=1000, testfun=list(t=comp.t(L), fc=comp.FC(L), sam=comp.SAM(L)),
                                tail=c("abs", "lower", "higher"),
                                distance=c("weuclid", "euclid"), adj=c("fdr", "adjp")) {
  ### sanity check
  tail <- match.arg(tail)
  distance <- match.arg(distance)
  adj <- match.arg(adj)
  newX <- deds.checkX(X, L, names(testfun), B)
  deds.checkothers(tail, distance, adj)
  options <- type2test(tail, distance, adj, L)

  ### initialization
  func.max=options$func.max
  adj.dist=options$adj.dist
  func.next.sample=options$func.next.sample
  func.compute.p=options$func.compute.p
  X <- newX$X
  L <- newX$L
  nT <- newX$nT
  B <- newX$B
  bD <- c()
  applyTest <- aggregateFun(testfun)
  
  ### print messages
  #cat("permutation of the dataset will be carried out... \n")

  ############### original data
  t.o <- applyTest(X)
  colnames(t.o) <- names(testfun)
  if (tail=="abs") t <- abs(t.o)
  E <- apply(t,2, func.max)
  cat("E of the original data is: ", E, "\n")
  if (adj.dist) {
    wval <- apply(t, 2, mad, na.rm=TRUE)
    wval <- 1/(wval^2)
  }
  else wval <- rep(1, nT)

  ################  start permutation 
  BT <- list(t)
  for (i in 1:B) {
    bX <- func.next.sample(X)
    tB <- applyTest(bX)
    
    if (tail=="abs") tB <- abs(tB)
    BT <- c(BT, list(tB))
    cat(i,'\t')
    if (i%%5==0) cat("\n")
  }

  ################## find E
  E <- apply(t(sapply(BT, function(z) apply(z, 2, func.max, na.rm=TRUE))), 2, func.max, na.rm=TRUE)
  cat("\nafter permutation, E is set at: (", E, ")\n")
  cat("distance to this extreme point will be measured ...\n")
  
  ################## compute distance to E
  geneOrder <- order(euclidean(t, E, wval))
  for (i in 1:(B+1)) 
    bD <- cbind(bD, euclidean(BT[[i]], E, wval))
  p <- func.compute.p(bD=bD[,-1], D=bD[,1], geneOrder-1)

  stats <- cbind(geneOrder=geneOrder, t.o[geneOrder,])
    
  res <- list(geneOrder=geneOrder, E=E, stats=stats, p=p, options=c("t", options))
  #Res <- new("DEDS", unclass(res))
  class(res) <- "DEDS"
  return(res)
}

deds.stat.linkC <- function(X, L, B=1000, tests=c("t", "fc", "sam"),
                                tail=c("abs", "lower", "higher"),
                                extras=NULL,
                                distance=c("weuclid", "euclid"), adj=c("fdr", "adjp"),
                                quick=TRUE) {
  tail <- match.arg(tail)
  distance <- match.arg(distance)
  adj <- match.arg(adj)
  newX <- deds.checkX(X, L, tests, B)
  deds.checkothers(tail, distance, adj)
  options <- c(tests, tail, distance, adj)
  if(quick) q<-1
  else q<-0
  if(is.null(extras)) extras <- deds.genExtra(newX$L, tests)
  
  res <- .C("get_deds_FDR",as.double(newX$X),as.integer(newX$nr),as.integer(newX$nc),as.integer(newX$L),as.character(options),as.single(extras), as.integer(q), as.integer(newX$nL), as.integer(newX$nT), as.integer(newX$B), E=double(newX$nT), R=integer(newX$nr), p=double(newX$nr),t=double(newX$nr*newX$nT), NAOK=TRUE, PACKAGE="DEDS")
  
  t.o <- matrix(res$t, byrow=FALSE, nc=newX$nT)
  colnames(t.o) <- tests

  geneOrder <- res$R+1
  p <- res$p
  stats <- cbind(geneOrder=geneOrder, t.o[geneOrder,])
   
  res <- list(E=res$E, geneOrder=geneOrder, stats=stats, p=p, options=c("t", options))
  class(res) <- "DEDS"
  return(res)
}


type2test <- function(tail=c("abs", "higher", "lower"),
                      distance=c("weuclid", "euclid"),
                      adj=c("fdr", "adjp"), L) {
  tail <- match.arg(tail)
  distance <- match.arg(distance)
  adj <- match.arg(adj)
  
  if(tail=="abs") func.max<-max
  if(tail=="higher") func.max<-max
  if(tail=="lower") func.max<-min

  if(distance=="weuclid") adj.dist<-TRUE
  if(distance=="euclid") adj.dist<-FALSE

  if(adj=="fdr") func.compute.p <- deds.calcFDR
  if(adj=="adjp") func.compute.p <- deds.calcAdjP

  L <- deds.checkclasslabel(L)
  func.next.sample=deds.next.sample(L)
  
  return(list(func.max=func.max, adj.dist=adj.dist, func.compute.p=func.compute.p, func.next.sample=func.next.sample))
}
  


aggregateFun <- function(...){
  funlist <- list(...)  #let the user supply a list
  if(length(funlist) == 1 && is.list(funlist[[1]]))  funlist <- funlist[[1]]
  if (length(funlist) == 0) stop("no tests available in input...")
  if (length(funlist) == 1) cat("warning: only one test selected, deds is better suited summarizing multiple tests, calculation will continue anyway...")
  function(X) {
    fval <- c()
    for (fun in funlist)   fval <- cbind(fval,fun(X))
    return(as.matrix(fval))
  }
}


### returns a column-permuted data matrix according to the class label
deds.next.sample <- function(L)  {
  nL <- length(unique(L))
  if (nL == 1) {
    function(X) {
      bL <- sample(c(-1,1), ncol(X), replace=TRUE)
      bX <- sweep(X, 2, bL, "*")
      return(bX)
    }
  }
  else {
    function(X) {
      bL <- sample(L)
      bX <- X[, order(bL)]
      return(bX)
    }
  }
}
        
# calculates FDR for class 'DEDS'
deds.calcFDR <- function(bD, D, R) {
  nc <- ncol(bD)
  nr <- nrow(bD)
  p <- .C("calc_FDR", as.single(bD), as.single(D), as.integer(R), as.integer(nr), as.integer(nc), p=single(nr), PACKAGE="DEDS")$p
  return(p)
}

# calculates adjusted p for class 'DEDS'
deds.calcAdjP <- function(bD, D, R) {
  nc <- ncol(bD)
  nr <- nrow(bD)
  p <- .C("calc_adjP", as.single(bD), as.single(D), as.integer(R), as.integer(nr), as.integer(nc), p=single(nr), PACKAGE="DEDS")$p
  return(p)
}


 
euclidean <- function(X, center, wval) {
  X <- if (is.vector(X))  matrix(X, ncol = length(X))
       else as.matrix(X)
  if (missing(wval)) wval <- rep(1, ncol(X))
  if (missing(center)) center <- colMeans(X, na.rm=TRUE)
  XC <- scale(X,center, scale=FALSE)
  res <- sweep(XC^2, 2, wval, "*")%*%rep(1,ncol(XC))
  return(sqrt(res))
}

comp.stat <- function(X, L, test=c("t","fc","sam","f","modt","modf","B"), extra=NULL){
  test <- match.arg(test)
  newX <- test.checkX(X, L, test)
  options <- test
  if(is.null(extra)) extra <- deds.genExtra(newX$L, test)
  res <- .C('get_stat',as.double(newX$X),as.integer(newX$nr),as.integer(newX$nc),
            as.integer(newX$L),t=single(newX$nr), as.character(options),
            as.single(extra), as.integer(newX$nL), NAOK=TRUE, PACKAGE="DEDS")$t
  return(res)
}


comp.unadjp <- function(X, L, B=1000, test=c("t","fc","sam","f"),tail=c("abs", "lower", "higher"), extra=NULL){
  tail <- match.arg(tail)
  test <- match.arg(test)
  newX <- deds.checkX(X, L, test, B)
  deds.checkothers(tail, distance="euclid", adj="fdr")
  options <- c(test, tail, "euclid", "fdr")
  if(is.null(extra)) extra <- deds.genExtra(newX$L, test)
  res <- .C('get_unadjp',as.double(newX$X),as.integer(newX$nr),as.integer(newX$nc),
           as.integer(newX$L),t=single(newX$nr), p=single(newX$nr),
            as.character(options), as.double(extra), as.integer(newX$nL), as.integer(newX$B), NAOK=TRUE, PACKAGE="DEDS")
  ret <- as.matrix(cbind(res$t,res$p))
  colnames(ret) <- c(test, "unadj.p")
  return(ret)
}

comp.adjp <- function(X, L, B=1000, test=c("t","fc","sam","f","modt","modf"),tail=c("abs", "lower", "higher"), extra=NULL){
  tail <- match.arg(tail)
  test <- match.arg(test)
  newX <- deds.checkX(X, L, test,  B)
  deds.checkothers(tail, distance="euclid", adj="fdr")
  options <- c(test, tail, "euclid", "fdr")
  if(is.null(extra)) extra <- deds.genExtra(newX$L, test)
  res <- .C('get_adjp',as.double(newX$X),as.integer(newX$nr),as.integer(newX$nc),
           as.integer(newX$L),t=single(newX$nr), p=single(newX$nr), adjp=single(newX$nr),
            r=integer(newX$nr), as.character(options), as.double(extra),
            as.integer(newX$nL), as.integer(newX$B), NAOK=TRUE, PACKAGE="DEDS")
  ret <- as.matrix(cbind(res$r+1,res$t,res$p,res$adjp))
  colnames(ret) <- c("order",test, "unadj.p", "adj.p")
  return(ret)
}

comp.fdr <- function(X, L, B=1000, test=c("t","fc","sam","f","modt","modf"),tail=c("abs", "lower", "higher"), extra=NULL){
  tail <- match.arg(tail)
  test <- match.arg(test)
  newX <- deds.checkX(X, L, test, B)
  deds.checkothers(tail, distance="euclid", adj="fdr")
  options <- c(test, tail, "euclid", "fdr")
  if(is.null(extra)) extra <- deds.genExtra(newX$L, test)
  res <- .C('get_fdr',as.double(newX$X),as.integer(newX$nr),as.integer(newX$nc),
           as.integer(newX$L),t=single(newX$nr), p=single(newX$nr), q=single(newX$nr),
            r=integer(newX$nr), as.character(options), as.double(extra), as.integer(newX$nL),
            as.integer(newX$B), NAOK=TRUE, PACKAGE="DEDS")
  ret <- as.matrix(cbind(res$r+1,res$t,res$p,res$q))
  colnames(ret) <- c("order", test, "unadj.p", "qvalues")
  return(ret)
}

deds.genExtra <- function(classlabel, tests) {
  ntests <- if(!missing(tests)) length(tests)
            else 0
  if(ntests==0) stop("please provide tests names in c('t','f','sam','fc','modt','modf','B')");
  nL <- length(unique(classlabel))
  extra <- c()
  for (i in 1:ntests) {
    test <- tests[i]
    if((!is.character(test))||(!is.vector(test))||(!any(test==c("t","f","sam","fc","modt","modf","B")))){
        mssg <- paste("your setting of test is",test,
                      "\nthe test needs to be a single character from c('t','f','sam','fc','modt','modf','B')")
        stop(mssg)
    }
    if(test=="t"||test=="modt"||test=="fc"||test=="f"||test=="modf")
        extra <- c(extra, nL)
    else if(test=="sam")
        extra <- c(extra, 0.5)
    else if(test=="B")
        extra <- c(extra, 0.01)
  }
  return(extra)
}
  
deds.checkX<-function(X, classlabel, tests, B){
  if((!is.matrix(X)) || !(is.numeric(X)))
     stop("X needs to be a matrix\n")
  if(ncol(X)!=length(classlabel))
    stop("the number of column of X needs to be the same as the length of classlabel\n")
  X <- X[, order(classlabel)]
  L <- deds.checkclasslabel(classlabel,tests)
  B <- deds.checkB(L, B)
  
  return(list(X=X, nr=nrow(X), nc=ncol(X), L=L, nL=length(unique(L)), nT=length(tests), B=B))
    
}


deds.checkothers<-function(tail="abs", distance="weuclid", adj="fdr")
{
  if((!is.character(tail))||(!is.vector(tail))||(length(tail)>1)
     ||(!any(tail==c("abs", "higher", "lower"))))
    stop(paste("the tail needs to be from c('abs', 'higher','lower')\n","your tail=",tail,"\n"))
  if((!is.character(distance))||(!is.vector(distance))||(length(distance)>1)
     ||(!any(distance==c("weuclid", "euclid"))))
    stop(paste("the distance argument needs to be from c('weuclid', 'euclid')\n","your distance=",distance,"\n"))
  if((!is.character(adj))||(!is.vector(adj))||(length(adj)>1)
     ||(!any(adj==c("fdr", "adjp"))))
    stop(paste("the adj argument needs to be from c('fdr', 'adjp')\n","your adj=",adj,"\n"))
}
 
deds.checkclasslabel <- function(classlabel, tests) {
  ntests <- if(!missing(tests)) length(tests)
            else 0
  newL <- as.integer(classlabel) 
  l <- length(unique(newL))
  if((!is.integer(newL)) ||(!is.vector(newL)))
     stop("classlabel needs to be just a vector of integers")
  l <- length(unique(newL))
  if(l == 1) newL <- rep(0, length(newL))
  else  newL <- rep(0:(l-1), table(newL))

  
  if (ntests>=1) {
    for (i in 1:ntests) {
      test <- tests[i]
      if((!is.character(test))||(!is.vector(test))||(!any(test==c("t","f","sam","fc","modt","modf","B","ebayes")))){
        mssg <- paste("your setting of test is",test,
                      "\nthe test needs to be a character string from c('t','f','sam','fc','modt','modf','B','ebayes')")
        stop(mssg)
      }
      if((test=="t")||test=="modt"||test=="B"||test=="ebayes"){
        K<-max(newL)
        if(K>2) {
           mssg <- paste("t test can not handel more than two groups\n",
                   "Your setting of classlabel is:\n",paste(newL, collapse=" "),"\n")
          stop(mssg)
         }
      }
      if((test=="f")||test=="modf"){
        K<-max(newL)
        if(K<1){
          mssg <- paste("in F test, we need at least two groups\n",
                   "Your setting of classlabel is:\n",paste(newL, collapse=" "),"\n")
          stop(mssg)
        }
        for(i in c(0:K))
          if(sum(newL==i)<2) {
            mssg <- paste("in F test, as the number of the groups is",K+1, "However, we found the number of objects in set with index",i, "has less than two objects\n", "the settings are:\n",paste(newL, collapse=" "),"\n")
            stop(mssg)
          }
      }
    }
  }
  return(newL)
}

#this functions finds the maximum number of permutation
#if the the initial B=0, or initial B greater than the maximum number of
#permutation maxB, it will return all possible of number of permutation.
deds.checkB<-function(classlabel, B, verbose=TRUE)
{
  if((length(B)>1) || !(is.integer(as.integer(B))) ||(!is.vector(B)))
     stop(paste("B needs to be just a integer\n","your B=",B,"\n"))
  if(B<0)
     stop(paste("the number of Permutations (B) needs to be positive\n, your B = ",B))
  
  n<-length(classlabel)
  k <- max(classlabel)+1
  nk <- table(classlabel)
   
  if (k > 1) {
    maxB <- 1
    rest <- n
    for (i in c(1:k)) {
      maxB <- maxB * choose(rest, nk[i])
      rest <- rest - nk[i]
    }
  }
  else maxB=2^n
   
  if((B>maxB) ||(B==0) && verbose ){
    cat("We'll do complete permutations, B = ", maxB, "\n")
    return(maxB)
  }
  else {
    cat("We'll do random permutations, B = ", B, "\n")
    return(B)
  }
}


deds.chooseTest <- function(L=NULL, tests=c("t","sam","fc"))  {
  if (is.null(L)) stop("need input of class label of data...")
  L <- deds.checkclasslabel(L, tests)
  nT <- length(tests)
  testfuncs <- c()
  for (i in 1:nT) {
    test <- tests[i]
    func<-switch(test,t=comp.t(L),
                      fc=comp.FC(L),
                      sam=comp.SAM(L),
                      f=comp.F(L),
                      modt=comp.modt(L),
                      modf=comp.modF(L),
                      B=comp.B(L))
    testfuncs <- c(testfuncs, list(func))
  }
  names(testfuncs) <- tests
  
  return(testfuncs)
}

topgenes <- function(obj,number=10,genelist=NULL, sort.by=c("deds", colnames(obj$stats[,-1])),
                     tail=c("abs", "lower", "higher")) {
#	Summary table of top genes
  if(data.class(obj)!="DEDS") stop("the argument obj must be of class DEDS ...")
  stats <- obj$stats[,-1]
    
  is.p <- obj$options[1]=="p"
  sort.by <- match.arg(sort.by, c("deds", colnames(stats)))
  tail <- match.arg(tail)
  
  if(sort.by=="deds") ord <- 1:number
  else if (!is.p) { #order by stats
    ord <- switch(tail,
                  abs=order(abs(stats[, sort.by]), decreasing=TRUE)[1:number],
                  lower=order(stats[, sort.by], decreasing=FALSE)[1:number],
                  higher=order(stats[, sort.by], decreasing=TRUE)[1:number])
  }
  else {# order by p
    ord <- order(stats[, sort.by], decreasing=FALSE)[1:number]
  }

  if(sort.by=="deds")
    statM <- cbind(geneOrder=obj$geneOrder[ord], DEDS=obj$p[ord], stats[ord,])
  else
    statM <- cbind(geneOrder=obj$geneOrder[ord], stats[ord, sort.by, drop=FALSE], DEDS=obj$p[ord], stats[ord, setdiff(colnames(stats), sort.by), drop=FALSE])
  
  if(is.null(genelist))
    tab <- data.frame(statM)
  else if(is.null(dim(genelist)))
    tab <- data.frame(Name=I((genelist[obj$geneOrder])[ord]),statM)
  else
    tab <- data.frame((genelist[obj$geneOrder])[ord,],statM)

  rownames(tab) <- as.character(1:number)
  tab
}



############# t test ###############

comp.t <- function(L=NULL, mu=0, var.equal=FALSE) {
  function(X) {
    if (is.vector(X)) X <- matrix(X, byrow=TRUE)
    newX <- test.checkX(X, L, "sam")
    L <- newX$L
    nL <- newX$nL
    nr <- newX$nr
    nc <- newX$nc
    if(nL==1) {
      df <- nc - 1
      stderr <- sqrt(apply(X,1,var,na.rm=TRUE)/nc)
      d <- rowMeans(X,na.rm=TRUE) - mu
      tstat <- d/stderr
      # pval <- 2 * pt(-abs(tstat), df)
      return(tstat)
    }
    else {
      G1 <- X[,L==0]
      G2 <- X[,L==1]
      d <- rowMeans(G1,na.rm=TRUE)-rowMeans(G2,na.rm=TRUE)
      v1 <- apply(G1,1,var,na.rm=TRUE)
      v2 <- apply(G2,1,var,na.rm=TRUE)
      n1 <- ncol(G1)
      n2 <- ncol(G2)
      if (var.equal) {
        df <- n1 + n2 - 2
        v <- ((n1 - 1) * v1 + (n2 - 1) * v2)/df
        stderr <- sqrt(v * (1/n1 + 1/n2))
      }
      else {
        stderr1 <- sqrt(v1/n1)
        stderr2 <- sqrt(v2/n2)
        stderr <- sqrt(stderr1^2 + stderr2^2)
        #df <- stderr^4/(stderr1^4/(n1 - 1) + stderr2^4/(n2 - 1))
      }
      tstat <- d/stderr
      # pval <- 2 * pt(-abs(tstat), df)
      return(tstat)
    }
  }
}


####### Fold Change ##################

comp.FC <- function(L=NULL, is.log=TRUE, FUN=mean) {
  function(X) {
    if (is.vector(X)) X <- matrix(X, byrow=TRUE)
    if(is.null(L)) L <- rep(0, ncol(X))
    newX <- test.checkX(X, L, "fc")
    L <- newX$L
    nL <- newX$nL
    nr <- newX$nr
    nc <- newX$nc

    if(nL==1) { # one-class
      fc <- apply(X, 1, FUN, na.rm=TRUE)
      return(fc)
    }
    else if(nL==2) { # two-class
      G1 <- X[, L==0]
      G2 <- X[, L==1]       
      m1 <- apply(G1,1,FUN,na.rm=TRUE)
      m2 <- apply(G2,1,FUN,na.rm=TRUE)
      if (is.log) fc <- m1-m2
      else fc <- m1/m2
      fc
    }
    else {  # multi-class
      m <- t(apply(X, 1, function(z) {tapply(z, L, FUN)}))
      if (is.log) 
         fc <-  apply(m, 1, function(z) {max(z) - min(z)})
      else
         fc <-  apply(m, 1, function(z) {max(z) /  min(z)})
      return(fc)
    }
  }
}

#########  SAM #########################

comp.SAM <- function(L=NULL,
                prob=0.5, B=200,
                stat.only=TRUE, verbose=FALSE, deltas,
                s.step=0.01, alpha.step=0.01, plot.it=FALSE) {
  function(X) {
    if(is.null(L)) L <- rep(0, ncol(X))
    newX <- test.checkX(X, L, "sam")
    L <- newX$L
    nL <- newX$nL
    
    sam.func <- if(nL==1) sam.oneclass.func
                else if(nL==2) sam.twoclass.func
                else sam.multiclass.func
    
    t <- sam.func(X,L,prob,B,stat.only,verbose,s.step=s.step,alpha.step=alpha.step,plot.it=plot.it)
    if(stat.only) return(t)
    else {
      n <- ncol(t)
      fdr.table <- sam.fdr(t[,2], t[,3:n], deltas)
      return(list(geneOrder=t[,1], sam=t[,2], fdr.table=fdr.table))
    }
  }
}

sam.oneclass.func <- function(X, L, prob=0.5, B=200, stat.only=TRUE, verbose=FALSE,
                              s.step=0.01, alpha.step=0.01, plot.it=FALSE) {
   n <- ncol(X)
   r <- rowMeans(X, na.rm=TRUE)
   s <- sqrt(rowSums(sweep(X, 1, r)^2, na.rm=TRUE) / (n*(n-1)))
   cat(prob)
   if(!is.null(prob))  s0 <- quantile(s,prob)
   else s0 <- sam.s0(r,s,s.step=s.step, alpha.step=alpha.step, plot.it=plot.it)
                
   if(verbose) cat("s0: ", s0, "\n")
   d <- r/(s+s0)
   if (stat.only) return(d)
   else {
     B <- deds.checkB(L, B)
     order.d <- order(d)
     sort.d <- sort(d)
     dB <- matrix(nr=nrow(X),nc=B)
     # start permutation
     cat("start permutation...\n")
     for (i in 1:B) {
       cat(i, "\t")
       if(i%%20 == 0) cat("\n")
       signs <- sample(c(1, -1), n, replace=TRUE)
       Xb <- sweep(X, 2, signs, "*")
       rb <- rowMeans(Xb,na.rm=TRUE)
       sb <- sqrt(sum(sweep(Xb, 1, rb)^2, na.rm=TRUE) / (n*(n-1)))
       db <- rb/(sb+s0)
       dB[,i] <- sort(db)
     }
     return(cbind(index=order.d, t=sort.d, dB))
   }
 }

sam.twoclass.func <- function(X, L, prob=0.5, B=200, stat.only=TRUE, verbose=FALSE,
                              s.step=0.01, alpha.step=0.01, plot.it=FALSE) {
  G1 <- X[, L==0]
  G2 <- X[, L==1]
  n1 <- ncol(G1)
  n2 <- ncol(G2)
  r <- rowMeans(G1,na.rm=TRUE) - rowMeans(G2,na.rm=TRUE)
  ss <- function(x) { sum((as.numeric(x)-mean(as.numeric(x),na.rm=TRUE))^2, na.rm=TRUE)}
  s<- sqrt((apply(G1,1,ss)+apply(G2,1,ss))*(1/n1+1/n2)/(n1+n2-2))
  if(!is.null(prob))   s0 <- quantile(s,prob)
  else s0 <- sam.s0(r,s,s.step=s.step, alpha.step=alpha.step, plot.it=plot.it)

  if(verbose) cat("s0: ", s0, "\n")
  d <- r/(s+s0)
  if (stat.only) return(d)
  else {
    order.d <- order(d)
    sort.d <- sort(d)
    B <- deds.checkB(L, B)
    dB <- matrix(nr=nrow(X),nc=B)
    # start permutation
    cat("start permutation...\n")
    for (i in 1:B) {
      cat(i, "\t")
      if(i%%20 == 0) cat("\n")
      id <- sample(c(1:(n1+n2)), n1+n2)
      G1 <- X[,id[1:n1]]
      G2 <- X[,id[(n1+1):(n1+n2)]]
      rb <- rowMeans(G2,na.rm=TRUE) - rowMeans(G1,na.rm=TRUE)
      sb<- sqrt((apply(G1,1,ss)+apply(G2,1,ss))*(1/n1+1/n2)/(n1+n2-2))
      db <- rb/(sb+s0)
      dB[,i] <- sort(db)
    }
    return(cbind(index=order.d, t=sort.d, dB))
  }
}

sam.multiclass.func <- function(X, L, prob=0.5, B=200, stat.only=TRUE, verbose=FALSE,
                                s.step=0.01, alpha.step=0.01, plot.it=FALSE) {
  n <- ncol(X)
  nks <- table(L)
  X.bar <- rowMeans(X, na.rm=TRUE)
  X.bar.L <- t(apply(X, 1, function(z) {tapply(z, L, mean)}))
  r <- sweep(X.bar.L, 1, X.bar)
  r <- rowSums(sweep(r^2, 2, nks, "*"), na.rm=TRUE)
  r <- sqrt((sum(nks) / prod(nks)) * r)
  s <- X - apply(X.bar.L,1, function(z) {rep(z, nks)})
  s <- sqrt((1/sum(nks-1)) * sum(1/nks) * sum(s^2, na.rm=TRUE))

  if(!is.null(prob))	s0 <- quantile(s,prob)
  else s0 <- sam.s0(r,s,s.step=s.step, alpha.step=alpha.step, plot.it=plot.it)
                
  if(verbose) cat("s0: ", s0, "\n")
  d <- r/(s+s0)
  if (stat.only) return(d)
  else {
    order.d <- order(d)
    sort.d <- sort(d)
    B <- deds.checkB(L, B)
    dB <- matrix(nr=nrow(X),nc=B)
     # start permutation
    cat("start permutation...\n")
    for (i in 1:B) {
           cat(i,"\t")
           if(i %% 20 ==0) cat("\n")
           L.b <- sample(L)
           X.b <- X[, order(L.b)]
           X.bar.b <- rowMeans(X.b, na.rm=TRUE)
           X.bar.L.b <- t(apply(X.b, 1, function(z) {tapply(z, L, mean)}))
           r.b <- sweep(X.bar.L.b, 1, X.bar.b)
           r.b <- rowSums(sweep(r.b^2, 2, nks, "*"), na.rm=TRUE)
           r.b <- sqrt((sum(nks) / prod(nks)) * r.b)
           s.b <- X.b - apply(X.bar.L.b,1, function(z) {rep(z, nks)})
           s.b <- sqrt((1/sum(nks-1)) * sum(1/nks) * sum(s.b^2,na.rm=TRUE))
           d.b <- r.b/(s.b+s0)
	   dB[,i] <- sort(d.b)
	}

	return(cbind(index=order.d, t=sort.d, dB))
  }
  } 
       
sam.s0 <- function(r, s, s.step=0.01, alpha.step=0.01, plot.it=FALSE)  {
  if (1/s.step >= length(s)) s.step <- 2*(1/length(s)) 
  q <- quantile(s,seq(0,1,by=s.step))
  while (any(duplicated(q)))  {
    s.step <- s.step*2
    q <- quantile(s,seq(0,1,by=s.step))
  }
  q.indices <- cut(s,q,labels=FALSE,right=FALSE)
  q.indices[which(is.na(q.indices))] <- 1/s.step
  #q.ranges <- cut(s,q,labels=NULL,right=FALSE)

  sam.d.alpha <- function(r,s,alpha) {
    s.alpha <- quantile(s,alpha)
    r/(s + s.alpha)
  }

  cv.alpha <- function(alpha) {
    d.alpha <- sam.d.alpha(r,s,alpha)

    v.j <- tapply(d.alpha,q.indices,mad)
    sd(v.j,na.rm=TRUE)/mean(v.j,na.rm=TRUE)
  }

  alpha.seq <- seq(0,1,by=alpha.step)
  cva <- sapply(alpha.seq,cv.alpha)
  if(plot.it) plot(cva,ty="l")
  alpha0 <- alpha.seq[which(cva == min(cva))]
  s0<-quantile(s,alpha0)
  i <- 2

  while (s0==0) {
    alpha0 <- alpha.seq[which(cva == cva[order(cva)[i]])]
    i <- i+1
    s0 <- quantile(s,alpha0)
  }
  s0
}



     
sam.fdr <- function(order.t, ordertB, deltas) {
  n <- length(order.t)
  ndelta <- length(deltas)
  tB <- rowMeans(ordertB)
  diff <- order.t - tB
  tmp  <- quantile(as.vector(ordertB),c(0.25,0.75), na.rm=TRUE)
  pi <- sum(order.t<tmp[2] & order.t>tmp[1], na.rm=TRUE) / (n/2)
  pi <- min(pi,1)
  table <- c()
  for (i in 1:ndelta) {
    delta <- deltas[i]
    if (sum((diff>delta)&(order.t>0))!=0) {
      pos <- min(which((diff>delta)&(order.t>0))):n
      n.pos <- length(pos)
    }
    else n.pos <- 0
    if (sum((diff<(-delta))&(order.t<0))!=0) {
      neg <-1: max(which((diff<(-delta))&(order.t<0)))
      n.neg <- length(neg)
    }
    else n.neg <- 0
    n.total <- n.pos + n.neg

    max <- ifelse(n.pos==0, Inf, order.t[n-n.pos+1])
    min <- ifelse(n.neg==0, -Inf, order.t[n.neg])
    
    median.fp <- median(apply(ordertB,2, function(z){
			sum(z>=max|z<=min)}), na.rm=TRUE)
    per90.fp <- quantile(apply(ordertB,2, function(z){
			sum(z>=max|z<=min)}),0.9, na.rm=TRUE)
    
    median.fp <- ifelse(pi*median.fp/n.total>1, 1, pi*median.fp/n.total)
    per90.fp <-  ifelse(pi*per90.fp/n.total>1, 1, pi*per90.fp/n.total)
	
    table <- rbind(table, c(delta, n.total, n.pos, n.neg, median.fp, per90.fp))
  }
  colnames(table) <- c("delta", "no.significance", "no.positive", "no.negative", "FDR(50%)", "FDR(90%)")
  
  return(table)
}


###########  F test  ######################
comp.F <- function(L=NULL){
  function(X) {
    if (is.vector(X)) X <- matrix(X, byrow=TRUE)
    newX <- test.checkX(X, L, "f")
    L <- newX$L
    X <- newX$X
     m1 <- apply(X, 1, function(z){
       (summary(lm(z ~ as.factor(L), na.action=na.omit))$f["value"])})
    return(m1)
  }
     
}

comp.modt <- function(L=NULL) {
  function(X) {
    if(is.null(L)) L <- rep(0,ncol(X))
    newX <- test.checkX(X, L, "modt")
    nc <- newX$nc
    nr <- newX$nr
    L <- newX$L
    nL <- newX$nL
    res <- .C("get_t_mod_stat",as.double(X), as.integer(nr), as.integer(nc), as.integer(L), t=single(nr), as.integer(nL), NAOK=TRUE, PACKAGE="DEDS")$t

    return(res)
  }
}

comp.B <- function(L=NULL, proportion=0.01) {
  function(X) {
    if(is.null(L)) L <- rep(0,ncol(X))
    newX <- test.checkX(X, L, "B")
    nc <- newX$nc
    nr <- newX$nr
    L <- newX$L
    nL <- newX$nL
    res <- .C("get_B",as.double(X), as.integer(nr), as.integer(nc), as.integer(L),as.integer(nL), B=single(nr), as.single(proportion), NAOK=TRUE, PACKAGE="DEDS")$B

    return(res)
  }
}

comp.ebayes <- function(L=NULL, proportion=0.01) {
  function(X) {
    if(is.null(L)) L <- rep(0,ncol(X))
    newX <- test.checkX(X, L, "ebayes")
    nc <- newX$nc
    nr <- newX$nr
    L <- newX$L
    nL <- newX$nL
    res <- .C("get_ebayes",as.double(X), as.integer(nr), as.integer(nc), as.integer(L),as.integer(nL), t=single(nr), b=single(nr), as.single(proportion), NAOK=TRUE, PACKAGE="DEDS")
    return(as.matrix(cbind(t=res$t, B=res$b)))
  }
}


comp.modF <- function(L=NULL){
  function(X) {
    newX <- test.checkX(X, L, "modf")
    nc <- newX$nc
    nr <- newX$nr
    L <- newX$L
    nL <- newX$nL
    res <- .C("get_f_mod_stat",as.double(X), as.integer(nr), as.integer(nc), as.integer(L), t=single(nr), as.integer(nL), NAOK=TRUE, PACKAGE="DEDS")$t
    return(res)
  }
}
      

test.checkX<-function(X, classlabel, test){
  if((!is.matrix(X)) || !(is.numeric(X)))
     stop("X needs to be a matrix\n")
  if(ncol(X)!=length(classlabel))
    stop("the number of column of X needs to be the same as the length of classlabel\n")
  X <- X[, order(classlabel)]
  L <- deds.checkclasslabel(classlabel,test)
  return(list(X=X, nr=nrow(X), nc=ncol(X), L=L, nL=length(unique(L))))
}

#########################################################
#### Plotting functions
########################################################3

### modified from John Fox

pairs.DEDS <- function(x, subset=c(1:nrow(x$stats)), labels=colnames(x$stats[,-1]), logit=FALSE,
                       diagonal=c("qqnorm", "boxplot", "density", "histogram", "none"),
                       lower=c("cor", "none"), groups.by.deds=TRUE, thresh=0.05,
                       reg.line=NULL, smooth=FALSE,
                       line.by.group=FALSE, diag.by.group=TRUE, lower.by.group=FALSE,
                       col=palette(), pch=1:n.groups, lwd=1,
                       legend.plot=length(levels(groups)) > 1, ...){
  stats <- x$stats[subset, -1]
  geneOrder <- x$geneOrder[subset]
  p <- x$p[subset] 
  if(logit) stats <- -log10(stats)
  if(groups.by.deds==TRUE) {
    groups <- rep("non-DE", length(geneOrder))
    if(thresh<1) noDE <- sum(p<=thresh)
    else noDE <- thresh
    sub <- 1:noDE
    groups[sub] <- "DE"
    groups<- factor(groups, levels=c("non-DE", "DE"))
  }
  else groups <- as.factor(rep(1, nrow(stats)))
  
  n.groups <- length(levels(groups))
  if (n.groups > length(col)) col <- rep(col, length=n.groups)
  if (n.groups > length(pch)) pch <- rep(pch, length=n.groups)
        
  panel.density<-function(x){
    par(new=TRUE)
    if(!diag.by.group) {
      xt <- x[!is.na(x)]
      plot(density(xt), axes=FALSE, main="", col="blue",lwd=2)
      points(xt, rep(0,length(xt)), pch="|", col="black")
    }
    else {
      xt <- x[!is.na(x)]
      gt <- groups[!is.na(x)]
      den.maxs <- c()
      for (i in 1:n.groups)  {
        subs <- gt==levels(gt)[i]
        den.maxs <- c(den.maxs, max(density(xt[subs])$y))
      }
      den.max <- max(den.maxs)
      plot(density(xt), type="n", ylim=c(0,den.max),axes=FALSE, main="")
      for (i in 1:n.groups) {
        subs <- gt==levels(gt)[i]
        lines(density(xt[subs]), col=palette()[i], lwd=2.5)
      }
   }
  }

  panel.histogram<-function(x){
    par(new=TRUE)
    hist(x, main="", axes=FALSE, col="blue")
  }

  panel.qqnorm<-function(x){
    par(new=TRUE)
    if(!diag.by.group){
      qqnorm(x, axes=FALSE, xlab="", ylab="", main="", col="blue")
      qqline(x, lty=2)
    }
    else {
      xt <- x[!is.na(x)]
      gt <- groups[!is.na(x)]
      tmp <- qqnorm(xt, axes=FALSE, plot=FALSE)
      plot(tmp$x,tmp$y, axes=FALSE, type="n")
      subs <- gt==levels(gt)[1]
      points(tmp$x[subs],tmp$y[subs], pch=20)
      for (i in 2:n.groups) {
        subs <- gt==levels(gt)[i]
        points(tmp$x[subs],tmp$y[subs], pch="*", cex=2, col=palette()[i])
      }
      qqline(x, lty=2)
   }
  }

  panel.boxplot<-function(x){
    par(new=TRUE)
    if (!diag.by.group) boxplot(x, axes=FALSE, main="", col="blue", horizontal=TRUE, boxwex=0.4)
    else boxplot(x~groups, axes=FALSE, main="", col=col, horizontal=FALSE, boxwex=0.4)
  }

  panel.blank<-function(x) NULL

  panel.group<-function(x, y, ...){
    for (i in 1:n.groups){
      subs<-groups==levels(groups)[i]
      points(x[subs], y[subs], pch=pch[i], col=col[i], cex=0.5+0.5*i)
      if (is.function(reg.line) & line.by.group)
        abline(reg.line(y[subs]~x[subs]), col=col[i], lwd=lwd)
      if (smooth & line.by.group) lines(lowess(x[subs],y[subs]), col=col[i],lwd=lwd, lty=2)
      if (is.function(reg.line) & (!line.by.group))
        abline(reg.line(y~x), col=col[4], lwd=lwd)
      if (smooth & (!line.by.group)) lines(lowess(x, y), col=col[4],lwd=lwd, lty=2)
    }
  }

  
  panel.cor <- function(x, y, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if (!lower.by.group) {
      r <- cor(x, y, use='c')
      text(0.5, 0.5, as.character(round(r,2)), cex = 2)
    }
    else {
      for (i in 1:n.groups) {
        subs <- groups==levels(groups)[i]
        r <- cor(x[subs], y[subs], use='c')
        text(0.5, i/(n.groups+1), as.character(round(r,2)), cex = 1.5, col=col[i])
      }
    }
  }

  which.fn<-match(match.arg(diagonal), c("density", "histogram", "boxplot","qqnorm", "none"))
  diag<-list(panel.density, panel.histogram, panel.boxplot, panel.qqnorm, panel.blank)[[which.fn]]
  which.fn<-match(match.arg(lower), c("cor", "none"))
  lower<-list(panel.cor, panel.group)[[which.fn]]
            
  pairs.default(stats, labels=labels, diag.panel=diag, panel=panel.group, lower.panel=lower)
  if(legend.plot) {
    frac<-1/ncol(stats)
    legend(1 - 0.95*frac, 0.8*frac, legend=levels(groups), pch=pch, col=col, 
          cex=cumprod(par("fin"))[2]*sqrt(frac)/(sqrt(n.groups)*20))
  }
}

qqnorm.DEDS <- function(y, subset=c(1:nrow(y$stats)),
                        xlab = "Quantiles of standard normal", thresh=0.05,
                        col=palette(), pch,... ) {
  stats <- y$stats[subset, -1]
  geneOrder <- y$geneOrder[subset]
  p <- y$p[subset]
  nT <- NCOL(stats)
  is.stat <- y$options[1]=="t"
  if(!is.stat) stop("DEDS not summarizing statistics, try the function hist.DEDS...")
  if(is.null(colnames(stats))) ylabs <- paste("Stat", c(1:nT), ": sample quantiles", sep="")
  else ylabs <- paste(colnames(stats), ": sample quantiles",sep="")

  groups <- rep("nsig", length(geneOrder))
  if(thresh<1) sub <- which(p<=thresh)
  else sub <- 1:thresh
  groups[sub] <- "sig"
  groups<- factor(groups, levels=c("nsig", "sig"))
  if(length(col)<2) col <- c(col, col+1)
  if(missing(pch)) pch <- c(20,8)
  else if(length(pch)<2) pch<- c(pch, pch+1)

  par(mfrow=c(2, floor((nT+1)/2)))
  for(i in 1:nT) {
    notNA <- !is.na(stats[,i])
    gt <- groups[notNA]
    xt <- stats[notNA,i]
    tmp <- qqnorm(xt, plot=FALSE, ...)
    plot(tmp$x, tmp$y, type="n", xlab = xlab, ylab=ylabs[i])
    points(tmp$x[gt=="nsig"],tmp$y[gt=="nsig"], pch=pch[1], col=col[1])
    points(tmp$x[gt=="sig"], tmp$y[gt=="sig"], pch=pch[2], col=col[2], cex=1.5) 
    qqline(xt, col="gray", lwd=2, lty=2)
  }
}

hist.DEDS <- function(x, subset=c(1:nrow(x$stats)), ...) {
  ps <- x$stats[subset, -1]
  nP <- NCOL(ps)
  is.p <- x$options[1]=="p"
  if(!is.p) stop("DEDS not summarizing p values, try the function qqnorm.DEDS...")
  if(is.null(colnames(ps))) xlabs <- paste("Model", c(1:nP), ": p values", sep="")
  else xlabs <- paste(colnames(ps), ": p values", sep="")

  args <- list(...)
  if("col" %in% names(args)) {
    col <- args$col
    args <- args[-which(names(args)=="col")]
  }
  else col <- "gray"
  if("border" %in% names(args)) {
    border <- args$border
    args <- args[-which(names(args)=="border")]
  }
  else border <- "gray"
  if("nclass" %in% names(args)) {
    nclass <- args$nclass
    args <- args[-which(names(args)=="nclass")]
  }
  else nclass <- 50
  if("main" %in% names(args)) {
    main <- args$main
    args <- args[-which(names(args)=="main")]
  }
  else main <- ""
  

  par(mfrow=c(2, floor((nP+1)/2)))
  for(i in 1:nP) {
    p <- ps[,i]
    do.call("hist", c(list(x=p, xlab=xlabs[i], nclass=nclass, col=col,
                           border=border, main=main), args))
    #hist(p, xlab = xlabs[i], nclass = nclass, col = col, border = border, cex = cex, main=main, ...)
    abline(h=nrow(p)/50, col="black", lwd=4, lty=2)
  }
}
