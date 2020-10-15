#--------------------------------------------------------------------------------------
# - plot2d:     plotta un dataset 2D, con la media. 
#               Se pca=T disegna gli assi delle componenti principali e l'ellissi 
#               indotta dalla matrice di covarianza.
#               Se gaussheat=T plotta la heatmap di una gaussiana sotto la nuvola di punti
#
# - qqhist:     dato il dataset X plotta istogramma vs curva gaussiana teorica,
#               qqnorm e qqline delle colonne di X. 
#               Se pca=T fa la stessa cosa per le componenti principali.
#               Se shapiro=T, plotta lo shapiro test per ogni vettore che viene plottato.
#
# - mahal.test: Effettua il test di goodness of fit sulla mahalanobis distance. Sotto H0
#               la mahal dist viene da una Chi2(p) e dunque i dati X sono gaussiani.
#            
# - mcshapiro.test: Metodo Monte Carlo per effettuare shapiro test su tutte le direzioni.
#               Approssima con MC la distribuzione di min(W) sotto H0: X gaussiano, e calcola
#               il pvalue della min(W) effettivamente osservata. 
# - cut.ouliers: Dato X, plotta le Mahalanobis dist delle righe di X colorando a vari livelli 
#                di frequenza (quantili emipirici) in modo da osservare visivamente gli outliers.
#                Restituisce poi d2.level con cui estrarre da X i dati non outlier.
#              e.g.: d2.lvl=cut.ouliers(X); X.nooutliers = X[which(d2<d2.lvl[5])]
#
#--------------------------------------------------------------------------------------


#-------------------
library(mvtnorm)
library(mvnormtest)
library(car)
#-------------------
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

plot2d = function(X, gaussheat=F, pca=T, ellipse=T){
    
    x11() 
    if(gaussheat){
    
        lb.x1 = min(X[,1]) - sd(X[,1]);
        ub.x1 = max(X[,1]) + sd(X[,1]);
        lb.x2 = min(X[,2]) - sd(X[,2]);
        ub.x2 = max(X[,2]) + sd(X[,2]);
        print(lb.x1)
        print(ub.x1)
        print(lb.x2)
        print(ub.x2)
    
        x.1 <- seq(lb.x1, ub.x1, length=100); 
        x.2 <- seq(lb.x2,ub.x2,length=100);
        w <- matrix(NA, length(x.1), length(x.2))

        for(i in 1:length(x.1)){  
            for(j in 1:length(x.2)){
                w[i,j] <- dmvnorm(c(x.1[i],x.2[j]),mu,sig)
            }
        }
        image(x.1, x.2, w,asp=1, main="Sample points")
        points(X[,1],X[,2],pch=20,cex=.75)
    }
    else{
        M <- colMeans(X)
        plot(X, asp=1, xlab='X.1', ylab='X.2',pch=1)
        points(M[1],M[2], pch=17, cex=2.5, col='red')
        segments(M[1],-1e6,M[1],M[2], col='red',cex=5,lty='dashed')        
        segments(-1e6,M[2],M[1],M[2], col='red',cex=5, lty='dashed')
        if(pca){
            
            S <- cov(X)
            ellipse(M, S, 1, add=T,lwd=1, col = 'blue')
            ev = eigen(S)$vectors
            abline(a= M[2]- (ev[2,1]/ev[1,1]) *M[1], b= ev[2,1]/ev[1,1], col = 'navyblue', lwd = 2)
            abline(a= M[2]- (ev[2,2]/ev[1,2]) *M[1], b= ev[2,2]/ev[1,2], col = 'navyblue', lwd = 2)
            
            
        }
    }
  
  
}

qqhist = function(X, shapiro=TRUE, pca=TRUE){
    p = dim(cbind(X))[2]

    # Lavoro sulle variabili originarie X1,...,Xp
    x11()
    par(mfrow=c(2,p))
    
    y.ub = max( unlist(apply( X, 2,  function(x) hist(x,plot=F)$density)) ) 
    x.dis = seq(min(X),max(X), by=min(abs(X))/10)
    
    for( i in 1:p){
      
      

        name = paste("X.",i,sep="")
        htit = paste("Histogram of", name)
        hist(X[,i], prob=T, ylab='density', xlab=name, main=htit ,ylim=c(0,y.ub))
        lines(x.dis, dnorm(x.dis,mean(X[,i]),sd(X[,i])), col='blue', lty=2)

        qqtit = paste("QQplot of",name)
        qqnorm(X[,i], main=qqtit ,xlab='theoretical quantiles', ylab='sample quantiles')
        qqline(X[,i])

        if(shapiro){
            pv=shapiro.test(X[,i])
            print( paste("Shapiro test", name, "pvalue =", toString(pv$p.value) )) 

        }
    }

    if(pca){
        # Lavoro sulle componenti principali 
        x11()
        par(mfrow=c(2,p))

        pc.X <- princomp(X,scores=T)
        y.ub = max( unlist(apply( pc.X$scores, 2,  function(x) hist(x,plot=F)$density)) ) 
        x.dis = seq(min(pc.X$scores),max(pc.X$scores), by=min(abs(pc.X$scores))/10)
        
        for( i in 1:p){

            name = paste("PC.",i,sep="")

            htit = paste("Histogram of", name)
            hist(pc.X$scores[,i], prob=T, ylab='density', xlab=name, main=htit ,ylim=c(0,y.ub))
            lines(x.dis, dnorm(x.dis,mean(pc.X$scores[,i]),sd(pc.X$scores[,i])), col='blue', lty=2)

            qqtit = paste("QQplot of",name)
            qqnorm(pc.X$scores[,i], main=qqtit ,xlab='theoretical quantiles', ylab='sample quantiles')
            qqline(pc.X$scores[,i])

            if(shapiro){
                pv=shapiro.test(pc.X$scores[,i])
                print( paste("Shapiro test", name, "pvalue =", toString(pv$p.value) )) 

            }
        }
    }
}




mcshapiro.test = function(X, devstmax=0.01, sim=ceiling(1/(4*devstmax^2))){
  library(mvnormtest)
  n   <- dim(X)[1]
  p   <- dim(X)[2]
  mu  <- rep(0,p)
  sig <- diag(p)
  W   <- NULL
  for(i in 1:sim)
  {
    Xsim <- rmvnorm(n, mu, sig)
    W   <- c(W, mshapiro.test(t(Xsim))$stat)
    # mshapiro.test(X): compute the statistics min(W) for the sample X
  }
  Wmin   <- mshapiro.test(t(X))$stat   # min(W) for the given sample
  pvalue <- sum(W < Wmin)/sim          # proportion of min(W) more extreme than the observed Wmin
  devst  <- sqrt(pvalue*(1-pvalue)/sim)
  list(Wmin = as.vector(Wmin), pvalue = pvalue, devst = devst, sim = sim)
}


mahal.test = function(X){
  p=dim(cbind(X))[2]
  n=dim(cbind(X))[1]
  d2 = mahalanobis(X, colMeans(X), cov(X))  #Mi aspetto sia distribuita come una Chi2(p)
  x.dis = seq(min(d2),max(d2), by=min(abs(d2))/10)
  
  
  x11()
  par(mfrow=c(1,2))
  hist(d2, prob=T, main='Histogram of the Mahalanobis dist.',xlab='d2',ylab='density', col='grey84')
  lines(x.dis, dchisq(x.dis), col='blue', lty=2, lwd=2)
  
  qqplot(qchisq((1:n - 0.5)/n, df = p), d2, main='QQplot of (sample) d2',xlab='theoretical quantiles Chi-sq(2)',
         ylab='sample quantiles')
  abline(0, 1)

  
  d2.class <- cut(d2, qchisq((0:10)/10, df = p))
  d2.freq  <- table(d2.class)
  
  chisq.test(x = d2.freq, p = rep(1/10, 10), simulate.p.value = T) # test sulla goodness of fit, H0: dati Chisq(p)
}


# ricerca degli outliers 

cut.outliers = function(X){
  
  M <- colMeans(X)
  S <- cov(X)
  d2 <- matrix(mahalanobis(X, M, S))
  d2.q= quantile(d2, probs = c(0,0.90,0.95,0.97,0.99,1))
  
  d2.level = cut(d2, breaks=d2.q, 
                labels=c("0-0.90","0.9-0.95","0.95-0.97","0.97-0.99","0.99-1"), include.lowest=TRUE)
  
  x11()
  colors <- c("#339900", "#F6CB00", "#ffcc00","#ff9966","#cc3300")
  colors <- colors[as.numeric(d2.level)]

  
  plot(d2, asp=1, main = "Mahal. dist", ylab = "Mahal. dist", col = colors, pch=16, cex=2)
  legend("topright", legend = levels(d2.level), 
         col=c("#339900", "#F6CB00", "#ffcc00","#ff9966","#cc3300"), pch=16, pt.cex=3)

  X11()
  plot(X, asp=1, pch=16, cex=1.6,col= ifelse( d2<d2.q[2] ,"#339900",
                                        ifelse(d2 < d2.q[3], "#F6CB00",
                                               ifelse(d2 < d2.q[4],"#ffcc00",
                                                      ifelse(d2 <d2.q[5],"#ff9966","#cc3300")
                                                     )
                                              )
                                            )
  )
  
  return( d2.q )
}
  

## Trasformazione dei dati tramite Box-Cox

# memo: i dati devono essere positivi. Traslarli prima di applicare B-C. 
# memo: lambda>1 
# lambda.x <- powerTransform(x)  // x univariata. Trova il miglior lambda di Box-Cox
# bc.x <- bcPower(x, lambda.x$lambda)  // Trasforma x con Box-Cox usando come lambda lambda.x$lambda

# multivariate case:
# lambda <- powerTransform(cbind(x,y))
# BC.x <- bcPower(x, lambda$lambda[1])
# BC.y <- bcPower(y, lambda$lambda[2])
# OCIO: trasformare su alcune direzioni mica basta! Serve la normalitá su tutte le direzioni! Ma piuttosto che niente...



# Inferenza sulla media

norm.mean.test = function(x, mu0, alpha=0.05, plot.region=T, plot.pval=T){
  #
  # TODO: plot.pval and plot.region
  #
  mu0=cbind(mu0)
  x = cbind(x);
  p = dim(x)[2];
  n = dim(x)[1];
  x.cov = var(x);
  x.invcov = solve(x.cov);
  x.mean = colMeans(x)
  if(dim(mu0)[1]!=p){print("ERROR: wrong length of mu0");}
  
  x.T2 = function(mu){(n-p)/(p*(n-1)) *n* t(x.mean-mu)%*%x.invcov%*%(x.mean-mu)}
  x.T2.0 = x.T2(mu0)
  
  test.alpha = x.T2.0 > qf(1-alpha, p,n-p)
  print( paste("Test level", alpha, ":", "x.T2.0 =", x.T2.0, ifelse(test.alpha, ">","<"), qf(1-alpha, p,n-p),"=", "qf(1-alpha, p,n-p). ",ifelse(test.alpha, "Evidence to REJECT H0. Accept mu!=mu0", "NO evidence to reject H0. Still accept mu=mu0 ") ) )
  
  pval = 1-pf( x.T2.0, p,n-p)
  
  print(paste( "Test p-value =", pval))
  
  ## confidence region
  print(paste( "Confidence region characterization: {mu in R^",p," : (n-p)/(p*(n-1)) *n* t(x.mean-mu)%*%x.invcov%*%(x.mean-mu) < qf(1-alpha, p,n-p)","}", sep=""))
  
  
  # MEMO: x'Ax = c ==> princ.dirs = eigvects(A), princ.lens = sqrt( c/eigvals(A) )
  ## in our case A=x.invcov, c = qf(1-alpha,p,n-p)*(p*(n-1))/(n-p)/n [rmk: we work con x.cov to prevent round-offs, but it's the same]
  CR.principal.dirs = eigen(x.cov)$vectors
  CR.principal.lens = sqrt( qf(1-alpha,p,n-p)*(p*(n-1))/(n-p)/n * eigen(x.cov)$val) 
  print(paste("info about the CR(mu) of confidence", 1-alpha,": "))
  print("  - Confidence region principal axes directions: ")
  print(CR.principal.dirs)
  print("  - Confidence region principal axes lengths: ")
  print(CR.principal.lens)
  
}


# Asymptotical mean test for the mean of non gaussian data 
asy.mean.test = function(x, mu0, alpha=0.05){
  
  mu0=cbind(mu0)
  x = cbind(x);
  p = dim(x)[2];
  n = dim(x)[1];
  x.cov = var(x);
  x.invcov = solve(x.cov);
  x.mean = matrix(sapply(x, mean))
  if(dim(mu0)[1]!=p){print("ERROR: wrong length of mu0");}
  if(p>n-10^p){print("WARNING: n not so large to use an asymptotical test...")}
  
  x.T2 = function(mu){ n*mahalanobis(t(mu),x.mean, x.invcov, inverted = TRUE)}
  x.T2.0 = x.T2(mu0)
  
  test.alpha = x.T2.0 > qchisq(1-alpha,df=p)
  print( paste("Asymptotical Test level", alpha, ":", "x.T2.0 =", x.T2.0, ifelse(test.alpha, ">","<"), qchisq(1-alpha,df=p),"=", "qchisq(1-alpha,df=p). ",ifelse(test.alpha, "Evidence to REJECT H0. Accept mu!=mu0", "NO evidence to reject H0. Still accept mu=mu0 ") ) )
  
  pval = 1 - pchisq(x.T2.0,df=p)
  
  print(paste( "Test p-value =", pval))
  
  ## confidence region
  print(paste( "Confidence region characterization: {mu in R^",p," : n* t(x.mean-mu)%*%x.invcov%*%(x.mean-mu) < qchisq(1-alpha,df=p)","}", sep=""))
  
  
  # MEMO: x'Ax = c ==> princ.dirs = eigvects(A), princ.lens = sqrt( c/eigvals(A) )
  ## in our case A=x.invcov, c = qf(1-alpha,p,n-p)*(p*(n-1))/(n-p)/n [rmk: we work con x.cov to prevent round-offs, but it's the same]
  CR.principal.dirs = eigen(x.cov)$vectors
  CR.principal.lens = sqrt( qchisq(1-alpha,df=p)/n * eigen(x.cov)$val) 
  print(paste("info about the asyCR(mu) of confidence", 1-alpha,": "))
  print("  - AsyConfidence region principal axes directions: ")
  print(CR.principal.dirs)
  print("  - AsyConfidence region principal axes lengths: ")
  print(CR.principal.lens)
  
}


norm.mean.oneAtTheTime.CI = function(x, dir, alpha=0.05, H1='='){
  x = cbind(x)
  x.mean = colMeans(x)
  x.cov = var(x)
  dir = cbind(dir)
  p = dim(x)[2]
  n = dim(x)[1]
  
  # Test stat: T2 = function(mu){return((t(dir)%*%x.mean - t(dir)%*%mu)/sqrt(t(dir)%*%x.cov%*%dir) * sqrt(n))}
  
  if(H1=='='){
    test.quantile = qt(1-alpha/2, df=n-1)
    
    CI  =  c( t(dir)%*%x.mean - test.quantile*sqrt(t(dir)%*%x.cov%*%dir)/sqrt(n),
              t(dir)%*%x.mean + test.quantile*sqrt(t(dir)%*%x.cov%*%dir)/sqrt(n))
    
    print(paste("CI One-at-the-Time at level", 1-alpha,"for t(dir)*mu = mu0 vs t(dir)*mu != mu0:"))
    print(CI)
  }
  if(H1=='>'){
    test.quantile = qt(1-alpha, df=n-1)
    
    CI  =  c( t(dir)%*%x.mean - test.quantile*sqrt(t(dir)%*%x.cov%*%dir)/sqrt(n), Inf)
    
    print(paste("CI One-at-the-Time at level", 1-alpha,"for t(dir)*mu = mu0 vs t(dir)*mu > mu0:"))
    print(CI)
  }
  if(H1=='<'){
    test.quantile = qt(1-alpha, df=n-1)
    
    CI  =  c( -Inf,t(dir)%*%x.mean + test.quantile*sqrt(t(dir)%*%x.cov%*%dir)/sqrt(n))
    
    print(paste("CI One-at-the-Time at level", 1-alpha,"for t(dir)*mu = mu0 vs t(dir)*mu < mu0:"))
    print(CI)
    
  }
  
  if(H1!='=' && H1!='>' && H1!='<'){
    print("ERROR: invalid H1")
    return(0)
  }
  return(CI)
}

norm.mean.sym.CI = function(x, dir, alpha=0.05){
  x = cbind(x)
  x.mean = colMeans(x)
  x.cov = var(x)
  dir = cbind(dir)
  p = dim(x)[2]
  n = dim(x)[1]
  
  # Test stat: T2 = function(mu){return((t(dir)%*%x.mean - t(dir)%*%mu)/sqrt(t(dir)%*%x.cov%*%dir) * sqrt(n))}
  
  test.quantile = qf(1-alpha,p,n-p)
  
  CI  =  c( t(dir)%*%x.mean - sqrt( test.quantile * p*(n-1)/(n-p) * t(dir)%*%x.cov%*%dir/n),
            t(dir)%*%x.mean + sqrt( test.quantile * p*(n-1)/(n-p) * t(dir)%*%x.cov%*%dir/n))
  
  print(paste("symCI at level", 1-alpha,"for t(dir)*mu:",CI))
  return(CI)
}

norm.mean.CI.onAxes = function(x){
  
  x=cbind(x)
  p=dim(x)[2]
  I=diag(p)
  DIRlist = NULL
  for(i in 1:p)DIRlist=append(DIRlist, list(I[,i]))
  bfCIlist = NULL
  symCIlist = NULL
  for(i in 1:p){
    bfCI = quiet(norm.mean.oneAtTheTime.CI(x,dir=DIRlist[[i]],alpha = alpha/p,H1 = '='))
    bfCIlist=append(bfCIlist,list(bfCI))
    
    symCI = quiet(norm.mean.sym.CI(x, dir = DIRlist[[i]],alpha = alpha))
    symCIlist = append(symCIlist,list(symCI))
  }
  
  print("Bonferroni CI on axes:")
  print(data.frame("lb"=matrix(unlist(bfCIlist),p,2, byrow = T)[,1],"ub"=matrix(unlist(bfCIlist),p,2,,byrow = T)[,2]))
  print("Sym CI on axes:")
  print(data.frame("lb"=matrix(unlist(symCIlist),p,2, byrow = T)[,1],"ub"=matrix(unlist(symCIlist),p,2,,byrow = T)[,2]))
  plot.CIlist(CIlist = bfCIlist, CI_longer.list = symCIlist, DIRlist =DIRlist )
  
  
}


plot.CIlist = function(CIlist, CI_longer.list=NULL, DIRlist=NULL, mu0=NULL, main='Confidence Intervals', xlab='Variables',ylab='Confidence intervals'){
  n.CI = length(CIlist)
  ylim=range(unlist(CIlist))
  if(!is.null(CI_longer.list) && n.CI==length(CI_longer.list)){
    ylim = range(c(unlist(CI_longer.list), unlist(CIlist)))
  }
  x11()
  matplot(1:n.CI,1:n.CI,pch='',ylim=ylim, xlim=c(0,n.CI+1),xlab=xlab, ylab=ylab,main=main)
  
  if(!is.null(CI_longer.list) && n.CI==length(CI_longer.list)){
    for(i in 1:n.CI) segments(i+0.1,CI_longer.list[[i]][1],i+0.1,CI_longer.list[[i]][2],lwd=2,col="navyblue")
    for(i in 1:n.CI) points(i+0.1, mean(CI_longer.list[[i]]) ,lwd=3,pch=16,col='orange')
  }
  
  for(i in 1:n.CI) segments(i,CIlist[[i]][1],i,CIlist[[i]][2],lwd=2,col=i)
  for(i in 1:n.CI) points(i, mean(CIlist[[i]]) ,lwd=3,pch=16,col='orange')
  
  if(!is.null(mu0) && !is.null(DIRlist)){
    for(i in 1:n.CI) points(i, t(cbind(mu0))%*%cbind(DIRlist[[i]]),lwd=3,pch=16,col='purple')
    
  }
  
}

norm.repeated.test= function(x, C, Cmu0, alpha=0.05){
  x11()
  matplot(t(x), type='l')
  
  Cmu0=cbind(Cmu0)
  x = cbind(x);
  C = cbind(C);
  
  q = dim(x)[2];
  n = dim(x)[1];
  x.mean = sapply(x, mean);
  x.cov = var(x);
  CcovC.inv = solve(C%*%x.cov%*%t(C))
  
  
  if(dim(Cmu0)[1]!=q-1){print("ERROR: wrong length of Cmu0");}
  
  x.T2 = function(Cmu){n* t(C%*%x.mean - Cmu)%*% CcovC.inv %*% (C%*%x.mean - Cmu) *(n-q+1)/((n-1)*(q-1))}
  x.T2.0 = x.T2(Cmu0)
  
  test.alpha = x.T2.0 > qf(1-alpha, q-1,n-q+1)
  print( paste("Contrast Test level", alpha, ":", "x.T2.0 =", x.T2.0, ifelse(test.alpha, ">","<"), qf(1-alpha, q-1,n-q+1),"=", "qf(1-alpha, q-1,n-q+1). ",ifelse(test.alpha, "Evidence to REJECT H0. Accept Cmu!=Cmu0", "NO evidence to reject H0. Still accept Cmu=Cmu0 ") ) )
  
  pval = 1-pf( x.T2.0, q-1,n-q+1)
  
  print(paste( "Test p-value =", pval))
  
  ## confidence region
  print(paste( "Confidence region characterization: {Cmu in R^",q-1," : (n-q+1)/((q-1)*(n-1)) *n* t(C*x.mean-Cmu)%*%solve(C*cov*C')%*%(C*x.mean-Cmu) <  qf(1-alpha, q-1,n-q+1)","}", sep=""))
  
  
  # MEMO: x'Ax = c ==> princ.dirs = eigvects(A), princ.lens = sqrt( c/eigvals(A) )
  ## in our case A=x.invcov, c = qf(1-alpha,p,n-p)*(p*(n-1))/(n-p)/n [rmk: we work con x.cov to prevent round-offs, but it's the same]
  CR.principal.dirs = eigen(x.cov)$vectors
  CR.principal.lens = sqrt(qf(1-alpha, q-1,n-q+1)*((n-1)*(q-1))/(n-q+1)/n * eigen(x.cov)$val) 
  print(paste("info about the CR(Cmu) of confidence", 1-alpha,": "))
  print("  - Confidence region principal axes directions: ")
  print(CR.principal.dirs)
  print("  - Confidence region principal axes lengths: ")
  print(CR.principal.lens)
  
}





