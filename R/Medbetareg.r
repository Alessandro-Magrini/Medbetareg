# rescaling procedure
rescale <- function(v,vscale=NA) {
  if(identical(vscale,NA)) {
    if(sum((v!=0&v!=1),na.rm=T)>0) stop("Invalid value",call.=F)
    vtil <- v
    } else {
    if(length(vscale)!=4 | (identical(sort(vscale),vscale)==F&identical(sort(vscale,T),vscale)==F)) {
      stop("Invalid scale",call.=F) 
      }
    if(identical(sort(vscale,T),vscale)) {
      vscale <- -vscale
      v <- -v
      }
    if(sum(v<vscale[1],na.rm=T)>0 | sum(v>vscale[4],na.rm=T)>0) stop("Invalid value",call.=F)
    vtil <- c()
    for(i in 1:length(v)) {
      if(!is.na(v[i])) {
        if(v[i]<vscale[2]) vtil[i] <- -1.5+(v[i]-vscale[1])/(vscale[2]-vscale[1])
        if(v[i]>=vscale[2]&v[i]<=vscale[3]) vtil[i] <- -0.5+(v[i]-vscale[2])/(vscale[3]-vscale[2])
        if(v[i]>vscale[3]) vtil[i] <- 0.5+(v[i]-vscale[3])/(vscale[4]-vscale[3])
        } else {
        vtil[i] <- NA
        }
      }
    }
  vtil
  }

# inverted rescaling procedure
invRescale <- function(v,vscale=NA) {
  if(identical(vscale,NA)) {
    if(sum((v!=0&v!=1),na.rm=T)>0) stop("Invalid value",call.=F)
    vtil <- v
    } else {
    if(sum(v< -1.5 | v>1.5)>0) stop("Invalid value",call.=F)
    if(length(vscale)!=4 | (identical(sort(vscale),vscale)==F&identical(sort(vscale,T),vscale)==F)) {
      stop("Invalid scale",call.=F) 
      }
    vtil <- c()
    for(i in 1:length(v)) {
      if(!is.na(v[i])) {
        if(v[i]< -0.5) vtil[i] <- vscale[1]+(v[i]+1.5)*(vscale[2]-vscale[1])
        if(v[i]>= -0.5&v[i]<=0.5) vtil[i] <- vscale[2]+(v[i]+0.5)*(vscale[3]-vscale[2])
        if(v[i]>0.5) vtil[i] <- vscale[3]+(v[i]-0.5)*(vscale[4]-vscale[3])
        } else {
        vtil[i] <- NA
        }
      }
    }
  vtil
  }

# check if a variable name is correctly spelled (internal use only)
checkName <- function(nome) {
  n <- nchar(nome)
  ifelse(n>20 | length(intersect(substr(nome,1,1),0:9))>0, ok<-F, ok<-T)
  if(ok==T) {
    for(i in 1:n) {
      if(length(intersect(substr(nome,i,i),c(0:9,letters,toupper(letters),".","_")))==0) {
        ok <- F
        break()
        }
      }
    }  
  ok
  }

# nu function (internal use only))
nuFun <- function(z) {
  r <- log(log(5)/log(2))/log(2)
  sign(z)*abs(z)^r
  }

# create an object of class 'mbr'
newPrior <- function(assess.code,nrep=100,digits=6) {
  if(!identical(round(nrep),nrep) | nrep<100) stop("Argument 'nrep' must be an integer number no less than 100",call.=F)
  if(!identical(round(digits),digits) | digits<0) stop("Argument 'digits' must be a non negative number",call.=F)
  auxmod <- strsplit(gsub("\n",";",assess.code),";")[[1]]
  auxmod <- auxmod[which(nchar(auxmod)>0)]
  if(!is.character(assess.code) || length(assess.code)!=1 || length(auxmod)<1) {
    stop("The model code has zero lines",call.=F)
    }
  cnam <- bnam <- auxord <- c()
  mlist <- list()
  for(i in 1:length(auxmod)) {
    iaux <- strsplit(auxmod[i]," ")[[1]]
    mlist[[i]] <- iaux[which(nchar(iaux)>0)]
    }
  y_code <- mlist[sapply(mlist,function(x){x[1]=="RESP"})] 
  cx_code <- mlist[sapply(mlist,function(x){x[1]=="CEV"})] 
  bx_code <- mlist[sapply(mlist,function(x){x[1]=="BEV"})] 
  ix_code <- mlist[sapply(mlist,function(x){x[1]=="INTER"})] 
  tau_code <- mlist[sapply(mlist,function(x){x[1]=="TAU"})] 
  if(length(y_code)!=1) stop("The model code must contain one and only one RESP statement",call.=F)
  if(length(y_code[[1]])!=3) stop("Invalid RESP statement: ",length(y_code[[1]])-1," arguments provided instead of 2",call.=F) 
  if(checkName(y_code[[1]][2])==F) stop("Invalid variable name in RESP statement: ",y_code[[1]][2],call.=F) 
  yscale <- eval(parse(text=paste("c",y_code[[1]][3],sep="")))
  if(!is.numeric(yscale)) stop("Non-numerical scale for variable ",y_code[[1]][2],call.=F)
  if(length(yscale)!=4 || !identical(yscale,sort(yscale))&&!identical(yscale,sort(yscale,T))) stop("Invalid scale for variable ",y_code[[1]][2],call.=F)
  if(length(tau_code)!=1) stop("The model code must contain one and only one TAU statement",call.=F)
  if(length(tau_code[[1]])!=3) stop("Invalid TAU statement: ",length(tau_code[[1]])-1," arguments provided instead of 2",call.=F)
  taupos <- eval(parse(text=tau_code[[1]][2]))
  if(taupos<0 || taupos>1) stop("Invalid relative position in TAU statement",call.=F)
  if((tau_code[[1]][3] %in% c("lp","n","hp"))==F) stop("Unknown range in TAU statement",call.=F)
  if(length(cx_code)==0 & length(bx_code)==0) stop("The model code must contain at least one CEV or BEV statement",call.=F)
  scales <- list()
  scales[[1]] <- yscale
  auxass <- c()
  if(length(cx_code)>0) {
    for(i in 1:length(cx_code)) {  
      if(checkName(cx_code[[i]][2])==F) stop("Invalid variable name in CEV statement: ",cx_code[[i]][2],call.=F) 
      if(length(cx_code[[i]])!=8) stop("Invalid CEV statement for variable ",cx_code[[i]][2],": \n  ",length(cx_code[[i]])-1," arguments provided instead of 7",call.=F)
      iscal <- eval(parse(text=paste("c",cx_code[[i]][3],sep="")))
      if(!is.numeric(iscal)) stop("Non-numerical scale for variable ",cx_code[[i]][2],call.=F)
      if(length(iscal)!=4 | (identical(sort(iscal),iscal)==F&identical(sort(iscal,T),iscal)==F)) {
        stop("Invalid scale for variable ",cx_code[[i]][2],call.=F)
        }
      scales[[i+1]] <- iscal
      irel1 <- eval(parse(text=cx_code[[i]][4]))
      irel2 <- eval(parse(text=cx_code[[i]][6]))
      iss <- eval(parse(text=cx_code[[i]][8]))
      if(irel1<0 || irel1>1 || irel2<0 || irel2>1) stop("Invalid relative position in CEV statement for variable ",cx_code[[i]][2],call.=F)
      if((cx_code[[i]][5] %in% c("lp","n","hp"))==F || (cx_code[[i]][7] %in% c("lp","n","hp"))==F) stop("Unknown range in CEV statement for variable ",cx_code[[i]][2],call.=F)
      if(iss<=0) stop("Invalid virtual cases in CEV statement for variable ",cx_code[[i]][2],call.=F)
      auxass <- rbind(auxass,c(irel1,cx_code[[i]][5],irel2,cx_code[[i]][7],iss))
      }
    cnam <- sapply(cx_code,function(x){x[2]})
    }
  if(length(bx_code)>0) {
    for(i in 1:length(bx_code)) {
      if(checkName(bx_code[[i]][2])==F) stop("Invalid variable name in BEV statement: ",bx_code[[i]][2],call.=F)
      if(length(bx_code[[i]])!=5) stop("Invalid BEV statement for variable ",bx_code[[i]][2],": \n  ",length(bx_code[[i]])-1," arguments provided instead of 4",call.=F)
      scales[[length(cx_code)+i+1]] <- NA
      irel <- eval(parse(text=bx_code[[i]][3]))
      iss <- eval(parse(text=bx_code[[i]][5]))
      if(irel<0 || irel>1) stop("Invalid relative position in BEV statement for variable ",bx_code[[i]][2],call.=F)
      if((bx_code[[i]][4] %in% c("lp","n","hp"))==F) stop("Unknown range in BEV statement for variable ",bx_code[[i]][2],call.=F)
      if(iss<=0) stop("Invalid virtual cases in BEV statement for variable ",bx_code[[i]][2],call.=F)  
      auxass <- rbind(auxass,c(NA,NA,irel,bx_code[[i]][4],iss))
      }
    bnam <- sapply(bx_code,function(x){x[2]})
    }
  if(length(cx_code)>0 || length(bx_code)>0) {
    auxord <- sapply(mlist[sapply(mlist,function(x){x[1]=="CEV"||x[1]=="BEV"})],function(x){x[2]})
    }
  auxdupl <- c(y_code[[1]][2],auxord)[duplicated(c(y_code[[1]][2],auxord))]
  if(length(auxdupl)>0) stop("Duplicated variable names: ",paste(auxdupl,collapse=", "),call.=F)
  names(scales) <- c(y_code[[1]][2],cnam,bnam)
  scales <- scales[c(y_code[[1]][2],auxord)]
  rownames(auxass) <- c(cnam,bnam)
  auxass <- auxass[auxord,]
  beta.hat <- xtil <- ytil <- q.hat <- c()
  for(i in 1:length(auxord)) {
    if(identical(scales[[auxord[i]]],NA)) {
      xtil[i] <- 1
      } else {
      xtil[i] <- ifelse(auxass[auxord[i],2]=="lp",-1.5,ifelse(auxass[auxord[i],2]=="n",-0.5,0.5))+as.numeric(auxass[auxord[i],1])
      }
    ytil[i] <- ifelse(auxass[auxord[i],4]=="lp",-1.5,ifelse(auxass[auxord[i],4]=="n",-0.5,0.5))+as.numeric(auxass[auxord[i],3])
    beta.hat[i] <- log((1.5+ytil[i])/(1.5-ytil[i]))/(log(5)*nuFun(xtil[i]))
    q.hat[i] <- as.numeric(auxass[auxord[i],5])
    }
  names(beta.hat) <- names(xtil) <- names(ytil) <- names(q.hat) <- auxord
  betaInt.hat <- intnam <- intess <- xtilInt <- c()
  if(length(ix_code)>0) {
    for(i in 1:length(ix_code)) {
      ilen <- length(ix_code[[i]])
      if(ilen<6) stop("Invalid INTER statement ",paste(ix_code[[i]],collapse=" "),":\n  ",ilen-1," arguments provided against a minimum of 5",call.=F)
      iev <- ix_code[[i]][2:(ilen-3)]
      iche <- setdiff(iev,auxord)
      if(length(iche)>0) stop("Unknown variable in INTER statement: ",paste(iche,collapse=", "),call.=F)
      irel <- eval(parse(text=ix_code[[i]][ilen-2]))
      iss <- eval(parse(text=ix_code[[i]][ilen]))
      if(irel<0 || irel>1) stop("Invalid relative position in INTER statement ",paste(iev,collapse=":"),call.=F)
      if((ix_code[[i]][ilen-1] %in% c("lp","n","hp"))==F) stop("Unknown range in INTER statement ",paste(iev,collapse=":"),call.=F)
      if(iss<=0) stop("Invalid virtual cases in INTER statement for ",paste(iev,collapse=":"),call.=F) 
      iran <- ix_code[[i]][ilen-1]
      iyint <- ifelse(iran=="lp",-1.5,ifelse(iran=="n",-0.5,0.5))+irel
      intess[i] <- iss
      intnam[i] <- paste(ix_code[[i]][2:(ilen-3)],collapse=":")
      xtilInt[i] <- prod(xtil[iev])
      betaInt.hat[i] <- log((1.5+iyint)/(1.5-iyint))/(log(5)*nuFun(xtilInt[i]))-sum(beta.hat[iev])
      }
    names(betaInt.hat) <- names(intess) <- names(xtilInt) <- intnam
    }
  q.hat <- c(q.hat,intess)
  xtil <- c(xtil,xtilInt)
  auxtau <- ifelse(tau_code[[1]][3]=="lp",-1.5,ifelse(tau_code[[1]][3]=="n",-0.5,0.5))+taupos
  myfun <- function(z) {
    qbeta(0.25,0.5*exp(z),0.5*exp(z))-(auxtau+1.5)/3
    }
  tau.hat <- exp(BBsolve(1,myfun,quiet=T)$par)                            
  BHAT <- c(beta.hat,betaInt.hat)
  n <- length(beta.hat)
  m <- length(betaInt.hat)
  k <- n+m
  G <- c()                
  for(i in 1:n) {                                          
    iGaux <- matrix(0,nrow=q.hat[names(beta.hat)[i]],ncol=k) 
    for(j in 1:nrow(iGaux)) {
      iGaux[j,i] <- xtil[names(beta.hat)[i]]
      }
    G <- rbind(G,iGaux)
    }
  if(m>0) {
    for(i in 1:m) {
      iauxnam <- names(betaInt.hat)[i]
      iauxev <- strsplit(iauxnam,":")[[1]]
      iGaux <- matrix(0,nrow=q.hat[iauxnam],ncol=k)                                          
      for(j in 1:nrow(iGaux)) {
        iGaux[j,n+i] <- xtil[iauxnam]        
        for(w in 1:length(iauxev)) {
          iGaux[j,which(names(BHAT)==iauxev[w])] <- xtil[iauxev[w]]             
          }
        }
      G <- rbind(G,iGaux)
      }  
    }
  auxG <- G
  bootEstim <- matrix(nrow=nrep,ncol=k+1)
  cat("  |",rep(" ",50),"| 0%",sep="")
  flush.console() 
  for(i in 1:nrep) {
    ysam <- c()
    for(j in 1:nrow(auxG)) {                           
      auxG[j,] <- nuFun(G[j,])
      auxd <- plogis(log(5)*t(auxG[j,])%*%BHAT)
      ysam[j] <- rbeta(1,auxd*tau.hat,(1-auxd)*tau.hat)
      }
    mod0 <- betareg(ysam~auxG-1)
    bootEstim[i,] <- c(coef(mod0)[1:k]/log(5),log(coef(mod0)[k+1]))
    indStar <- round(i*50/nrep)
    cat('\r',"  |",rep("*",indStar),rep(" ",50-indStar),"| ",indStar*2,"%",sep="")
    flush.console() 
    }
  cat("\n")
  mu <- c(BHAT,'(log-precision)'=log(tau.hat))
  S <- cov(bootEstim)
  rownames(S) <- colnames(S) <- names(mu)
  res <- list(scales=scales,mean=round(mu,digits),vcov=round(S,digits))
  class(res) <- "mbr"
  res
  }

# print method for class 'mbr'
print.mbr <- function(x,...) {
  cat("Bayesian Elicitation of Dependence Relationships among Medical Variables","\n")
  cat("Response: ",names(x$scales)[1],sep="","\n","\n")
  print(x[c("mean","vcov")])
  }

# summary method for class 'mbr'
summary.mbr <- function(object,...) {
  mu <- object$mean
  S <- object$vcov
  q99 <- -qnorm(0.005)
  q95 <- -qnorm(0.025)
  res <- matrix(nrow=length(mu),ncol=5)
  nomi <- names(mu)
  for(i in 1:length(nomi)) {
    res[i,] <- mu[nomi[i]]+c(-q99,-q95,0,q95,q99)*sqrt(S[nomi[i],nomi[i]])
    }
  rownames(res) <- nomi
  colnames(res) <- c("0.5%","2.5%","50%","97.5%","99.5%")
  res["(log-precision)",] <- exp(res["(log-precision)",])
  rownames(res) <- c(setdiff(rownames(res),"(log-precision)"),"(precision)")
  res
  }

# plot method for class 'mbr'
plot.mbr <- function(x,...) {
  mu <- x$mean
  S <- x$vcov
  k <- length(mu)-1
  nomi <- names(mu)[1:k]
  parLen <- round((k+1)/3+0.49)
  par(mfrow=c(round((k+1)/parLen+0.49),parLen))
  for(i in 1:k) {
    if(identical(nomi,paste("b",1:k,sep=""))) {
      curve(dnorm(x,mu[i],sqrt(S[i,i])),xlim=mu[i]+c(-3,3)*sqrt(S[i,i]),type="n",las=1,
        xlab="",ylab="",main=substitute(paste(beta[ind],sep=""),list(ind=i)),cex.main=1.4)
        } else {
      curve(dnorm(x,mu[i],sqrt(S[i,i])),xlim=mu[i]+c(-3,3)*sqrt(S[i,i]),type="n",las=1,
        xlab="",ylab="",main=nomi[i],cex.main=1)
      }
    grid(col="grey70",lty=2)
    curve(dnorm(x,mu[i],sqrt(S[i,i])),add=T)
    abline(v=mu[i],lty=2)
    box()
    }
  tauLim <- qlnorm(c(.005,.995),mu[k+1],sqrt(S[k+1,k+1]))
  curve(dlnorm(x,mu[k+1],sqrt(S[k+1,k+1])),xlim=tauLim,type="n",las=1,
    xlab="",ylab="",main=expression(paste(tau)),cex.main=1.4)
  grid(col="grey70",lty=2)
  curve(dlnorm(x,mu[k+1],sqrt(S[k+1,k+1])),add=T)
  abline(v=exp(mu[k+1]-S[k+1,k+1]),lty=2)
  box()
  par(mfrow=c(1,1))   
  }

# compute predictive distributions
predictive <- function(x,cnfg,nrep=500,title=NULL) {
  if(class(x)!="mbr") stop("The first argument must be an object of class 'mbr'",call.=F)
  if(!identical(round(nrep),nrep) | nrep<500) stop("Argument 'nrep' must be an integer number no less than 500",call.=F)
  scales.list <- x$scales
  yscale <- scales.list[[1]]
  nomiInt <- names(x$mean)[grepl(":",names(x$mean))]
  nomiFact <- setdiff(names(x$mean),nomiInt)
  if(length(cnfg)!=(length(scales.list)-1)) stop("Invalid configuration",call.=F)
  parCnfg <- c()              
  for(i in 1:length(cnfg)) {
    iscal <- scales.list[[i+1]]
    if(identical(iscal,NA)) {
      parCnfg[i] <- cnfg[i]
      } else {
      parCnfg[i] <- rescale(cnfg[i],iscal)
      }
    names(parCnfg)[i] <- names(scales.list)[i+1]
    }                                             
  if(length(nomiInt)>0) {
    for(i in 1:length(nomiInt)) {
      iaux <- strsplit(nomiInt[i],":")[[1]]
      parCnfg <- c(parCnfg,prod(parCnfg[iaux]))
      }
    }   
  pmean <- x$mean
  S <- x$vcov
  k <- length(pmean)-1
  ypred <- ytil <- c()
  cat("  |",rep(" ",50),"| 0%",sep="")    
  flush.console() 
  for(i in 1:nrep) {
    theta <- pmean+c(rnorm(k+1)%*%chol(S))
    relmu <- plogis(log(5)*nuFun(parCnfg)%*%theta[1:k])
    taU <- exp(theta[k+1])
    ytil[i] <- -1.5+3*rbeta(1,relmu*taU,(1-relmu)*taU)
    ypred[i] <- invRescale(ytil[i],yscale)
    indStar <- round(i*50/nrep)
    cat('\r',"  |",rep("*",indStar),rep(" ",50-indStar),"| ",indStar*2,"%",sep="")
    flush.console() 
    }
  cat("\n")
  if(is.null(title)) {
    cfglab <- paste(paste(nomiFact,"=",cnfg,sep=""),collapse=", ")
    } else {
    cfglab <- title
    }
  xquart <- quantile(ypred,prob=c(0.25,0.5,0.75))
  ydens <- density(ypred)
  limy <- c(0,max(ydens$y))
  limx <- c(yscale[1],yscale[4])          
  plot(ydens,cex.axis=1.1,cex.lab=1.1,type="n",cex.main=1.2,main=cfglab,cex.main=1.2,las=1,
    xlab=paste("N = ",nrep," (bw: ",round(ydens$bw,6),
    ")    Median: ",round(xquart[2],6),"    IQR: ",round(xquart[3]-xquart[1],6),sep=""),
    ylab="",yaxs="i",xaxs="i",ylim=limy,xlim=limx)
  grid(col="grey70",nx=0,ny=NULL,lty=2)
  lines(ydens,lwd=2)
  abline(v=yscale[2:5],lty=1)
  abline(v=c((yscale[1]+yscale[2])/2,(yscale[2]+yscale[3])/2,(yscale[3]+yscale[4])/2),lty=2)
  box()
  }

# plot the marginal conditional expected value function      
eValFun.plot <- function(B=1,yscale=NULL,xscale=NULL) {
  if(length(B)!=1 || !is.numeric(B)) stop("Argument 'B' must be a numerical value of length 1",call.=F)
  y.restrict <- x.restrict <- "no"
  xstd <- ystd <- F
  if(is.null(yscale)) {
    yscale <- c(-1.5,-0.5,0.5,1.5)
    ystd <- T
    } else {
    if(length(yscale)!=4 | (identical(sort(yscale),yscale)==F&identical(sort(yscale,T),yscale)==F)) {
      stop("Invalid scale for the response variable",call.=F) 
      }
    if(yscale[1]==yscale[2]) y.restrict <- "hyper"
    if(yscale[3]==yscale[4]) y.restrict <- "hypo"
    }  
  if(is.null(xscale)) {
    xscale <- c(-1.5,-0.5,0.5,1.5)
    xstd <- T
    } else {
    if(length(xscale)!=4 | (identical(sort(xscale),xscale)==F&identical(sort(xscale,T),xscale)==F)) {
      stop("Invalid scale for the explanatory variable",call.=F)
      }
    if(xscale[1]==xscale[2]) x.restrict <- "hyper"
    if(xscale[3]==xscale[4]) x.restrict <- "hypo"
    }                                     
  if(y.restrict=="hyper") yscale[1] <- NA
  if(y.restrict=="hypo") yscale[4] <- NA
  if(x.restrict=="hyper") xscale[1] <- NA
  if(x.restrict=="hypo") xscale[4] <- NA
  par(mar=c(5.1,6.2,4.1,2.8))
  plot(0,type="n",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),yaxs="i",xaxs="i",yaxt="n",xaxt="n",
    xlab="",ylab="",mgp=c(2,0.75,0),main=substitute(paste(beta[i]," = ",ind,sep=""),list(ind=B)),cex.main=1.4)
  ynewlab <- c(yscale[1],(yscale[1]+yscale[2])/2,yscale[2],(yscale[2]+yscale[3])/2,yscale[3],(yscale[3]+yscale[4])/2,yscale[4])
  axis(2,labels=ynewlab,at=c(-1.5,-1,-0.5,0,0.5,1,1.5),cex.axis=1.2,cex.lab=1.3,tck=0)
  xnewlab <- c(xscale[1],(xscale[1]+xscale[2])/2,xscale[2],(xscale[2]+xscale[3])/2,xscale[3],(xscale[3]+xscale[4])/2,xscale[4])
  axis(1,labels=xnewlab,at=c(-1.5,-1,-0.5,0,0.5,1,1.5),cex.axis=1.2,cex.lab=1.3,tck=0)
  abline(h=c(-0.5,0.5),v=c(-0.5,0.5),lty=1,col="grey50")
  abline(h=c(-1,0,1),v=c(-1,0,1),lty=2,col="grey70")
  if(ystd==T & xstd==T) {  
    title(ylab=expression(paste("E [ ",tilde(Y)," ]",sep="")),mgp=c(4.5,2,0),cex.lab=1.4)
    title(xlab=expression(paste(tilde(X)[i])),mgp=c(4.2,2,0),cex.lab=1.4)
    } else {
    title(ylab="E [ Y ]",mgp=c(4.5,2,0),cex.lab=1.4)    
    title(xlab=expression(paste(X[i])),mgp=c(4.2,2,0),cex.lab=1.4) 
    }
  xaux <- seq(-1.5,1.5,length=1000)
  yaux <- -1.5+3/(1+exp(-log(5)*nuFun(xaux)*B))
  if(y.restrict=="hyper") yaux[which(yaux< -0.5)] <- -0.5
  if(y.restrict=="hypo") yaux[which(yaux>0.5)] <- 0.5
  if(x.restrict=="hyper") xaux[which(xaux< -0.5)] <- NA
  if(x.restrict=="hypo") xaux[which(xaux>0.5)] <- NA  
  lines(xaux,yaux,lwd=2)  
  text(0,-1.84,expression(italic("n-range")),xpd=T,cex=1.2)
  if(x.restrict!="hyper") text(-1,-1.84,expression(italic("lp-range")),xpd=T,cex=1.2)
  if(x.restrict!="hypo") text(1,-1.84,expression(italic("hp-range")),xpd=T,cex=1.2)
  text(-1.86,0,expression(italic("n-range")),xpd=T,cex=1.2,srt=90)
  if(y.restrict!="hyper") text(-1.86,-1,expression(italic("lp-range")),xpd=T,cex=1.2,srt=90)
  if(y.restrict!="hypo") text(-1.86,1,expression(italic("hp-range")),xpd=T,cex=1.2,srt=90)
  box()
  }