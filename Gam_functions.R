#-----------------------------------------------------------------------------------------------
# Script: gam for multivariate data: functions
# Project: METACOM 1
# Date: 25/03/2021
#-----------------------------------------------------------------------------------------------

# Latent variables
gll.hat.fun <- function(sp, family = "binomial", max.lv = 3){
  # determine number of latent variables; max.lv = maximum number of variables
  lv <- as.list(1:max.lv)
  for(i in 1:length(lv)){
    lv[[i]] <- gllvm(sp,family = family, num.lv = i)
  }
  # select model
  aic <- unlist(lapply(lv, AIC))
  num.lv <- which(aic == min(aic))
  gll.hat <- as.data.frame(predict.gllvm(lv[[num.lv]], type = "response"))
  return(gll.hat)
}

#----------------------------------------------------------------------------------------------------------------------------------------

# r.squared using GAM for multivariate data (non-parallel computing. Not used in the final code)
gam.multi <- function(sp, pred, family = "quasibinomial", select = FALSE, method  = "REML",...){
  #This function obtains the r squared of GAMs for multivariate data
  #sp=species matrix
  #pred=predictor matrix. For example, environmental matrix
  #family=distribution and link in fitting GAMs,  for the GAM function in mgcv
  #select=if FALSE, it does not remove terms in the model,  for the GAM function in mgcv
  #method=smoothing parameter estimation method, for the GAM function in mgcv
  r.sq.vec <- rep(0, ncol(sp)) #Empty vector for r squared of each species
  fm.pred <- paste("s(", colnames(pred), ")", sep = "", collapse = " + ") #Term structure for GAM formula
  for(j in 1:ncol(sp)){ #loop to perform GAM for each species. 
    fm <- as.formula(paste("sp[,",j,"] ~ ", fm.pred, sep = "")) #Formula for each species model
    r.sq.vec[j] <- summary(gam(fm, data = pred,
                               family = family, select = select, method = method, ...))$r.sq #GAM for each species and extraction of r2
  } 
  r.sq <- sum(r.sq.vec * apply(sp, 2, var)) / sum(apply(sp, 2, var)) #multivariate r2
  return(r.sq)
}

#--------------------------------------------------------------------------------------------------------------------------------------

# r.squared using GAM for multivariate data (parallel computing, used in the final code for faster analysing)
gam.multi.par <- function(sp, pred, family = "quasibinomial", select = FALSE, method  = "REML",...){
  #This function obtains the r squared of GAMs for multivariate data
  #sp=species matrix
  #pred=predictor matrix. For example, environmental matrix
  #family=distribution and link in fitting GAMs,  for the GAM function in mgcv
  #select=if FALSE, it does not remove terms in the model,  for the GAM function in mgcv
  #method=smoothing parameter estimation method, for the GAM function in mgcv
  fm.pred <- paste("s(", colnames(pred), ")", sep = "", collapse = " + ") #Term structure for GAM formula
  cores <- detectCores() #Number of cores in computer
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl) #Register number of cores to be used
  r.sq.vec <- tryCatch(foreach(j = 1:ncol(sp), .combine = c) %dopar% {#loop to perform GAM for each species.
    fm <- as.formula(paste("sp[,",j,"] ~ ", fm.pred, sep = "")) #Formula for each species model
    summary(mgcv::gam(fm, data = cbind(sp[,j], pred),
                      family = family, select = select, method = method))$r.sq #GAM for each species and extraction of r2
  },  error=function(e) {return(NA)})
  parallel::stopCluster(cl) #Stop parallel computing
  r.sq <- sum(r.sq.vec * apply(sp, 2, var)) / sum(apply(sp, 2, var)) #multivariate r2
  return(r.sq)
}

#-------------------------------------------------------------------------------------------------------------------------------------

# forward selecion for multi.gam (non-parallel computing)

fw.sel.gam <- function(sp, pred, r2a.thresh, alpha = 0.05, ...){
  #sp=species matrix
  #pred=predictor matrix. For example, environmental matrix
  #r2a.thresh=adjusted r2 for forward selection threshold
  #alpha=significance level for forwarsd selection threshold
  fw.tab <- as.data.frame(matrix(0, nrow = ncol(pred), ncol = 7)) #empty table for relevant values of selection
  colnames(fw.tab) <- c("variable", "order", "R2", "F", "pvalue", "R2cum", "adjR2cum") #table colnames. 
  #variable=name of variable
  #order=position of the variable in the predictor matrix
  #R2=multivariate r2 of GAMs using a given variable
  #F=statistic parameter F
  #pvalue=significance of F statistic
  #R2cum=accumulated r2 with a given variable, and the variables above its row
  #adjR2cum=accumulated adjusted r2 (Ezequiel formula) with a given variable, and the variables above its row
  
  fw.tab[,1] <- colnames(pred)#Predictor names, in Variable column
  fw.tab[,2] <- seq(1, nrow(fw.tab), 1) #number of the column of each variable in the pred matrix, for Order column
  for(i in 1:ncol(pred)){ #loop for calculating GAM multivariate r2 for each variable
    fw.tab[i,3] <- gam.multi(sp = sp, pred = pred[,i, drop = FALSE], ...)
  }
  fw.tab <- fw.tab[order(fw.tab[,3], decreasing = TRUE),] #Order Variables in decreasing R2
  n <- nrow(sp) #number of species
  df2 <- (n - 1 - 1) #degrees of freedom (n + number of predictors - 1)
  fw.tab[,4] <- (fw.tab[,3]*df2) / (1 - fw.tab[,3]) #Calculating F statistic (based on adespatial forwardsel function)
  fw.tab[,5] <- round(pf(fw.tab[,4], 1, df2, lower.tail = FALSE), 4) #Calculating pvalue
  r2cum <- fw.tab[1,3] #r2cum of first variabe as object
  adjr2cum <- 1 - (1 - r2cum) * (n - 1)/(n - 1 - 1) #Ezequiel formula for adjusting r2
  for(i in 1:ncol(pred)){#loop for double-selection criterion (Blanchet)
    fw.tab[i,6] <- r2cum
    fw.tab[i,7] <- adjr2cum
    if(fw.tab[i,5] < alpha & fw.tab[i,7] <= r2a.thresh){
      if(i == ncol(pred)){
        print("all selected")
        break
      } else {
        r2cum <- gam.multi(sp = sp, pred = pred[,fw.tab$variable[1:(i+1)]], ...)
        adjr2cum <- 1 - (1 - r2cum) * (n - 1)/(n - (i+1) - 1)  
      }
    } else {
      if(fw.tab[i,5] > alpha) {
        print(paste("variable with a p.value = ", fw.tab[i,5], " > 0.05", sep = ""))
      }
      if(fw.tab[i,7] > r2a.thresh){
        print(paste("next variable sum an adj.R.squared = ", fw.tab[i,7], " > ",
                    r2a.thresh, sep = ""))
      }
      if(i == 1){
        fw.tab <- fw.tab[0, ]  
      } else {
        fw.tab <- fw.tab[1:(i-1), ]
      }
      break
    }
  }
  return(fw.tab)
}

#-----------------------------------------------------------------------------------------------------------------------------

# forward selecion for multi.gam.par (parallel computing)
fw.sel.gam.par<- function(sp, pred, r2a.thresh, alpha = 0.05, ...){
  #sp=species matrix
  #pred=predictor matrix. For example, environmental matrix
  #r2a.thresh=adjusted r2 for forward selection threshold
  #alpha=significance level for forwarsd selection threshold
  fw.tab <- as.data.frame(matrix(0, nrow = ncol(pred), ncol = 7)) #empty table for relevant values of selection
  colnames(fw.tab) <- c("variable", "order", "R2", "F", "pvalue", "R2cum", "adjR2cum") #table colnames. 
  #variable=name of variable
  #order=position of the variable in the predictor matrix
  #R2=multivariate r2 of GAMs using a given variable
  #F=statistic parameter F
  #pvalue=significance of F statistic
  #R2cum=accumulated r2 with a given variable, and the variables above its row
  #adjR2cum=accumulated adjusted r2 (Ezequiel formula) with a given variable, and the variables above its row
  
  fw.tab[,1] <- colnames(pred)#Predictor names, in Variable column
  fw.tab[,2] <- seq(1, nrow(fw.tab), 1) #number of the column of each variable in the pred matrix, for Order column
  for(i in 1:ncol(pred)){ #loop for calculating GAM multivariate r2 for each variable
    fw.tab[i,3] <- gam.multi(sp = sp, pred = pred[,i, drop = FALSE], ...)
  }
  fw.tab <- fw.tab[order(fw.tab[,3], decreasing = TRUE),] #Order Variables in decreasing R2
  n <- nrow(sp) #number of species
  df2 <- (n - 1 - 1) #degrees of freedom (n + number of predictors - 1)
  fw.tab[,4] <- (fw.tab[,3]*df2) / (1 - fw.tab[,3]) #Calculating F statistic (based on adespatial forwardsel function)
  fw.tab[,5] <- round(pf(fw.tab[,4], 1, df2, lower.tail = FALSE), 4) #Calculating pvalue
  r2cum <- fw.tab[1,3] #r2cum of first variabe as object
  adjr2cum <- 1 - (1 - r2cum) * (n - 1)/(n - 1 - 1) #Ezequiel formula for adjusting r2
  for(i in 1:ncol(pred)){#loop for double-selection criterion (Blanchet)
    fw.tab[i,6] <- r2cum
    fw.tab[i,7] <- adjr2cum
    if(fw.tab[i,5] < alpha & fw.tab[i,7] <= r2a.thresh){
      if(i == ncol(pred)){
        print("all selected")
        break
      } else {
        r2cum <- gam.multi.par(sp = sp, pred = pred[,fw.tab$variable[1:(i+1)]], ...)
        adjr2cum <- 1 - (1 - r2cum) * (n - 1)/(n - (i+1) - 1)  
      }
    } else {
      if(fw.tab[i,5] > alpha) {
        print(paste("variable with a p.value = ", fw.tab[i,5], " > 0.05", sep = ""))
      }
      if(fw.tab[i,7] > r2a.thresh){
        print(paste("next variable sum an adj.R.squared = ", fw.tab[i,7], " > ",
                    r2a.thresh, sep = ""))
      }
      if(i == 1){
        fw.tab <- fw.tab[0, ]  
      } else {
        fw.tab <- fw.tab[1:(i-1), ]
      }
      break
    }
  }
  return(fw.tab)
}

#-------------------------------------------------------------------------------------------------------------------------------------

#INCLUDING TIME!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Spatio-temporal neighbours (listw)
xyz2links <- function(names, coords){
  #names=names of the sites and sampling season for each sample
  #coord=X and Y coordinates and sampling season for each sample
  require(spdep)
  names <- as.data.frame(lapply(names, factor)) #names as factor
  coords <- as.data.frame(lapply(coords, as.numeric)) #coords as numeric
  n.site <- nlevels(names[,1]) #number of sites
  n.date <- nlevels(names[,2]) #number of sampling season
  # spatial links
  coord.list <- split(coords[,1:2], names[,2]) #split spatial coordinates for each sampling season
  nb.spa <- lapply(coord.list, function(x) graph2nb(gabrielneigh(as.matrix(x)), sym = TRUE)) #Gabriel graph of each season
  nb.spa <- nb.spa[which(sapply(nb.spa, function(x) sum(card(x))) != 0)] 
  nb.spa <- lapply(nb.spa, function(x) which(unlist(nb2mat(x)) != 0, arr.ind = TRUE))#links between sites for each sampling season
  pos <- sapply(1:length(coord.list), function(x) as.numeric(rownames(coord.list[[x]])[1])) #Combine all the links in different seasons
  pos <- pos[sapply(coord.list, nrow) > 1]
  nb.spa <- lapply(1:length(nb.spa), function(x) nb.spa[[x]] + pos[x] - 1)
  nb.spa <- do.call("rbind", nb.spa)
  # temporal links
  tim.list <- split(1:nrow(names), names[,1]) #list of samples for each site
  nb.tim <- lapply(tim.list, function(x) cbind(x[-length(x)], x[-1])) #Links between times for each site
  nb.tim <- do.call("rbind", nb.tim)
  nb.tim <- rbind(nb.tim, nb.tim[,c(2,1)])
  # spatiotemporal links
  nb.st <- rbind(nb.spa, nb.tim) #Combine spatial and temporal links
  mat.st <- array(0, rep(nrow(names), 2)) #Links as matrix
  mat.st[nb.st] <- 1
  listw.st <- mat2listw(mat.st, style = "W") #Compute weight matrix
  nb <- list("st" = nb.st, "spa" = nb.spa, "tim" = nb.tim) #Save as list of neighbours
  return(list("listw" = listw.st, "nb" = nb)) #Return list for MEMs
}

#-------------------------------------------------------------------------------------

#MSR with spatio-temporal links
xmsrfun3 <- function(env,names, coords, nrepet = 99){
  #env=environmental matrix
  #nrep=number of random environmental matrices
  list<-xyz2links(names = names,coords = coords)#Obtain spatio-temporal links 
  Xmsr <- msr(x = env, listwORorthobasis = list$listw, nrepet = nrepet, simplify  = FALSE) #MSR
  for(i in seq_along(Xmsr)) {#Give names for each variable of the random matrices for subsequent GAMs
    colnames(Xmsr[[i]]) <- paste("V", 1:ncol(env), sep="")
  }
  return(Xmsr)
}

#-----------------------------------------------------------------------------------------------------------

#Adjusting multivariate r2 of environment with MSR correction
adjR2X1.3<- function(sp, env, names, coords, nrepet=10, family = "quasibinomial", ...){
  #sp=species matrix
  #env=environmental matrix
  #xy=spatio-temporal matrix
  #names=names of the sites and sampling season for each sample
  #coord=X and Y coordinates and sampling season for each sample
  #number of permutations in MSR
  #family=family or link function in GAMS
  r.sq.msr<-rep(0,nrepet)#Empty matrix for number of permutations
  Xmsr <- xmsrfun3(env=env, names=names, coords=coords, nrepet = nrepet)
  Xmsr<-lapply(Xmsr, as.data.frame)
  for(i in 1:length(r.sq.msr)){print(i) #loop 
    r.sq.msr[i]<-gam.multi.par(sp = sp, pred = Xmsr[[i]],select=FALSE, method="REML",family = "quasibinomial")
  }
  R2aX1 <- 1 - (1-R2X1) / (1-mean(r.sq.msr, na.rm = TRUE))
  return(R2aX1)
}


adjR2X1X2.3<-function(sp, env, spa, names,coords, nrepet=999, family = "quasibinomial", ...){
  R2X1X2.msr<-rep(0,nrepet)
  Xmsr <- xmsrfun3(env=env,names=names,coords=coords, nrepet = nrepet)
  Xmsr<-lapply(Xmsr, data.frame)
  WXmsr <- lapply(Xmsr, cbind, spa)
    for(i in 1:length(R2X1X2.msr)){print(i) #loop 
      R2X1X2.msr[i]<-gam.multi.par(sp = sp, pred = WXmsr[[i]],select=FALSE, method="REML",family = "quasibinomial")
  }
  return(R2X1X2.msr)
}

#-------------------------------------------------------------------------------------------------------------------------------------

xmsrfun <- function(env, xy, nrepet = 99){
  xyir <- as.matrix(xy)
  nbgab <- graph2nb(gabrielneigh(xyir), sym = TRUE)
  listw <- nb2listw(nbgab)
  Xmsr <- msr(x = env, listwORorthobasis = listw, nrepet = nrepet, simplify  = FALSE)
  for(i in seq_along(Xmsr)) {
    colnames(Xmsr[[i]]) <- paste("V", 1:ncol(env), sep="")
  }
  return(Xmsr)
}

adjR2X1<- function(sp, env, xy, nrepet=999, family = "quasibinomial", ...){
  r.sq.msr<-rep(0,nrepet)
  Xmsr <- xmsrfun(env=env, xy=xy, nrepet = nrepet)
  Xmsr<-lapply(Xmsr, as.data.frame)
  for(i in 1:length(r.sq.msr)){print(i) #loop 
    r.sq.msr[i]<-gam.multi.par(sp = sp, pred = Xmsr[[i]],select=FALSE, method="REML",family = "quasibinomial")
  }
  R2aX1 <- 1 - (1-R2X1) / (1-mean(r.sq.msr, na.rm = TRUE))
  return(R2aX1)
}

adjR2X1X2<-function(sp, env, spa, xy, nrepet=999, family = "quasibinomial", ...){
  R2X1X2.msr<-rep(0,nrepet)
  Xmsr <- xmsrfun(env=env, xy=xy, nrepet = nrepet)
  Xmsr<-lapply(Xmsr, as.data.frame)
  WXmsr <- lapply(Xmsr, cbind, spa)
  for(i in 1:length(R2X1X2.msr)){print(i) #loop 
    R2X1X2.msr[i]<-gam.multi.par(sp = sp, pred = WXmsr[[i]],select=FALSE, method="REML",family = "quasibinomial")
  }
  return(R2X1X2.msr)
}



