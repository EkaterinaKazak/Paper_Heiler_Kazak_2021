##############################################################
#Heiler and Kazak (2021): Re-analysis of Connors et al. (1996)
#using overlap robust inference: #############################
##############################################################
library(sandwich)
library(caret)
library(foreach)
library(Matrix)

#1) Extract and relabel data:
{rhc <- read.csv("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.csv")[,-c(1,63)]
  rhc[,"death"] <- as.numeric(rhc[,"death"]) -1
  attach(rhc)
  weight_dum = (rhc[,"wtkilo1"] == 0)
  rhc2 = cbind(rhc,weight_dum)
  n <- dim(rhc2)[1]
}
{rhc3 <- rhc2[,c("death","weight_dum","wtkilo1","age","sex","race","edu","income","ninsclas",
                 "cat1","das2d3pc","dnr1","ca","surv2md1","aps1","scoma1",
                 "temp1","meanbp1","resp1","hrt1","pafi1","paco21","ph1","wblc1", 
                 "hema1","sod1","pot1","crea1","bili1","alb1","resp","card","neuro", 
                 "gastr","renal","meta","hema","seps","trauma","ortho","cardiohx", 
                 "chfhx","dementhx","psychhx","chrpulhx","renalhx","liverhx", 
                 "gibledhx","malighx","immunhx","transhx","amihx","swang1")]}

#2.1) Hirano & Imbens (2001) models for pscore and potential outcomes:
{m1 = formula(swang1 ~ weight_dum  + I(1-weight_dum):wtkilo1 +
                age + sex + race + edu + income + ninsclas + 
                cat1 + das2d3pc + dnr1 + ca + surv2md1 + aps1 + scoma1 + 
                temp1 + meanbp1 + resp1 + hrt1 + pafi1 + paco21 + ph1 + wblc1 + 
                hema1 + sod1 + pot1 + crea1 + bili1 + alb1 + resp + card + neuro + 
                gastr + renal + meta + hema + seps + trauma + ortho + cardiohx + 
                chfhx + dementhx + psychhx + chrpulhx + renalhx + liverhx + 
                gibledhx + malighx + immunhx + transhx + amihx)
  
  mr1 = formula(I(1-death) ~ weight_dum  + I(1-weight_dum):wtkilo1 +
                  age + sex + race + edu + income + ninsclas + 
                  cat1 + das2d3pc + dnr1 + ca + surv2md1 + aps1 + scoma1 + 
                  temp1 + meanbp1 + resp1 + hrt1 + pafi1 + paco21 + ph1 + wblc1 + 
                  hema1 + sod1 + pot1 + crea1 + bili1 + alb1 + resp + card + neuro + 
                  gastr + renal + meta + hema + seps + trauma + ortho + cardiohx + 
                  chfhx + dementhx + psychhx + chrpulhx + renalhx + liverhx + 
                  gibledhx + malighx + immunhx + transhx + amihx) 
}
#2.2) Declare functions for estimation and inference:
#Inputs: 
{
#q = coarseness parameter for grid as in Bickel and Sakov (2008)
#B = Number of Bootstrap Iterations 
#MC = Number of Monte Carlo simulations 
}
#Functions:
{
muav_drN <- function(y,d,p,mu1,mu0){ 
  DR <- vector(mode = "list", length = 2)
  DR$estimate <- mean(d*(y-mu1)/p - (1-d)*(y-mu0)/(1-p) + mu1 - mu0)
  DR$variance <- mean(d*(y - mu1)^2/p^2) + mean((1-d)*(y - mu0)^2/(1-p)^2) +
    mean((mu1 - mu0 - DR$estimate)^2)
  return(DR) 
}
m_BS_DR1reest <-function(rhc3,q,B,MC){
  #global: mr1,m1
  n=dim(rhc3)[1]
  mgrid=ceiling(n*q^seq(0,40,by=1)); mgrid=mgrid[mgrid>(log(n)/log(66))^2*66];mgrid=unique(mgrid);
  lm=length(mgrid);
  mboot <-  foreach (b = 1:MC, .packages = c("sandwich","caret"),
                     .export = c("mr1","m1","muav_drN"), .combine='c') %dopar% {
    tb=matrix(0,B,lm)
    for (j in 1:lm){
      m=mgrid[j]
      for (bb in 1:B){
        J = ceiling(runif(m)*n)
        { #redefine data for reestimation:
        rhc_bo <- (rhc3[J,]) 
        ui <- dummyVars(mr1,data=rhc_bo,fullRank=TRUE)
        rhc_b <- cbind(1-rhc_bo[,"death"],data.frame(predict(ui,newdata = rhc_bo)))
        colnames(rhc_b)[1] <- "death"
        rhc_b1 <- (rhc_b[rhc_bo[,"swang1"]=="RHC",])
        rhc_b0 <- (rhc_b[rhc_bo[,"swang1"]!="RHC",]) 
        fitb1 <- lm(death~.,data= rhc_b1)
        fitb0 <- lm(death~.,data=rhc_b0)
        mu1 <- predict(fitb1,rhc_b,type="response") 
        mu0 <- predict(fitb0,rhc_b,type="response")
        ui3 <- dummyVars(m1,data=rhc_bo,fullRank=TRUE)
        rhc_b3 <- cbind(rhc_bo[,"swang1"],data.frame(predict(ui3,newdata = rhc_bo)))
        colnames(rhc_b3)[1] <- "swang1"
        glmQuadb <- glm(swang1~., family  = binomial(link = "probit"),data  = rhc_b3) 
        p <- predict(glmQuadb,type="response")
        }
        DR <- muav_drN(rhc_bo[,1],(rhc_bo[,"swang1"]=="RHC"),p,mu1,mu0)
        tb[bb,j] <- sqrt(m)*DR$estimate/sqrt(DR$variance)
      }
      
    }
    rho=matrix(0,lm-1,1)
    for (i in 1:(lm-1)){
      out<- ks.test(tb[,i+1],tb[,i]) 
      rho[i]=out$statistic
    }
    iopt <- apply(rho,2,which.min)
    mgrid[iopt] 
  }
  return(mboot)
}
muav_ipwNreest <- function(y,d,p,S){ 
  IPW <- vector(mode = "list", length = 2)
  m <- length(y)
  mu1 <- mean(y*d/p)/mean(d/p)
  mu0 <- mean(y*(1-d)/(1-p))/mean((1-d)/(1-p))
  k1 <- d*(y-mu1)/p
  k0 <- (1-d)*(y-mu0)/(1-p)
  {
  ###############################################################################################
  #workaround against numerical error for inverse that has small (but pos) eigenvalues############
  #Use scale-normalization for Thik. correction/remove intercept part of score:###################
  #corresponds to ridge with lambda = O(1), or adding O(1/n) on E[SiSi'] -> difference is Op(1/n)
  SS <- (t(S)%*%(S))
  if (rankMatrix(SS)[1] != dim(SS)[1]){
    
    SS <- (t(S)%*%(S)) + diag(c(0,apply(S[,-c(1)],2,var)))  
  } 
  ################################################################################################
  }
  SSi <-  chol2inv(chol(SS)) 
  b1 = colSums(matrix(as.numeric(k1),nrow = m, ncol=dim(S)[2])*S)%*%SSi
  b0 = colSums(matrix(as.numeric(k0),nrow = m, ncol=dim(S)[2])*S)%*%SSi
  kap1 = mean((k1 - S%*%t(b1))^2)
  kap0 = mean((k0 - S%*%t(b0))^2)
  omega =b1%*%SS%*%t(b0)/m 
  
  IPW$estimate = mu1 - mu0
  IPW$variance = as.numeric(kap1 + kap0 - 2*omega)
  return(IPW) 
}
m_BS_IPW2reest <-function(rhc3,q,B,MC){
  n=dim(rhc3)[1]
  mgrid=ceiling(n*q^seq(0,40,by=1)); mgrid=mgrid[mgrid>(log(n)/log(66))^2*66];mgrid=unique(mgrid);
  lm=length(mgrid); 
  mboot<-  foreach (b = 1:MC, .packages = c("sandwich","caret","Matrix"), 
                    .export = c("m1","muav_ipwNreest"), .combine='c') %dopar% {
    tb=matrix(0,B,lm)
    for (j in 1:lm){
      m=mgrid[j]
      for (bb in 1:B){ 
        J = ceiling(runif(m)*n)
        {
        ind <- sample.int(n,m,replace=TRUE)
        rhc_bo <- (rhc3[J,]) 
        db = as.matrix((rhc_bo[,"swang1"] == "RHC"))
        ui3 <- dummyVars(m1,data=rhc_bo,fullRank=TRUE)
        rhc_b3 <- cbind(db,data.frame(predict(ui3,newdata = rhc_bo)))
        colnames(rhc_b3)[1] <- "swang1"
        glmQuadb <- glm(swang1~., family  = binomial(link = "probit"),data  = rhc_b3) 
        pscoreb <- predict(glmQuadb,type="response")
        yb = 1-as.matrix(rhc_bo[,"death"])
        Sb <- estfun(glmQuadb,length="m")
        }
        IPW <- muav_ipwNreest(yb,db,pscoreb,Sb)
        tb[bb,j] <- sqrt(m)*IPW$estimate/sqrt(IPW$variance)
      }
    }
    rho=matrix(0,lm-1,1)
    for (i in 1:(lm-1)){
      out<- ks.test(tb[,i+1],tb[,i]) 
      rho[i]=out$statistic
    }
    iopt <- apply(rho,2,which.min)
    mgrid[iopt] 
    }
  return(mboot)
}
}
#2.3) Estimates:
{
glmFull <- glm(m1,  family  = binomial(link = "probit"), data  = rhc2)
pscore <- vector(mode = "list", length = 1) 
pscore$m1 <- predict(glmFull,type="response")
S <- estfun(glmFull,length="n")
fit11 <- lm(mr1,data=rhc2[(rhc2[,"swang1"]=="RHC"),])
fit10 <- lm(mr1,data=rhc2[(rhc2[,"swang1"]!="RHC"),])
pot11 <- predict(fit11,rhc2,type="response")
pot10 <- predict(fit10,rhc2,type="response")
d = as.matrix((rhc2[,"swang1"] == "RHC"))
y = 1-as.matrix(rhc2[,"death"])

avipw <- muav_ipwNreest(y,d,pscore$m1,S)
avdr <- muav_drN(y,d,pscore$m1,pot11,pot10)
ateIPW <- avipw$estimate  
ateDR <- avdr$estimate  
}
#2.4) Creates Figure 7.1: Propensity Score Distributions:
{bk <- seq(0,1,l=16)
#setEPS()
#postscript("g_RHC_pscoreHIST2.eps")
par(mfrow=c(1,2),family = 'Latin Modern Math',pty="s")
hist(pscore$m1[(rhc2[,"swang1"]=="RHC")], xlim=c(0,1),breaks=bk,font.main=1,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     ylim=c(0,600),freq=TRUE,ylab = "", xlab="",main = "Treatment, n = 2184",col="gainsboro")
hist(pscore$m1[(rhc2[,"swang1"]!="RHC")], xlim=c(0,1),breaks=bk,font.main=1,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     ylim=c(0,600),freq=TRUE,ylab = "", xlab="",main = "Control, n = 3551",col="gainsboro")
#dev.off()
}
#Print point estimates for average treatment effect:
cbind(ateIPW,ateDR)

#3) Determine optimal m-choice: 
{
#HERE: Monte Carlo commented out -> read .csv for results
#cl <- parallel::makeCluster(4)
#doParallel::registerDoParallel(cl)
#msDR <- m_BS_DR1reest(rhc3,0.75,1000,250)
#msIPW <- m_BS_IPW2reest(rhc3,0.75,1000,250)
#parallel::stopCluster(cl)

ms <-  vector(mode = "list", length = 2) 
ms$dr <- read.csv("drMCHOICE.csv",header=TRUE)[,2]
ms$ipw <- read.csv("ipwMCHOICE.csv",header=TRUE)[,2] 
} 

#4.1) Inference using asymptotic normality, percentile bootstrap,
#and overlap robust moon using the selected optimal m:
#NOTE: Bootstrap loops require some computation time!
{B <- 5000 #Number of Bootstrap Iterations
ate = matrix(0,B,2)
sig = matrix(0,B,2)
ateBS = matrix(0,B,2)
sigBS = matrix(0,B,2)

####moon######################
#IPW
m = round(mean(ms$ipw))
for (b in (1:B)) {
  m = round(mean(ms$ipw))
  ind <- sample.int(n,m,replace=TRUE)
  rhc_bo <- (rhc3[ind,]) 
  
  ui3 <- dummyVars(m1,data=rhc_bo,fullRank=TRUE)
  rhc_b3 <- cbind(as.numeric(rhc_bo[,"swang1"])-1,data.frame(predict(ui3,newdata = rhc_bo)))
  colnames(rhc_b3)[1] <- "swang1"
  glmQuadb <- glm(swang1~., family  = binomial(link = "probit"),data  = rhc_b3) 
  pscoreb <- predict(glmQuadb,type="response")
  
  db = as.matrix((rhc_bo[,"swang1"] == "RHC"))
  yb = 1-as.matrix(rhc_bo[,"death"])
  Si <- estfun(glmQuadb,length="m")

  IPW_b <- muav_ipwNreest(yb,db,pscoreb,Si)
  ate[b,1] <- IPW_b$estimate
  sig[b,1] <- sqrt(IPW_b$variance/m)
  print(b/B)
}
#DR
m = round(mean(ms$dr))
for (b in (1:B)) {
    
    ind <- sample.int(n,m,replace=TRUE)
    rhc_bo <- (rhc3[ind,]) 
    
    ui <- dummyVars(mr1,data=rhc_bo,fullRank=TRUE)
    rhc_b <- cbind(1-rhc_bo[,"death"],data.frame(predict(ui,newdata = rhc_bo)))
    colnames(rhc_b)[1] <- "death"
    
    rhc_b1 <- (rhc_b[rhc_bo[,"swang1"]=="RHC",]) 
    rhc_b0 <- (rhc_b[rhc_bo[,"swang1"]!="RHC",]) 
    fitb1 <- lm(death~.,data= rhc_b1)
    fitb0 <- lm(death~.,data=rhc_b0)
    potb1 <- predict(fitb1,rhc_b,type="response") 
    potb0 <- predict(fitb0,rhc_b,type="response")
    
    ui3 <- dummyVars(m1,data=rhc_bo,fullRank=TRUE)
    rhc_b3 <- cbind(rhc_bo[,"swang1"],data.frame(predict(ui3,newdata = rhc_bo)))
    colnames(rhc_b3)[1] <- "swang1"
    glmQuadb <- glm(swang1~., family  = binomial(link = "probit"),data  = rhc_b3) 
    pscoreb <- predict(glmQuadb,type="response")
    
    db = as.matrix((rhc_bo[,"swang1"] == "RHC"))
    yb = 1-as.matrix(rhc_bo[,"death"])
    
    DR_b <- muav_drN(yb,db,pscoreb,potb1,potb0)
    ate[b,2] <- DR_b$estimate 
    sig[b,2] <- sqrt(DR_b$variance/m) 
    print(b/B)
  }
####nonparametric Bootstrap####
m <- n
for (b in (1:B)) {
    
    ind <- sample.int(n,m,replace=TRUE)
    rhc_bo <- (rhc3[ind,]) 
    ui <- dummyVars(mr1,data=rhc_bo,fullRank=TRUE)
    rhc_b <- cbind(1-rhc_bo[,"death"],data.frame(predict(ui,newdata = rhc_bo)))
    colnames(rhc_b)[1] <- "death"
    rhc_b1 <- (rhc_b[rhc_bo[,"swang1"]=="RHC",]) 
    rhc_b0 <- (rhc_b[rhc_bo[,"swang1"]!="RHC",]) 
    fitb1 <- lm(death~.,data= rhc_b1)
    fitb0 <- lm(death~.,data=rhc_b0)
    potb1 <- predict(fitb1,rhc_b,type="response") 
    potb0 <- predict(fitb0,rhc_b,type="response")
    
    ui3 <- dummyVars(m1,data=rhc_bo,fullRank=TRUE)
    rhc_b3 <- cbind(rhc_bo[,"swang1"],data.frame(predict(ui,newdata = rhc_bo)))
    colnames(rhc_b3)[1] <- "swang1"
    glmQuadb <- glm(swang1~., family  = binomial(link = "probit"),data  = rhc_b3) 
    pscoreb <- predict(glmQuadb,type="response")
    
    db = as.matrix((rhc_bo[,"swang1"] == "RHC"))
    yb = 1-as.matrix(rhc_bo[,"death"])
    Si <- estfun(glmQuadb,length="m")

    IPW_b <- muav_ipwNreest(yb,db,pscoreb,Si)
    DR_b <- muav_drN(yb,db,pscoreb,potb1,potb0)
    ateBS[b,2] <- DR_b$estimate 
    sigBS[b,2] <- sqrt(DR_b$variance/m)
    ateBS[b,1] <- IPW_b$estimate
    sigBS[b,1] <- sqrt(IPW_b$variance/m)
    print(b/B)
}
}
#4.2) Calculate asymptotic p-values (Table 7.1):
{
#I) asymptotic Normality, 
#II) bootstrap (percentile),
#III) moon(robust):
pm = matrix(0,3,2)
tIPW <- (sqrt(n)*(avipw$estimate/sqrt(avipw$variance)))
tDR <- (sqrt(n)*(avdr$estimate/sqrt(avdr$variance)))
#asymtptotic normality
pm[1,1] <- 2*(1-pnorm(abs(tIPW)))
pm[1,2] <- 2*(1-pnorm(abs(tDR)))

#conventional bootstrapp percentile+data-centered 
tBS_IPW <- (ateBS[,1]-avipw$estimate)/sigBS[,1]
pm[2,1] <- mean(abs(tBS_IPW) > abs(tIPW)) 
tBS_DR <- (ateBS[,2]-avdr$estimate)/sigBS[,2]
pm[2,2] <- mean(abs(tBS_DR) > abs(tDR))

#moon overlap robust p-values
tmoon_IPW <- (ate[,1]-avipw$estimate)/sig[,1]
pm[3,1] <- mean(abs(tmoon_IPW) > abs(tIPW)) 
tmoon_DR <- (ate[,2]-avdr$estimate)/sig[,2]
pm[3,2] <- mean(abs(tmoon_DR) > abs(tDR))
} 
#Print estimates + p-values:
rbind(cbind(ateIPW,ateDR),pm)
#Additional Values for Table 7.1: Symmetric 95% and 99% confidence intervals:
{
# CI1_IPWnormal <- avipw$estimate + sqrt(avipw$variance/n)*cbind(qnorm(0.005),qnorm(0.995))
# CI1_IPWboot <- avipw$estimate + sqrt(avipw$variance/n)*c(-1,1)*quantile(abs(tBS_IPW),0.99)
# CI1_IPWmoon <- avipw$estimate + sqrt(avipw$variance/n)*c(-1,1)*quantile(abs(tmoon_IPW),0.99)
# 
# CI5_IPWnormal <- avipw$estimate + sqrt(avipw$variance/n)*cbind(qnorm(0.025),qnorm(0.975))
# CI5_IPWboot <- avipw$estimate + sqrt(avipw$variance/n)*c(-1,1)*quantile(abs(tBS_IPW),0.95)
# CI5_IPWmoon <- avipw$estimate + sqrt(avipw$variance/n)*c(-1,1)*quantile(abs(tmoon_IPW),0.95)
# 
# CI1_DRnormal <- avdr$estimate + sqrt(avdr$variance/n)*cbind(qnorm(0.005),qnorm(0.995))
# CI1_DRboot <- avdr$estimate + sqrt(avdr$variance/n)*c(-1,1)*quantile(abs(tBS_DR),0.99)
# CI1_DRmoon <- avdr$estimate + sqrt(avdr$variance/n)*c(-1,1)*quantile(abs(tmoon_DR),0.99)
# 
# CI5_DRnormal <- avdr$estimate + sqrt(avdr$variance/n)*cbind(qnorm(0.025),qnorm(0.975))
# CI5_DRboot <- avdr$estimate + sqrt(avdr$variance/n)*c(-1,1)*quantile(abs(tBS_DR),0.95)
# CI5_DRmoon <- avdr$estimate + sqrt(avdr$variance/n)*c(-1,1)*quantile(abs(tmoon_DR),0.95)
}






