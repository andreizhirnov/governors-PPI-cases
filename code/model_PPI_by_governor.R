## Replication materials for
## Shvetsova, Olga, Andrei Zhirnov, Frank Giannelli, Michael Catalano, and Olivia Catalano. 2021. 
##"Governor's Party, Policies, and COVID-19 Outcomes: Further Evidence of an Effect."  

library(splines)
library(foreign)
library(rjags)
library(mcmcplots)
library(MASS)


rm(list=ls())
options(stringsAsFactors = FALSE)

dtdir <- "./data"
psdir <- "./posterior_samples"
outdir <- "./output"

start.date <- as.Date("2020-3-1")
end.date <- as.Date("2020-11-30")

load(file.path(dtdir,"PM_data.RData"))

# states vector
df <- cs.df[order(cs.df$w_rep,cs.df$state),]
state <- df$state
m_tot <- length(unique(df[,"state"]))
m_dem <- length(unique(df[df$w_rep==0,"state"]))


# population
lnoff <- log(df$pop) - log(100000)

# party
repp <- df$w_rep

# other variables
reg <- df$w_region
regl <- lapply(unique(reg), function(x) {
  as.numeric(reg==x) 
})
regl <- do.call("cbind", regl)
colnames(regl) <- paste0("region.", unique(reg))

othvar <- setdiff(colnames(df), c("state","pop","w_rep","w_region"))
mat <- as.matrix(df[,othvar])
# mat <- cbind(mat, regl)
mat <- scale(mat)
mat <- cbind(int=rep(1,nrow(mat)), mat)
xvar <- colnames(mat)

# compute splines
tdf <- subset(csts.df, date >= start.date & date <= end.date)
tdf$tid <- as.numeric(tdf$date - start.date) + 1
uda <- sort(unique(tdf$tid))
knots <- uda[which(uda %% 14 == 0)]
spl <- ns(uda, knots=knots, intercept=TRUE)
splnam <- paste0("sp.",1:ncol(spl))
colnames(spl) <- splnam
spl_base <- data.frame(tid=uda,spl)

# expand the state-level dataset
ndf <- data.frame(state,repp,mat,lnoff)
ndf <- merge(ndf, tdf, by="state")
ndf <- merge(ndf, spl_base, by="tid")

# generate datasets by governor's party
frm <- as.formula(paste0("ppi ~ ",paste(splnam, collapse="+")))
m <- lm(frm, data=ndf)
coef <- coef(m)[splnam]

## Prepare the data list  
dt <- ndf
dt$sid <-  match(dt$state,state)
A.prep <- unique(dt[c("tid",splnam)])
A <- as.matrix(A.prep[order(A.prep$tid),splnam])
mydata <- list(y = dt$ppi,
               sid = dt$sid,
               tid = dt$tid,
               m_dem = m_dem,
               m_tot = m_tot, 
               A=A,
               ones.beta = diag(rep(1, ncol(A))))   

# model
model <- "
data{
  n <- length(y)
  da <- dim(A) 
  zeros.beta <- rep(0, da[2])
  df.beta <- da[2]+1
}
## Likelihood  
model{
  for (i in 1:n) {
     y[i] ~ dnorm(yhat[i],tau)
     yhat[i] <- Ab[tid[i],sid[i]] 
  }
 
  for (j in 1:m_dem) {
     Ab[1:da[1],j] <- A %*% b[1:da[2],j]
     b[1:da[2],j] ~ dmnorm(beta.D, Omega.b)
  }
  for (j in (m_dem+1):m_tot) {
     Ab[1:da[1],j] <- A %*% b[1:da[2],j]
     b[1:da[2],j] ~ dmnorm(beta.R, Omega.b)
  }
   
## Priors
  beta.D ~ dmnorm(zeros.beta, Omega.beta)
  beta.R ~ dmnorm(zeros.beta, Omega.beta)
  
  Omega.beta ~ dwish(0.1*ones.beta, df.beta)
  Omega.b ~ dwish(0.1*ones.beta, df.beta)
  tau ~ dgamma(0.001,0.001)
}"


# estimation 
## list parameters
bayes.mod.params <- c("beta.R", "beta.D", "b","tau")

mod.fit <- jags.model(textConnection(model), data=mydata, n.chains = 3, n.adapt=200)
update(mod.fit,10000)
postsam.ppi <- coda.samples(mod.fit, variable.names=bayes.mod.params, thin=10, n.iter=10000)

postsam.ppi.re <- postsam.ppi
for (i in seq_along(postsam.ppi.re)) {
  for (j in seq_along(splnam)) {
    unam <- list(D=paste0("b[", j, ",", 1:m_dem, "]"), R=paste0("b[", j, ",", (m_dem+1):m_tot, "]"))
    for (k in c("D","R")) {
      u <- postsam.ppi.re[[i]][,unam[[k]]]
      meanu <- rowMeans(u)
      postsam.ppi.re[[i]][,paste0("beta.",k,"[",j,"]")] <- meanu
      postsam.ppi.re[[i]][,unam[[k]]] <- postsam.ppi.re[[i]][,unam[[k]]] - meanu 
    }
  }
}


gelman.diag(postsam.ppi.re, multivariate=FALSE)
mcmcplot(postsam.ppi.re)

saveRDS(postsam.ppi.re, file=file.path(psdir, "postsam.ppi.rds"))