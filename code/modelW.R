## Replication materials for
## Shvetsova, Olga, Andrei Zhirnov, Frank Giannelli, Michael Catalano, and Olivia Catalano. 2021. 
##"Governor's Party, Policies, and COVID-19 Outcomes: Further Evidence of an Effect."  


library(R2jags)


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

# party
repp <- df$w_rep

# other variables
reg <- df$w_region
regl <- lapply(unique(reg), function(x) {
  as.numeric(reg==x) 
})
regl <- do.call("cbind", regl)
colnames(regl) <- paste0("region.", unique(reg))

othvar <- c("w_rural","w_hispanic","w_black","w_senior","w_obesity","w_smoke","w_asthma","w_heart","w_actph","w_density","w_poverty")
mat <- as.matrix(df[,othvar])
# mat <- cbind(mat, regl)
mat <- scale(mat)
xvar <- colnames(mat)

# time varying variables
tdf <- subset(csts.df, date >= start.date & date <= end.date)
tdf$tid <- as.numeric(tdf$date - start.date) + 1

tdf.l3 <- within(csts.df, {
  tid <- as.numeric(date - start.date) + 4
  ppi.l3 <- ppi
})[c("state","tid","ppi.l3")]

tdf.l1 <- within(csts.df, {
  tid <- as.numeric(date - start.date) + 2
  cases.l1 <- cases
  ncases.l1 <- ncases
})[c("state","tid","cases.l1","ncases.l1")]

# expand the state-level dataset
dt <- data.frame(state,repp,mat,pop=df$pop)
dt <- merge(tdf, dt, by="state")
dt <- merge(dt, tdf.l3, by=c("state","tid"), all.x=TRUE)
dt <- merge(dt, tdf.l1, by=c("state","tid"), all.x=TRUE)
lg <- paste0("G.",1:21)
G.li <- lapply(seq_along(lg), function(x) {
  y <-  data.frame(
    state = csts.df$state,
    tid = as.numeric(csts.df$date - start.date) + 1 + x,
    val = csts.df$ncases 
  ) 
  colnames(y)[3] <- lg[x]
  y
})
G.prep <- Reduce(function(x,y) merge(x,y), G.li)
dt <- merge(dt, G.prep, by=c("state","tid"), all.x=TRUE)

dt <- within(dt, {
  lns <- log(pop-cases.l1) -log(pop)
})
 
dt$sid <- match(dt$state,state)
A <- cbind(int=1, ppi=dt$ppi.l3)
G <- as.matrix(dt[,lg]) 

X.prep <- unique(dt[c("sid","pop",xvar)])
pop <- X.prep[order(X.prep$sid),"pop"]

lz <- paste0("Z.",state)
Z.prep <- within(dt, {
  Z <- ncases.l1
})

  
dt$mid <- as.numeric(format(dt$date, "%m"))
months <- sort(unique(dt$mid))
dt$mid <- match(dt$mid,months)

cap.prep <- subset(dt, cases.l1 >0)
cap <- min(cap.prep$pop/cap.prep$cases.l1)


ones.gamma <- matrix(0,ncol(G),ncol(G))
ind <- which(lower.tri(ones.gamma, diag=TRUE), arr.ind=TRUE)
ones.gamma[ind] <- 1

mydata <- list(y = dt$ncases,
               cu= dt$cases/100000,
               sid = dt$sid,
               pop= pop/100000,
               A=A,
               G=G,
               cap=cap,
               ones.beta = diag(rep(1, ncol(A)))) 
 
# model
cat("
data{
  n <- length(y) 
  da <- dim(A)
  dg <- dim(G)
  zeros.beta <- rep(0, da[2])
  df.beta <- da[2]+1
  zeros.gamma <- rep(0, dg[2])
  df.gamma <- dg[2]+1
  for (i in 1:length(y)) {
  lnoff[i] <- log(pop[sid[i]] - cu[i]) - log(pop[sid[i]]) 
  }
}
## Likelihood  
model{
  for (i in 1:n) {
     y[i] ~ dnegbin(p[i],r)
     p[i] <- r/(r+lambda[i]) 
     log(lambda[i]) <-lnoff[i]  + A[i,] %*% beta + log(xi + Ghc[i])
  }
  
  Ghc <- G %*% hc
  
## Priors
  for (j in 1:dg[2]) {
   hc[j] <- exp(-gamma*(j-1)) 
  }
  
  beta ~ dmnorm(zeros.beta, Omega.beta)
  Omega.beta ~ dwish(0.1*ones.beta, df.beta)
  xi ~ dlnorm(5, 0.001)
  gamma ~ dgamma(0.01,0.01)
  r ~ dgamma(0.01,0.01)
}", file="jags-model.jag")

# estimation 
## list parameters
bayes.mod.params <- c("beta","gamma","xi","r")

### run mcmc 
mod.fit <- jags.parallel(data=mydata, parameters.to.save = bayes.mod.params, n.iter=30000, n.chains=3, model.file="jags-model.jag", n.thin=6)
samples<- as.mcmc(mod.fit)

saveRDS(samples, file=file.path(psdir,"postsam.ppi2casesW.rds"))



