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
df$sid <- 1:nrow(df)
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

xvar <- c("w_rural","w_hispanic","w_black","w_senior","w_obesity","w_smoke","w_asthma","w_heart","w_actph","w_density","w_poverty")
uvar <- c("negtests.pc","w_popdoc")
mat <- as.matrix(df[,c(xvar,uvar)])
mat <- scale(mat)

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
dt <- data.frame(df[c("state","sid","pop")], mat[,c(uvar,xvar)])
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

A <- cbind(int=1, ppi=dt$ppi.l3)
G <- as.matrix(dt[,lg]) 

X.prep <- unique(dt[c("sid","pop",xvar,uvar)])
maxc <- aggregate(cases.l1 ~ sid, dt, max)
X.prep <- within(X.prep, {
  cap <- pop/maxc$cases.l1[match(sid,maxc$sid)]
})
pop <- X.prep[order(X.prep$sid),"pop"]
U <- cbind(int=1,X.prep[order(X.prep$sid),uvar])
cap <- X.prep[["cap"]]
X <- X.prep[order(X.prep$sid),xvar]

mydata <- list(y = dt$ncases,
               cu= dt$cases/100000,
               sid = dt$sid,
               pop= pop/100000,
               A=A,
               G=G,
               U=U,
               X=X,
               cap=cap,
               ones.coef = diag(rep(1, ncol(A) + ncol(X) + ncol(U))),
               n = nrow(dt),
               s = nrow(U),
               al =ncol(A),
               xl = ncol(X),
               ul = ncol(U),
               gl = ncol(G),
               zeros.coef = rep(0, ncol(A) + ncol(X) + ncol(U)),
               df.coef = ncol(A) + ncol(X) + ncol(U) + 1
               ) 
 
# model
cat("model
{
for (i in 1:n) {
     y[i] ~ dnegbin(p[i],r)
     p[i] <- r/(r+lambda[i]) 
     log(lambda[i]) <- log(pop[sid[i]] - mu[sid[i]]*cu[i]) - log(pop[sid[i]]) + Ab[i] + Xb[sid[i]] + log(xi + Ghc[i])
  }
  
  Ghc <- G %*% hc
  
## Priors
  for (j in 1:gl) {
    hc[j] <- exp(-gamma*(j-1)) 
    }
  
  for (j in 1:s) {
    mu[j] ~ dlnorm(exp(Ud[j]), tau)T(,cap[j])
    }
  
  Ud <- U %*% delta
  Ab <- A %*% beta[1:al]
  Xb <- X %*% beta[(al+1):(al+xl)]

  beta <- coef[1:(al+xl)]
  delta <- coef[(al+xl+1):(al+xl+ul)]
  coef ~ dmnorm(zeros.coef, Omega.coef)
  Omega.coef ~ dwish(0.1*ones.coef, df.coef)
  xi ~ dlnorm(5, 0.001)
  gamma ~ dgamma(0.01,0.01)
  tau ~ dgamma(0.01,0.01)
  r ~ dgamma(0.01,0.01)
}", file="jags-model.jag")

# estimation 
## list parameters
bayes.mod.params <- c("beta","gamma","delta","xi","r","mu")

### run mcmc
mod.fit <- jags.parallel(data=mydata, parameters.to.save = bayes.mod.params, n.iter=100000, n.chains=3, model.file="jags-model.jag", n.thin=20)

samples<- as.mcmc(mod.fit)

saveRDS(samples, file=file.path(psdir,"postsam.ppi2casesZ.rds"))



