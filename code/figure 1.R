## Replication materials for
## Shvetsova, Olga, Andrei Zhirnov, Frank Giannelli, Michael Catalano, and Olivia Catalano. 2021. 
##"Governor's Party, Policies, and COVID-19 Outcomes: Further Evidence of an Effect."  

library(ggplot2)

rm(list=ls())
options(stringsAsFactors = FALSE)

dtdir <- "./data"
psdir <- "./posterior_samples"
outdir <- "./output"

theme_set(theme_bw() + theme(strip.text.x =element_text(size=12, family="sans"), 
                             strip.text.y =element_text(size=12, family="sans"),  
                             axis.title.x =element_text(size=12, family="sans"), 
                             axis.title.y =element_text(size=12, family="sans"), 
                             axis.text=element_text(size=12, family="sans"),
                             legend.text=element_text(size=12, family="sans")                             
))
 

start.date <- as.Date("2020-3-1")
end.date <- as.Date("2020-11-30")

load(file.path(dtdir,"PM_data.RData"))

## sample from the posterior samples
sam <- readRDS(file.path(psdir, "postsam.ppi2casesZ.rds")) 
sam <- do.call("rbind", sam)
part <- sam[sample(1:nrow(sam), size=1000),]

# scale state-level variable as in the estimation sample
xvar <- c("w_rural","w_hispanic","w_black","w_senior","w_obesity","w_smoke","w_asthma","w_heart","w_actph","w_density","w_poverty")
df <- cs.df[order(cs.df$w_rep,cs.df$state),c("state","pop",xvar)]
df$sid <- 1:nrow(df)
df <- data.frame(df[,c("sid","state","pop")], scale(df[,xvar]))

# grid for simulations
starts0 <- data.frame(state="NY", new.end.date = as.Date(c("2020-3-31","2020-7-31")), example=paste0("Period ",1:2))

bulk <- merge(csts.df, starts0, by="state")
bulk <- merge(bulk, df, by="state")
bulk <- within(bulk, tid <- as.numeric(date - new.end.date))

examples <- unique(bulk$example)
names(examples) <- examples

focus <- lapply(examples, function(j) { 
  y <- subset(bulk, example==j)
  y <- within(y, ncases.re <- ncases*1000000/pop) 
  ks <- ksmooth(x=y$tid, y=y$ncases.re,kernel="normal", bandwidth=7)
  a <- data.frame(tid=ks$x, ncases.sm=ks$y)
  z <- merge(y,a, all.x=TRUE)
  subset(z, tid >=-30 & tid <= 90)
})

# initial values
starts <- lapply(focus, function(y) {
  z <- subset(y, tid <=0)
  list(cu = max(z$cases),
       pop = unique(z$pop),
       g = z$ncases[match(0:-20,z$tid)],
       sid = unique(z$sid),
       xval = as.vector(as.matrix(z[1,xvar])) 
       )
})

# loop over estimates and days
levs <- c(lo=0.2, hi=0.7)
results <- lapply(names(levs), function(p) {
  w <- lapply(examples, function(j) {
    x <- starts[[j]]
    y <- apply(part, 1, function(r) {
      z <- numeric()
      h <- exp(-r["gamma"]*(0:20))
      pop <- x$pop
      mu <- r[paste0("mu[",x$sid,"]")]
      cu <- mu*x$cu
      g <- x$g
      xb <- sum(r[paste0("beta[",1:(2+length(x$xval)),"]")]*c(int=1,levs[p],x$xval)) 
      for (t in 1:90) {
        lambda <- log(pop-cu) - log(pop) + xb + log(r["xi"] + sum(g*h))
        prob <- r["r"]/(r["r"] + exp(lambda))
        z[t] <- rnbinom(1, r["r"], prob)
        cu <- cu + mu*z[t]
        g <- c(z[t], g[1:20])
        }
      z*1000000/pop
    })
    v <- apply(y, 1, quantile, probs=c(0.025, 0.25,0.75, 0.975))
    rownames(v) <- c("lb","lbm","ubm","ub")
    sim <- data.frame(sid=x$sid, lev.ppi=format(levs[p], nsmall=2), tid = 1:ncol(v), t(v), pe=rowMeans(y))
    merge(focus[[j]][c("sid","tid","example","state","date","ncases.re","ncases.sm")], sim, by=c("sid","tid"), all=TRUE)
  })
  do.call("rbind", w)
})
results <- do.call("rbind", results)

biweekly <- as.Date("2020-4-5") + seq(0,300,by=10)
results <- within(results, {
  marked <- (2*lbm + ubm)/3
  marked[which(!date %in% biweekly)] <- NA
})


## plotting
pic <- ggplot(results, aes(x=date)) + 
  geom_ribbon(aes(ymin=lbm, ymax=ubm, group=interaction(example,state,lev.ppi)), fill="grey50", alpha=0.7) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=interaction(example,state,lev.ppi)), fill="grey50", alpha=0.3) +  
  geom_line(mapping=aes(y=ncases.re, color="Observed cases")) +
  geom_line(mapping=aes(y=ncases.sm, color="Observed cases"), linewidth=1) +
  geom_point(aes(y=marked, shape=lev.ppi), size=3) +
  scale_x_date(date_labels="%b%e") +
  scale_y_continuous("New cases per million", trans="log", breaks=c(1,10,100,1000,10000), labels=scales::label_comma(accuracy=1)) +
  scale_color_manual(breaks="Observed cases", values="black") +
  labs(x=element_blank(), color=element_blank(), shape="Assumed Protective Policy Index") + 
  scale_shape_manual(breaks=c("0.20","0.70"), values=c(0,15)) +
  facet_wrap(vars(example), ncol=2, scales="free_x") + theme(legend.position="bottom") 
ggsave(file.path(outdir,"figure-1.tiff"), pic, height=6, width=11, compression="lzw")

## for display
pic <- ggplot(results, aes(x=date)) + 
  geom_ribbon(aes(ymin=lbm, ymax=ubm, group=interaction(example,state,lev.ppi)), fill="grey50", alpha=0.7) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=interaction(example,state,lev.ppi)), fill="grey50", alpha=0.3) +  
  geom_line(mapping=aes(y=ncases.re, color="Observed cases")) +
  geom_line(mapping=aes(y=ncases.sm, color="Observed cases"), linewidth=1) +
  geom_point(aes(y=marked, shape=lev.ppi), size=3) +
  scale_x_date(date_labels="%b%e") +
  scale_y_continuous("New cases per million", trans="log", breaks=c(1,10,100,1000,10000), labels=scales::label_comma(accuracy=1)) +
  scale_color_manual(breaks="Observed cases", values="darkblue") +
  labs(x=element_blank(), color=element_blank(), shape="Assumed Protective Policy Index") + 
  scale_shape_manual(breaks=c("0.20","0.70"), values=c(0,15)) +
  facet_wrap(vars(example), ncol=2, scales="free_x") + theme(legend.position="bottom") 
ggsave("images/figure-1.png", pic, height=6, width=11)