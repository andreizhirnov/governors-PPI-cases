## Replication materials for
## Shvetsova, Olga, Andrei Zhirnov, Frank Giannelli, Michael Catalano, and Olivia Catalano. 2021. 
##"Governor's Party, Policies, and COVID-19 Outcomes: Further Evidence of an Effect." 

library(ggplot2)
library(splines)

rm(list=ls())
options(stringsAsFactors = FALSE)

theme_set(theme_bw() + theme(strip.text.x =element_text(size=12, family="sans"), 
                             strip.text.y =element_text(size=12, family="sans"),  
                             axis.title.x =element_text(size=12, family="sans"), 
                             axis.title.y =element_text(size=12, family="sans"),
                             plot.caption =element_text(size=12, family="sans"), 
                             axis.text=element_text(size=12, family="sans"),
                             legend.text=element_text(size=12, family="sans"),
                             legend.title=element_text(size=12, family="sans")
                             
))

dtdir <- "./data"
psdir <- "./posterior_samples"
outdir <- "./output"
  

start.date <- as.Date("2020-3-1")
end.date <- as.Date("2020-11-30")

load(file.path(dtdir,"PM_data.RData"))
sam <- readRDS(file=file.path(psdir, "postsam.ppi.rds"))

# states vector
states <- sort(cs.df$state)

sample <- csts.df
sample$tid <- as.numeric(sample$date - start.date + 1)
 
## for smoothing
bulk <- lapply(states, function(j) {
  y <- subset(sample, state==j, select=c("state","tid","ncases"))
  ks <- ksmooth(x=y$tid, y=y$ncases,kernel="normal", bandwidth=7)
  data.frame(state=j, tid=ks$x, ncases.sm=ks$y)
})
bulk <- do.call("rbind", bulk)
sample <- merge(sample, bulk, by=c("state","tid"))

## PPi at peaks
peaks <- aggregate(ncases.sm ~ state, data=subset(sample, date >= start.date & date <= end.date), FUN=max)
peaks <- merge(sample,peaks)
peaks <- merge(peaks, cs.df, by="state")
peaks <- within(peaks, {
  ncases.sm.re <- ncases.sm*10000/pop
  gov <- ifelse(w_rep==1, "Republican","Democrat")
  })
peaks <- subset(peaks, !is.na(ppi), select=c("state","gov","ppi","ncases.sm.re"))

peaks.comb <- reshape(peaks, direction="long", v.names="value", timevar="type", times=c("ppi","ncases.sm.re"), varying=match(c("ppi","ncases.sm.re"),colnames(peaks)))

## changes start here
(means <- aggregate(list(mean=peaks.comb$value), by=peaks.comb[c("gov","type")], FUN=mean))

pic <- ggplot(peaks.comb, aes(x=value, y=gov))+ 
  geom_boxplot(width=0.25)+
  geom_point(aes(x=mean, y=gov), data=means, shape=15, size=5, color="black") +
  labs(x=element_blank(), y="Governor", caption="Black squares indicate average values.") +
  facet_wrap(vars(type), scales="free_x", ncol=2, labeller=labeller(type=c("ppi"="Protective Policy Index at the peak of the pandemic", 
                                                          "ncases.sm.re"="The highest 7-day average number of cases per 10,000"))) +
  theme(legend.position="none", axis.ticks.y=element_blank(), axis.text.y=element_text(angle=90, hjust=0.5) )
ggsave(file.path(outdir, "figure-3.tiff"), pic, height=4, width=10, compression="lzw")
##

## work with the posterior samples
tdf <- subset(csts.df, date >= start.date & date <= end.date)
tdf$tid <- as.numeric(tdf$date - start.date) + 1
uda <- sort(unique(tdf$tid))
knots <- uda[which(uda %% 14 == 0)]
spl <- ns(uda, knots=knots, intercept=TRUE)

## background -- state-specific trajectories
bg <- merge(tdf, cs.df, by="state")
bg <- within(bg, gov <- ifelse(w_rep==1,"rep","dem"))
bg <- bg[c("state","gov","date","ppi")]
bg$type <- "Average Protective Policy Index"

## the part with predictions
post <- lapply(sam, function(z) {
  w <- lapply(list(dem="D",rep="R"), function(j) {
       z[,paste0("beta.",j,"[",1:ncol(spl),"]")] %*% t(spl)
      })
  w[["diff"]] <- w[["dem"]]-w[["rep"]]
  w
})
pick <- lapply(c(dem="dem", rep="rep", diff="diff"), function(j) {
  w <- do.call("rbind", c(lapply(post, "[[", j), list(make.row.names=FALSE)))
  s <- apply(w, 2, quantile, probs=c(0.025,0.25,0.75,0.975))
  m <- data.frame(date=start.date-1+1:ncol(s), t(s))
  colnames(m) <- c("date", "lb","lbm", "ubm", "ub")
  m
}) 

pick$dem$gov <- "dem"
pick$rep$gov <- "rep"
pick$dem$type <- pick$rep$type <- "Average Protective Policy Index"
pick$diff$gov <- "diff"
pick$diff$type <- "Difference in average Protective Policy Index"
pick <- do.call("rbind", pick) 

pick <- within(pick, {
  marked <- (ubm+lbm)/2
  marked[which(gov=="diff")] <- NA
  marked[which(!date %in% as.Date(paste0("2020-",4:11,"-1")))] <- NA
  marker <- gov
  marker[which(is.na(marked))] <- NA
  })

## draw
reflines <- data.frame(type="Difference in average Protective Policy Index", yint=0)

pic <- ggplot(pick, aes(x=date, group=gov)) +
  geom_line(aes(group=interaction(gov,state), y=ppi, linetype=gov), data=bg, alpha=0.5, color="gray") +
  geom_point(aes(y=marked, shape=marker), size=4) +
  geom_ribbon(aes(ymin=lb, ymax=ub), color=NA, fill="grey40") +
  scale_linetype_manual("Protective Policy Index by state",breaks=c("dem","rep"), values=c("solid","longdash"),labels=c("Governor is a Democrat","Governor is a Republican")) + 
  scale_shape_manual("Average Protective Policy Index",breaks=c("dem","rep"), values=c(6,16),labels=c("Governor is a Democrat","Governor is a Republican")) + 
  scale_x_date(date_labels="%b %d") +
  labs(x=element_blank(), y=element_blank()) +
  geom_hline(aes(yintercept=yint), data=reflines, color="black") + 
  facet_wrap(vars(type), scales="free_y", ncol=1)
ggsave(file.path(outdir,"figure-2.tiff"), pic, height=6, width=10, compression="lzw")
