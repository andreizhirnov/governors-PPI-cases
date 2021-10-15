## Replication materials for
## Shvetsova, Olga, Andrei Zhirnov, Frank Giannelli, Michael Catalano, and Olivia Catalano. 2021. 
##"Governor's Party, Policies, and COVID-19 Outcomes: Further Evidence of an Effect."  

library(coda)

rm(list=ls())
options(stringsAsFactors = FALSE)

dtdir <- "./data"
psdir <- "./posterior_samples"
outdir <- "./output"
  
load(file.path(dtdir,"PM_data.RData"))

# states vector
state <- cs.df[order(cs.df$w_rep,cs.df$state),"state"]

template <- read.csv("./code/template_for_tables.csv")

guide <- c("V","W","X","Y","Z","Z_2","Z_3")
names(guide) <- guide

est.list <- lapply(guide, function(j) {
  readRDS(file=paste0(psdir, "/postsam.ppi2cases", j, ".rds")) 
  }) 
est.list$V <- est.list$V[c(1,3)] 

  
mat.list <- lapply(est.list, function(x) do.call("rbind",x))

DIC <- lapply(mat.list, function(x) {
  dev <- x[,"deviance"] 
  mean(dev) + var(dev)/2
})
names(DIC) <- paste0("val.", names(DIC))
dic.df <- data.frame(
  param="DIC",type=3, DIC
)

param.list <- lapply(guide, function(x) {
  cm <- colMeans(mat.list[[x]])
  qs <- apply(mat.list[[x]], 2, quantile, probs=c(0.025, 0.975))
  means <- data.frame(param = names(cm), type = 1,
                      val=format(round(cm, 3), trim=TRUE, 4,3))
  qsd <- data.frame(param = colnames(qs), type = 2,
                      val=paste0("[",format(round(qs[1,], 3), trim=TRUE, 4,3),";",
                                 format(round(qs[2,], 3), trim=TRUE, 4,3),"]"))
  combo <- rbind(means,qsd)
  colnames(combo) <- c("param","type",paste0("val.",x))
  subset(combo, param != "deviance")
})

param <- Reduce(function(x,y) merge(x,y,by=c("param","type"), all=TRUE), param.list)
export <- rbind(param, dic.df)
export <- merge(template,export)
export <- export[order(export$sort),]

write.csv(export, file.path(outdir,"est_cases.csv"), na="")

## underreporting 
mus <- data.frame(state=state, param=paste0("mu[",seq_along(state),"]"))
mus <- merge(mus, param.list[["Z"]])
export <- reshape(mus, direction="wide", idvar="state", v.names="val.Z", timevar="type")
export <- merge(cs.df[c("state","negtests.pc","w_popdoc")],export)
write.csv(export, file.path(outdir,"est_mus.csv"), na="")

## governor to PPI
sam <- readRDS(file=file.path(psdir, "postsam.ppi.rds"))
temp <- do.call("rbind", sam)
stat.mean <- colMeans(temp)
stat.q <- apply(temp, 2, quantile, probs=c(0.025, 0.975))
stat.rhat <- gelman.diag(sam, multivariate=FALSE)[["psrf"]][,1]
stat.tab <- data.frame(coef=names(stat.mean), 
                       mean=stat.mean, lb=stat.q[1,], 
                       ub=stat.q[2,], rhat=stat.rhat)
stat.nam <- list(beta.D=grep("^beta.D", stat.tab$coef, value=TRUE),
                 beta.R=grep("^beta.R", stat.tab$coef, value=TRUE),
                 b=grep("^b\\[", stat.tab$coef, value=TRUE))

stat.nam[["rest"]] <- setdiff(stat.tab$coef, c(stat.nam$beta.D,
                                               stat.nam$beta.R,
                                               stat.nam$b))

stat.nam <- lapply(names(stat.nam), function(x) {
  val <- stat.nam[[x]]
  if (x=="b") {
    parts <- strsplit(val, ",")
    beg <- sapply(parts, "[", 1)
    end <- sapply(parts, "[", 2)
    sta <- as.numeric(as.character(substring(end, 1, nchar(end)-1)))
    prep <- paste0(beg, ",", state[sta], "]")
    data.frame(coef=val, param=prep)
  } else {
    data.frame(coef=val, param=val)
  }
})

stat.nam <- do.call("rbind", stat.nam)
stat.nam$sort <- 1:nrow(stat.nam)
stat.tab <- merge(stat.nam, stat.tab, by="coef")
stat.tab <- stat.tab[order(stat.tab$sort),]
write.csv(stat.tab, file.path(outdir,"est_ppi.csv"), na="")
