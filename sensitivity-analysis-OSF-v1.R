rm(list=ls())
library("lavaan")

## download data from OSF into the same working directory as this R script,
## then unzip into "osfstorage-archive" folder
# load data
df2 <- read.csv("osfstorage-archive/Hoard_Study2Data.csv")
colnames(df2)
nrow(df2)

#standardize variables for simple slopes and mediation
df2 <- within(df2, {
  policyz <- as.numeric(scale(policy_c))
  condz <- as.numeric(scale(condition))
  ssmz <- as.numeric(scale(ssm_c))
  overthoardz <- as.numeric(scale(overthoard_c))
})

# Check SES is standardized with mean 0 and sd 1
mean(df2$sescomp,na.rm=TRUE)
sd(df2$sescomp,na.rm=TRUE)

# with sescomp as covariate
model.m <- '
 ssmz~a*condz+sescomp
 overthoardz~b*ssmz+c*condz+sescomp
 ie := a*b
 de := c
'

# restrict to high SES: at least 1 sd above mean
df2.h <- subset(df2,sescomp>=1)
fit.h <- lavaan::sem(model=model.m, data=df2.h)
subset(parameterEstimates(fit.h),subset=label %in% c("ie","de"))
rm(fit.h)

# restrict to low SES: at least 1 sd below mean
df2.l <- subset(df2,sescomp<=-1)
fit.l <- lavaan::sem(model=model.m, data=df2.l)
subset(parameterEstimates(fit.l),subset=label %in% c("ie","de"))
rm(df2.l,fit.l)

# sensitivity analysis to unmeasured M-Y confounding
model.sens <- '
 ssmz~a*condz+sescomp
 overthoardz~b*ssmz+c*condz+sescomp
 ie := a*b
 de := c
 ## correlated residuals
 overthoardz ~~ cov_ym*ssmz
 ## residual variances
 overthoardz ~~ var_y*overthoardz
 ssmz ~~ var_m*ssmz
 var_y > 0
 var_m > 0
 ## fixed correlations
 cov_ym/sqrt(var_y*var_m) == rho
'
rho.vals <- seq(from=-0.7,to=0.7,by=0.05)
RES.sens <- lapply(rho.vals, function(rho) {
  model.rho <- model.sens
  model.rho <- gsub(pattern="rho",replacement=rho,x=model.rho)
  #high SES
  fit.h <- lavaan::sem(model=model.rho, data=df2.h)
  res.rho <- cbind(
    "rho"=rho,
    subset(parameterEstimates(fit.h),subset=label %in% c("ie","de")))
  return(res.rho)
})
RES.sens <- do.call(rbind,RES.sens)

# make plot
RES.ie <- subset(RES.sens,label=="ie")
RES.de <- subset(RES.sens,label=="de")
pdf(paste0("plot-sensitivity-ses1.pdf"),width=6,height=4)
plot(RES.de$rho-0.005,RES.de$est,
     pch=18, cex=0.9, 
     ylim=range(c(RES.de$ci.lower,RES.de$ci.upper,
                  RES.ie$ci.lower,RES.ie$ci.upper)),
     xlim=range(rho.vals),
     xlab=expression(rho), ylab="Effect Estimate", 
     main="Sensitivity Analysis")
points(RES.ie$rho+0.005,RES.ie$est,
       pch=20, cex=0.9, 
       col=2)
abline(h=0,v=0, lty=3, lwd=0.5)
# 95% CIs
for (j in 1:nrow(RES.de)) {
  exclude0.de <- sign(RES.de[j,"ci.lower"])==sign(RES.de[j,"ci.upper"])
  lines(rep(RES.de[j,"rho"]-0.005,2),
        RES.de[j,c("ci.lower","ci.upper")],
        col=1,
        lwd=ifelse(exclude0.de,0.5,1),
        lty=ifelse(exclude0.de,2,1))
  exclude0.ie <- sign(RES.ie[j,"ci.lower"])==sign(RES.ie[j,"ci.upper"])
  lines(rep(RES.ie[j,"rho"]+0.005,2),
        RES.ie[j,c("ci.lower","ci.upper")],
        col=2,
        lwd=ifelse(exclude0.ie,0.5,1),
        lty=ifelse(exclude0.ie,2,1))
}
legend("bottom",
       legend= c("Indirect","Direct"),
       col=c(2:1),pch=c(20,18),bty="n",cex=0.9)
dev.off()
