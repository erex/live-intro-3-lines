## ----setup, include=FALSE-------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
options(width=100)


## ---- echo=TRUE, eval=TRUE, message=FALSE---------------------------------------
library(Distance)
data("LTExercise")
# Check that it has been imported correctly
head(LTExercise, n=3)


## -------------------------------------------------------------------------------
LTExercise[100:102, ]


## ---- message=FALSE-------------------------------------------------------------
conversion.factor <- convert_units("meter", "kilometer", "square kilometer")
# Truncate at 20metres
lt.hn.t20m <- ds(data=LTExercise, key="hn", adjustment=NULL, truncation=20, 
                convert_units=conversion.factor)
summary(lt.hn.t20m)


## ---- message=FALSE-------------------------------------------------------------
# Truncate 10% of largest distances
lt.hn.t10per <- ds(data=LTExercise, key="hn", adjustment=NULL, truncation="10%", 
                convert_units=conversion.factor)
summary(lt.hn.t10per)


## ---- echo=TRUE, eval=TRUE------------------------------------------------------
# Divide plot window
par(mfrow=c(1,2))
plot(lt.hn.t20m, main="Truncation 20m")
plot(lt.hn.t10per, main="Truncation 10%")


## ---- echo=TRUE, eval=TRUE, message=FALSE---------------------------------------
# Fit a few different models
# Half normal model, no adjustments, no truncation
lt.hn <- ds(data=LTExercise, key="hn", adjustment=NULL, convert_units=conversion.factor)
# Half normal model, cosine adjustments, truncation at 20m
lt.hn.cos.t20m <- ds(data=LTExercise, key="hn", adjustment="cos", truncation=20, 
                     convert_units=conversion.factor)
# Uniform model, cosine adjustments, truncation at 20m
lt.uf.cos.t20m <- ds(data=LTExercise, key="unif", adjustment="cos", 
                     truncation=20, convert_units=conversion.factor)
# Hazard rate model, no adjustments, truncation at 20m
lt.hr.t20m <- ds(data=LTExercise, key="hr", adjustment="poly", truncation=20,
                 convert_units=conversion.factor)


## ---- echo=FALSE, eval=TRUE-----------------------------------------------------
# This block of code is quite complex, but not because of performing
#   distance sampling analysis.  Instead it is used to make the tables
#   for the solution look attractive.
lt.tab <- data.frame(DetectionFunction=c("Half-normal",
                                         "Half-normal","Uniform","Hazard rate"),
                     Adjustments=c("None","Cosine","Cosine","Polynomial"), 
                     Terms=c(0,0,1,0), Truncation=c(35.8,20,20,20), AIC=rep(NA,4), Pa=rep(NA,4), Density=rep(NA,4), D.CV=rep(NA,4), Lower.CI=rep(NA,4), Upper.CI=rep(NA,4))

get.results.f <- function(fit.model) {   
  return(c(AIC=summary(fit.model$ddf)$aic,
         Pa=fit.model$dht$individuals$average.p,
         D=fit.model$dht$individuals$D$Estimate,
         D.CV=fit.model$dht$individuals$D$cv,
         lCL=fit.model$dht$individuals$D$lcl,
         uCL=fit.model$dht$individuals$D$ucl))
}

lt.tab[1,5:10] <- get.results.f(lt.hn)
lt.tab[2,5:10] <- get.results.f(lt.hn.cos.t20m)
lt.tab[3,5:10] <- get.results.f(lt.uf.cos.t20m)
lt.tab[4,5:10] <- get.results.f(lt.hr.t20m)


## ---- echo=FALSE, eval=TRUE-----------------------------------------------------
# Print results
knitr::kable(lt.tab, digits=3,
             caption="Results for simulated data with differing truncation and detection functions.")


## ---- echo=TRUE, eval=TRUE------------------------------------------------------
# Divide plot window
par(mfrow=c(2,2))
# Plot detection functions
plot(lt.hn, main="HN, no truncation")
plot(lt.hn.cos.t20m, main="HN, truncation at 20m")
plot(lt.uf.cos.t20m, main="Uniform, truncation at 20m")
plot(lt.hr.t20m, main="HR, truncation at 20m")


## ---- fig.height=4, fig.width=4, message=FALSE----------------------------------
data(capercaillie)
head(capercaillie, n=3)
conversion.factor <- convert_units("meter", "kilometer", "hectare")
# Fit a half normal model with no adjustments and no truncation
caper.hn <- ds(data=capercaillie, key="hn", adjustment=NULL, 
               convert_units=conversion.factor)
# Plot with lots of bins, each of width 2m
plot(caper.hn, nc=40)


## ---- message=FALSE-------------------------------------------------------------
# Half normal model 
caper.hn.cos <- ds(data=capercaillie, key="hn", adjustment="cos",
                   convert_units=conversion.factor)
# Hazard rate model  
caper.hr.cos <- ds(data=capercaillie, key="hr", adjustment="cos",
                   convert_units=conversion.factor)
# Uniform model  
caper.uf.cos <- ds(data=capercaillie, key="unif", adjustment="cos",
                   convert_units=conversion.factor)


## ---- echo=TRUE, eval=TRUE, results="hide"--------------------------------------
# Divide plot window
par(mfrow=c(3,2))
par(mar=c(4,4,.2,.1))
plot(caper.hn.cos, main="Half normal")
gof_ds(caper.hn.cos)
plot(caper.hr.cos, main="Hazard rate")
gof_ds(caper.hr.cos)
plot(caper.uf.cos, main="Uniform")
gof_ds(caper.uf.cos)


## -------------------------------------------------------------------------------
knitr::kable(summarize_ds_models(caper.hn.cos, caper.hr.cos, caper.uf.cos, output="plain"),
               caption="Summary of results of Capercaillie analysis.", digits = 3)


## ---- echo=FALSE, eval=TRUE-----------------------------------------------------
# Harvest results
caper.tab <- data.frame(DetectionFunction=c("Half-normal","Hazard rate",
                                            "Uniform"), 
                        AIC=rep(NA,3), Pa=rep(NA,3), Density=rep(NA,3), 
                        D.CV=rep(NA,3), Lower.CI=rep(NA,3), Upper.CI=rep(NA,3))
caper.tab[1,2:7] <- get.results.f(caper.hn.cos)
caper.tab[2,2:7] <- get.results.f(caper.hr.cos)
caper.tab[3,2:7] <- get.results.f(caper.uf.cos)

# Print results
knitr::kable(caper.tab, digits=3,
             caption="Capercaillie point estimates of density and associated measures of precision.")


## ---- echo=TRUE, eval=TRUE, fig.width=5, fig.height=4, message=FALSE------------
# Specify (uneven) cutpoint for bins
bins <- c(0, seq(from=7.5, to=67.5, by=10), 80)
# Check bins
bins
# Specify model with binned distances
caper.hn.bin <- ds(data=capercaillie, key="hn", adjustment="cos", cutpoints=bins,
                   convert_units=conversion.factor)
# Plot
plot(caper.hn.bin, main="Capercaillie, binned distances")
# Summarise results
caper.hn.bin$dht$individuals$summary
caper.hn.bin$dht$individuals$D[1:6]

