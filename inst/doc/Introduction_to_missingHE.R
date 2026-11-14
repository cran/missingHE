## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
has_jags <- isTRUE(tryCatch({ loadNamespace("rjags"); TRUE }, error = function(e) FALSE))

knitr::opts_chunk$set(
  eval = has_jags,
  prompt = TRUE,
  highlight = FALSE,
  background = '#FFFFFF',
  collapse = TRUE,
  comment = "#>"
)

library(missingHE)
library(bookdown)
library(ggplot2)

options(width = 300)
set.seed(1014)

## ----menss, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MenSS$u.0 <- round(MenSS$u.0, digits = 3)
n <- nrow(MenSS)
n1 <- nrow(MenSS[MenSS$trt=="1",])
n2 <- nrow(MenSS[MenSS$trt=="2",])
MenSS_outcomes <- MenSS[, c("sex_inst","sti","e","c","trt")]
MenSS$trt <- factor(MenSS$trt)
levels(MenSS$trt) <- c("SoC", "MenSS")
MenSS_outcomes$trt <- factor(MenSS_outcomes$trt)
levels(MenSS_outcomes$trt) <- c("SoC", "MenSS")
MenSS_ec <- MenSS[,c("e","c","trt")]
MenSS_ec_cc <- MenSS_ec[complete.cases(MenSS_ec),]
n_ec_cc <- nrow(MenSS_ec_cc)
n1_ec_cc <- nrow(MenSS_ec_cc[MenSS_ec_cc$trt=="1",])
n2_ec_cc <- nrow(MenSS_ec_cc[MenSS_ec_cc$trt=="2",])
n_e_ones <- nrow(MenSS_ec_cc[MenSS_ec_cc$e==1,])
n1_e_ones <- nrow(MenSS_ec_cc[MenSS_ec_cc$e==1&MenSS_ec_cc$trt=="1",])
n2_e_ones <- nrow(MenSS_ec_cc[MenSS_ec_cc$e==1&MenSS_ec_cc$trt=="2",])
n_c_zeros <- nrow(MenSS_ec_cc[MenSS_ec_cc$c==0,])
n1_c_zeros <- nrow(MenSS_ec_cc[MenSS_ec_cc$c==0&MenSS_ec_cc$trt=="1",])
n2_c_zeros <- nrow(MenSS_ec_cc[MenSS_ec_cc$c==0&MenSS_ec_cc$trt=="2",])

hist_e <- ggplot(MenSS_outcomes, aes(x=e, fill=trt)) + geom_histogram(color="black", binwidth = 0.05) + facet_wrap(~trt, labeller=label_parsed) + scale_fill_manual(values=c("grey","grey")) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme(panel.grid.major = element_blank(), legend.key = element_rect(fill = "white"), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position="none",
axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))+ ylab("Frequency") + xlab("QALYs") + scale_x_continuous(breaks=seq(0.6,1,0.2))
hist_c <- ggplot(MenSS_outcomes, aes(x=c, fill=trt)) + geom_histogram(color="black", binwidth = 100) + facet_wrap(~trt) + scale_fill_manual(values=c("grey","grey")) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme(panel.grid.major = element_blank(), legend.key = element_rect(fill = "white"), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position="none",
axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))+ ylab("Frequency") + xlab("Tcosts") + scale_x_continuous(breaks=seq(0,1000,1000))

## ----tbl1-hide, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# rbind(head(MenSS), tail(MenSS))

## ----tbl1, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(rbind(head(MenSS), tail(MenSS)), row.names = FALSE)

## ----fig1e, echo=FALSE, eval=TRUE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Histograms of the observed data distributions for the QALYs in the SoC and MenSS arm.", out.width='100%'--------------------------------------------------------------------------------------------------------------
hist_e

## ----fig1c, echo=FALSE, eval=TRUE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Histograms of the observed data distributions for the Total costs in the SoC and MenSS arm.", out.width='100%'--------------------------------------------------------------------------------------------------------
hist_c

## ----sel1, echo=FALSE, message=FALSE, warning=FALSE, results='hide'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sm1_mar <- selection(data = MenSS, dist_e = "norm", dist_c = "norm",
         model.eff = e ~ trt + u.0, model.cost = c ~ trt + e,
         model.me = me ~ u.0, model.mc = mc ~ 1, 
         type = "MAR", n.iter = 1000, ref = 2)

## ----pat1, echo=FALSE, message=FALSE, warning=FALSE, results='hide'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pm1_mar <- pattern(data = MenSS, dist_e = "logis", dist_c = "norm", 
         model.eff = e ~ trt + u.0, model.cost = c ~ trt + e, 
         type = "MAR", restriction = "CC", n.iter = 1000, ref = 2)

## ----hur1, echo=FALSE, message=FALSE, warning=FALSE, results='hide'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
hm1_mar <- hurdle(data = MenSS, dist_e = "beta", dist_c = "gamma", 
         model.eff = e ~ trt + u.0, model.cost = c ~ trt + e, 
         model.se = se ~ 1, model.sc = sc ~ 1, 
         se = 1, sc = 0, type = "SCAR", 
         n.iter = 1000, ref = 2)

## ----diag1, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# diagnostic(x = sm1_mar, type = "traceplot", param = "sd.e")

## ----fig2, echo=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Checking convergence using the diagnostic function for a model fitted in missingHE, for example through inspection of the traceplots.", out.width='100%'--------------------------------------------------------------------------
diagnostic(x = sm1_mar, type = "traceplot", param = "sd.e")

## ----ppc1, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ppc(x = sm1_mar, type = "dens_overlay", outcome = "effects", ndisplay = 20)

## ----fig3, echo=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Posterior predictive graphical checks using the ppc function for a model fitted in missingHE, for example through inspection of overlayed densities.", out.width='100%'-----------------------------------------------------------
ppc(x = sm1_mar, type = "dens_overlay", outcome = "effects", ndisplay = 20)

## ----pic1, echo=TRUE, message=FALSE, warning=FALSE, results='hide'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pic_sm1_mar <- pic(x = sm1_mar, criterion = "waic", cases = "cc")

## ----pic1show, echo=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pic_sm1_mar$waic

## ----names, echo=TRUE, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# names(sm1_mar)

## ----namesshow, echo=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
names(sm1_mar)

## ----namesmodel, echo=TRUE, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# names(sm1_mar$model_output$`model`)

## ----namesmodelshow, echo=FALSE-------------------------------------------------------------------
options(width = 100)
names(sm1_mar$model_output$`model`)

## ----postsum, echo=TRUE, eval=FALSE---------------------------------------------------------------
# postsum(sm1_mar)

## ----postsumshow, echo=FALSE----------------------------------------------------------------------
postsum(sm1_mar)

## ----print, echo=TRUE, eval=FALSE-----------------------------------------------------------------
# print(sm1_mar)

## ----printshow, echo=TRUE, eval=FALSE-------------------------------------------------------------
# print(sm1_mar)

## ----coef, echo=TRUE, eval=FALSE------------------------------------------------------------------
# coef(sm1_mar)

## ----coefshow, echo=FALSE-------------------------------------------------------------------------
coef(sm1_mar)

## ----plot1, echo=TRUE, eval=FALSE, warning=FALSE--------------------------------------------------
# plot(sm1_mar, class = "scatter", outcome = "effects", trt = "MenSS")

## ----figplot, echo=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Scatter plots of observed and imputed effectiveness data in the MenSS arm under normal distributions using a selection model.", out.width='100%'----
plot(sm1_mar, class = "scatter", outcome = "effects", trt = "MenSS")

## ----summary, echo=TRUE, eval=FALSE---------------------------------------------------------------
# summary(sm1_mar)

## ----summaryshow, echo=FALSE----------------------------------------------------------------------
diff_e_m1 <- round(mean(unlist(sm1_mar$cea$delta_e)), digits = 3)
diff_c_m1 <- round(mean(unlist(sm1_mar$cea$delta_c)), digits = 3)
icer_m1 <- round(sm1_mar$cea$ICER, digits = 3)
summary(sm1_mar)

## ----sumincr, echo=TRUE, eval=FALSE---------------------------------------------------------------
# summary(sm1_mar, incremental = TRUE)

## ----bcea, echo=TRUE, eval=FALSE------------------------------------------------------------------
# BCEA::ceplane.plot(sm1_mar$cea, graph = "ggplot2")
# BCEA::ceac.plot(sm1_mar$cea, graph = "ggplot2")

## ----figplotcea1, echo=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Inspecting probabilistic sensitivity analysis using BCEA built-in functions, for example in terms of cost-effectiveness plane.",out.width='100%'----
cep_sm1_mar <- BCEA::ceplane.plot(sm1_mar$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
ceac_sm1_mar <- BCEA::ceac.plot(sm1_mar$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
cep_sm1_mar

## ----figplotcea2, echo=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Inspecting probabilistic sensitivity analysis using BCEA built-in functions, for example in terms of cost-effectiveness acceptability curve.",out.width='100%'----
ceac_sm1_mar

