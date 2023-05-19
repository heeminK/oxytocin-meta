# 0. Setup ----
rm(list = ls())
options(scipen = 999)

if (getwd() != "~/Dropbox/PhD/OxyFlex/meta-nonsocial-oxytocin/oxytocin-meta") {
  setwd("~/Dropbox/PhD/OxyFlex/meta-nonsocial-oxytocin/oxytocin-meta")
}

# Load necessary packages
library(metafor)
library(dplyr)
library(tidyr)
#devtools::install_github("dsquintana/metameta")
library(metameta)
install.packages(c("TOSTER", "RoBMA"))
library(TOSTER)
library(RoBMA)
library(tibble)

# Functions to be used in analyses
# Calculates descriptive summary
descriptive <- function(data) {
  tmp_within <- data[data$EXP_DESIGN == "within", ] |>
    distinct(
      STUDY_ID, 
      .keep_all = TRUE
    )
  mean_within <- tmp_within |>
    summarise(
      within_mean_age = mean(as.numeric(MEAN_AGE_OXT), 
                             na.rm = TRUE),
      within_sd_age = sd(as.numeric(MEAN_AGE_OXT), 
                         na.rm = TRUE),
      within_n = sum(N_OXT)
    )
  tmp_between <- data[data$EXP_DESIGN == "between", ] |>
    distinct(STUDY_ID, 
             .keep_all = TRUE
    ) |>
    rbind(data[data$X == 6, ])
  mean_between <- tmp_between |>
    summarise(
      between_mean_age = mean(as.numeric(MEAN_AGE_OXT, MEAN_AGE_PLC), 
                              na.rm = TRUE),
      between_sd_age = sd(as.numeric(MEAN_AGE_OXT, MEAN_AGE_PLC), 
                          na.rm = TRUE),
      between_n = sum(N_OXT, N_PLC)
    )
  
  return(data.frame(mean_within, mean_between))}

# Equivalence test (Based on Lakens, 2017)
equivalence <- function(model_fit, SESOI) {
  et <- TOSTmeta(
    ES = model_fit$b,
    se = model_fit$se,
    low_eqbound_d = (-1) * SESOI,
    high_eqbound_d = SESOI,
    alpha = 0.05,
    plot = FALSE,
    verbose = TRUE
  )
  et_L <- et[['LL_CI_TOST']]
  et_U <- et[['UL_CI_TOST']]
  et <- list(et_L, et_U)
  names(et) = c("et_L", "et_U")
  df = data.frame(as.list(et))
  return(df)
}

# Run equivalence test for each study
equivalence_study <- function(rows, SESOI) {
  et <- TOSTmeta(
    ES = as.numeric(rows[['yi']]),
    se = as.numeric(rows[['vi']]),
    low_eqbound_d = (-1) * SESOI,
    high_eqbound_d = SESOI,
    alpha = 0.05,
    plot = FALSE,
    verbose = FALSE
  )
  et_L <- et[['LL_CI_TOST']]
  et_U <- et[['UL_CI_TOST']]
  et <- list(et_L, et_U)
  names(et) = c("et_L", "et_U")
  df = data.frame(as.list(et))
  return(df)
}

# Moderator analysis
meta_mod <- function(data, mods) {
  df <- rma.mv(yi,
               vi,
               random = ~ 1 | STUDY_ID/X,
               data = data,
               mods = ~ mods)
}

# Making label for forest plot
mlabfun <- function(text, res) {
  W <- diag(1/res$vi)
  X <- model.matrix(res)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  
  I2.b <- 100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P))) # Between: 41.85
  I2.w <- 100 * res$sigma2 / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P))) # Within: 0.00
  
  list(bquote(paste(.(text),
                    " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                    ", df = ", .(res$k - res$p),
                    ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                    "Between-", I^2, " = ", .(formatC(I2.b, digits=2, format="f")), "%, ",
                    "Within-", I^2, " = ", .(formatC(I2.w[2], digits=2, format="f")), "%, ",
                    tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")")))}

# 1. Load data ----
# Read in calculated effect size file and assign it to a variable
dat <- read.csv("ES.csv")

# Excluding columns with missing effect estimate
dat <- dat |>
  filter(!is.na(yi))

# Replacing size of placebo group for within-subject study to 0
dat$N_PLC[is.na(dat$N_PLC)] <- 0

# Reorder rows in the descending order
dat <- dat |> arrange(-yi)
#str(dat)
#View(dat)

# 2. Descriptive statistics ----
summary_data <- descriptive(dat)
summary_data$within_n + summary_data$between_n # total_n = 640
mean(dat[dat$PER_FEMALE != 0 & dat$PER_FEMALE != 100,]$PER_FEMALE) # 48.82%
sd(dat[dat$PER_FEMALE != 0 & dat$PER_FEMALE != 100,]$PER_FEMALE) # SD = 16.85

# 3. Multilevel random-effects meta-analysis ----
# Random-effects meta-analysis clustered by STUDY_ID
res <- rma.mv(
  yi = yi,
  V = vi, 
  random = ~ 1 | STUDY_ID/X, 
  data = dat,
  slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")")
) # ES = 0.10 [-0.05, 0.24] p = 0.21
  # Test for Heterogeneity (Q(df=18) = 30.82, p-val = 0.03)
#summary(res)

# I2-statistics (Codes from https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
W <- diag(1/res$vi)
X <- model.matrix(res)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

I2.b <- 100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P))) # Between: 41.85
I2.w <- 100 * res$sigma2 / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P))) # Within: 0.00

# 3-1.Equivalence test ----
# Equivalence bounds: standardized mean difference = ±0.25,
# the effect size to detect from 20 studies with 20 participants
# each in treatment and placebo control condition with power = 0.80 (α=0.05) 
# (Valentine et al., 2010).
# Convert Hedge's g to Cohen's d before running equivalence test
dat_d <- dat |>
  mutate(ES_d = case_when(
    EXP_DESIGN == "within" ~ yi / (1 - 3 / (4 * (N_OXT + N_OXT) - 9)),
    EXP_DESIGN == "between" ~ yi / (1 - 3 / (4 * (N_OXT + N_PLC) - 9))
  ))

res_d <- rma.mv(
  yi = ES_d,
  V = vi,
  # standard error is the same as Hedges' g
  random = ~ 1 | STUDY_ID/X,
  data = dat_d,
  slab = paste0(gsub(", .*$", "", AUTHORS), " et al. (", YEAR, ")")
) # d = 0.10 [-0.05, 0.25] p = 0.20

# Apply equivalence_study function for calculating 90% CI
for (s in c(0.1, 0.2, 0.25)) {
  # for-loop: equivalence test in a range of SESOIs
  tmp <- as.data.frame(apply(dat_d, 1, equivalence_study, SESOI = s))
  tmp <- matrix(tmp, ncol = 2, byrow = TRUE)
  tmp <- as.data.frame(tmp)

  if (s == 0.1) {
  dat[, ncol(dat) + 1] <- matrix(round(as.numeric(tmp$V1), 2), nrow = 19)
  dat[, ncol(dat) + 1] <- matrix(round(as.numeric(tmp$V2), 2), nrow = 19)
  names(dat)[ncol(dat) - 1] <- paste0("UL_TOST")
  names(dat)[ncol(dat)] <- paste0("LL_TOST")
  }
}

## Equivalence test with preregistered SESOI (g = 0.25)
et_25 <- equivalence(
  model_fit = res_d, 
  SESOI = 0.25
) # TOST: reject null equivalence hypothesis (Z = -1.94 (p = 0.03))

## Equivalence tests with smaller SESOIs (g = 0.1, 0.2)
et_10 <- 
  equivalence(
    model_fit = res_d, 
    SESOI = 0.1
  ) # TOST: don't reject null equivalence hypothesis (Z = -0.01 (p-val = 0.50))

et_20 <- equivalence(
  model_fit = res_d, 
  SESOI = 0.2
) # TOST: don't reject null equivalence hypothesis (Z = -1.30 (p = 0.10))

# 3-2. Forest plot using metafor ----
# Forest Plot with Multiple CIs
png("figures/forest_metafor_CIs0215.png", units = "in", height = 10, width = 15, res = 800)
forest_CIs <- with(dat, metafor::forest(yi, vi, xlim=c(-7,4), ylim=c(-1.5,nrow(dat)+3),
                 psize=1, efac=1, slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")"),
                 ilab = EXP_OUTCOME_MEASURE,
                 ilab.xpos = -4,
                 at = c(-2, -1, 0, 1, 2),
                 xlab = "Hedges' g",
                 order = yi,
                 cex = 1.15,
                 header = "Author(s) (Year)"))
op <- par(cex=1.15, font=2)
text(-4, res$k+2, "Outcome measure")
par(op)
par(new=TRUE)
with(dat, metafor::forest(yi, vi, xlim=c(-7,4), ylim=c(-1.5,nrow(dat)+3), level=90, lwd=3,
                 psize=NA, efac=NA, annotate=FALSE, refline=NA, slab=NA,
                 xaxt="n", xlab="", lty=c("solid","blank"),
                 order = yi))
#addpoly(res, row = -1, level = 95, mlab = mlabfun("Multilevel RE Model", res), annotate = FALSE, border = FALSE)
rect(-0.1, 0, 0.1, 20, 
     col = rgb(red = 0.2, blue = 0.2, green = 0.2, alpha = 0.3),
     border = "transparent")
rect(-0.2, 0, 0.2, 20, 
     col = rgb(red = 0.2, blue = 0.2, green = 0.2, alpha = 0.2),
     border = "transparent")
rect(-0.25, 0, 0.25, 20, 
     col = rgb(red = 0.2, blue = 0.2, green = 0.2, alpha = 0.1),
     border = "transparent")

rect(-0.1, -1.5, 0.1, -0.5, 
     col = rgb(red = 0.2, blue = 0.2, green = 0.2, alpha = 0.3),
     border = "transparent")
rect(-0.2, -1.5, 0.2, -0.5, 
     col = rgb(red = 0.2, blue = 0.2, green = 0.2, alpha = 0.2),
     border = "transparent")
rect(-0.25, -1.5, 0.25, -0.5, 
     col = rgb(red = 0.2, blue = 0.2, green = 0.2, alpha = 0.1),
     border = "transparent")

abline(h=0)

segments(res$ci.lb, -1, res$ci.ub, -1, lty = 1, col = "red", lwd=1.5)
segments(-0.028, -1, 0.227, -1, lty = 1, col = "red", lwd=3)

points(x= res$b, y = -1, col = "red", pch = 23, bg="white", cex = 1.5)
points(x= res$ci.lb, y = -1, col = "red", pch = "|", cex = 1.25)
points(x= res$ci.ub, y = -1, col = "red", pch = "|", cex = 1.25)
text(4, -1, pos = 2, cex=1.15, bquote(paste(
                                            .(formatC(res$b, digits = 2, format = "f")), 
                                           " [", 
                                           .(formatC(res$ci.lb, digits = 2, format = "f")), 
                                           ",  ",
                                           .(formatC(res$ci.ub, digits = 2, format = "f")),
                                           "]")))

W <- diag(1/res$vi)
X <- model.matrix(res)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

I2.b <- 100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P))) # Between: 41.85
I2.w <- 100 * res$sigma2 / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P))) # Within: 0.00

text(-7, -1, pos=4, cex=1.15, bquote(paste("Multilevel RE Model",
                                          " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                                          ", df = ", .(res$k - res$p),
                                          ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                                          "Between-", I^2, " = ", .(formatC(I2.b, digits=2, format="f")), "%, ",
                                          "Within-", I^2, " = ", .(formatC(I2.w[2], digits=2, format="f")), "%, ",
                                          tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")"))
)

dev.off()

# 4. Moderator analysis ----
# png("figures/gross_moderator.png")
# par(mfrow = c(3, 1))
# Between vs. Within
png("figures/fp_design0215.png", units = "in", height = 12, width = 15, res = 800)
fp_design <- with(dat, metafor::forest(yi, vi, ilab = EXP_OUTCOME_MEASURE, ilab.xpos = -3, xlim = c(-6, 4), 
                  ylim = c(-1, 29), at=c(-2, -1, 0, 1, 2),
                  mlab = mlabfun("RE Model for All Studies", res),
                  psize = 1, header = "Author(s) (Year)",
                  rows = c(3:13, 18:25),
                  order = EXP_DESIGN,
                  slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")"),
                  xlab = substitute(paste("Hedges\'", italic("g")))))
### set font expansion factor (as in forest() above) and use a bold font
op <- par(cex = 1, font=2)

### add additional column headings to the plot
text(-3, 28, "Outcome measure")

### switch to bold italic font
par(font = 4)

### add text for the subgroups
text(-6, c(26,14), pos=4, c("Within",
                            "Between"))

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.b <- rma.mv(yi, vi, subset=(EXP_DESIGN=="between"), data=dat, random = ~ 1 | STUDY_ID / X)
res.w <- rma.mv(yi, vi, subset=(EXP_DESIGN=="within"),     data=dat, random = ~ 1 | STUDY_ID / X)

### add summary polygons for the three subgroups
addpoly(res.b, row=16.5, mlab=mlabfun("RE Model for Subgroup", res.b))
addpoly(res.w, row= 1.5, mlab=mlabfun("RE Model for Subgroup", res.w))

### fit meta-regression model to test for subgroup differences
res_design <- meta_mod(dat, dat$EXP_DESIGN)
#(res$sigma2[1] - res_design$sigma2[1]) / res$sigma2[1] # 1% variance can be explained

#op <- par(cex=1, font=2)
#text(-3, res$k+2, "Outcome measure")
#par(op)
par(new=TRUE)
addpoly(res, row = -1, level = 95, annotate = TRUE, border = FALSE)


### add text for the test of subgroup differences
text(-6, -1.8, pos=4, cex=1, bquote(paste("Test for Subgroup Differences: ",
                                          Q[M], " = ", .(formatC(res_design$QM, digits=2, format="f")), ", df = ", .(res_design$p - 1),
                                          ", p = ", .(formatC(res_design$QMp, digits=2, format="f")))))
abline(h = 0)
dev.off()

# intranasal vs other
png("figures/fp_admin0215.png", units = "in", height = 12, width = 15, res = 800)
fp_admin <- with(dat, metafor::forest(yi, vi, ilab = EXP_OUTCOME_MEASURE, ilab.xpos = -3, xlim = c(-6, 4), 
                               ylim = c(-1, 29), at=c(-2, -1, 0, 1, 2),
                               mlab = mlabfun("RE Model for All Studies", res),
                               psize = 1, header = "Author(s) (Year)",
                               rows = c(3:19, 24:25),
                               order = OT_ROUTE,
                               slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")"),
                               xlab = substitute(paste("Hedges\'", italic("g")))))

### set font expansion factor (as in forest() above) and use a bold font
op <- par(cex = 1, font=2)

### add additional column headings to the plot
text(-3, 28, "Outcome measure")

### switch to bold italic font
par(font = 4)

### add text for the subgroups
text(-6, c(26, 20), pos = 4, c("Oral administration",
                              "Intranasal administration"))

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.o <- rma.mv(yi, vi, subset=(OT_ROUTE=="oral"), data=dat, random = ~ 1 | STUDY_ID / X)
res.i <- rma.mv(yi, vi, subset=(OT_ROUTE=="intranasal"),     data=dat, random = ~ 1 | STUDY_ID / X)

### add summary polygons for the three subgroups
addpoly(res.o, row=22.5, mlab=mlabfun("RE Model for Subgroup", res.o))
addpoly(res.i, row= 1.5, mlab=mlabfun("RE Model for Subgroup", res.i))

### fit meta-regression model to test for subgroup differences
res_OT <- meta_mod(dat, dat$OT_ROUTE)
# (res$sigma2[1] - res_OT$sigma2[1]) / res$sigma2[1] # 91% variance of total heterogeneity can be accounted for by ioncluding moderator in the mdel

### add text for the test of subgroup differences
text(-6, -1.8, pos=4, cex=1, bquote(paste("Test for Subgroup Differences: ",
                                          Q[M], " = ", .(formatC(res_OT$QM, digits=2, format="f")), ", df = ", .(res_OT$p - 1),
                                          ", p = ", .(formatC(res_OT$QMp, digits=2, format="f")))))
par(new=TRUE)
addpoly(res, row = -1, level = 95, annotate = TRUE, border = FALSE)
abline(h=0)

dev.off()


# primary vs secondary
png("figures/fp_primary0215.png", units = "in", height = 12, width = 15, res = 800)
fp_primary <- with(dat, metafor::forest(yi, vi, ilab = EXP_OUTCOME_MEASURE, ilab.xpos = -3, xlim = c(-6, 4), 
                               ylim = c(-1, 29), at=c(-2, -1, 0, 1, 2),
                               mlab = mlabfun("RE Model for All Studies", res),
                               psize = 1, header = "Author(s) (Year)",
                               rows = c(3:8, 13:25),
                               order = EX_NONSOCIAL,
                               slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")"),
                               xlab = substitute(paste("Hedges\'", italic("g")))))

### set font expansion factor (as in forest() above) and use a bold font
op <- par(cex = 1, font=2)

### add additional column headings to the plot
text(-3, 28, "Outcome measure")

### switch to bold italic font
par(font = 4)

### add text for the subgroups
text(-6, c(26, 9), pos = 4, c("Primary",
                              "Non-primary"))

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.p <- rma.mv(yi, vi, subset=(EX_NONSOCIAL==TRUE), data=dat, random = ~ 1 | STUDY_ID / X)
res.n <- rma.mv(yi, vi, subset=(EX_NONSOCIAL==FALSE),     data=dat, random = ~ 1 | STUDY_ID / X)

### add summary polygons for the three subgroups
addpoly(res.p, row=11.5, mlab=mlabfun("RE Model for Subgroup", res.p))
addpoly(res.n, row= 1.5, mlab=mlabfun("RE Model for Subgroup", res.n))

### fit meta-regression model to test for subgroup differences
res_PRIMARY <- meta_mod(dat, dat$EX_NONSOCIAL)
#(res$sigma2[1] - res_PRIMARY$sigma2[1]) / res$sigma2[1] # 70%

### add text for the test of subgroup differences
text(-6, -1.8, pos=4, cex=1, bquote(paste("Test for Subgroup Differences: ",
                                          Q[M], " = ", .(formatC(res_PRIMARY$QM, digits=2, format="f")), ", df = ", .(res_PRIMARY$p - 1),
                                          ", p = ", .(formatC(res_PRIMARY$QMp, digits=2, format="f")))))
par(new=TRUE)
addpoly(res, row = -1, level = 95, annotate = TRUE, border = FALSE)
abline(h=0)

dev.off()

# working memory vs flexibility vs inhibition
png("figures/fp_EF0215.png", units = "in", height = 12, width = 15, res = 800)
fp_EF <- with(dat, metafor::forest(yi, vi, ilab = EXP_OUTCOME_MEASURE, ilab.xpos = -3, xlim = c(-6, 4), 
                          ylim = c(-1, 33), at=c(-2, -1, 0, 1, 2),
                          mlab = mlabfun("RE Model for All Studies", res),
                          psize = 1, header = "Author(s) (Year)",
                          rows = c(3:6, 11:17, 22:29),
                          order = EXP_EF_DOMAIN,
                          slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")"),
                          xlab = substitute(paste("Hedges\'", italic("g")))))

### set font expansion factor (as in forest() above) and use a bold font
op <- par(cex = 1, font=2)

### add additional column headings to the plot
text(-3, 32, "Outcome measure")

### switch to bold italic font
par(font = 4)

### add text for the subgroups
text(-6, c(30, 18, 7), pos = 4, c("Working memory",
                                  "Inhibition",
                                  "Flexibility"))

### set par back to the original settings
par(op)

### fit random-effects model in the three subgroups
res.wm <- rma.mv(yi, vi, subset=(EXP_EF_DOMAIN=="working memory"), data=dat, random = ~ 1 | STUDY_ID / X)
res.i <- rma.mv(yi, vi, subset=(EXP_EF_DOMAIN=="inhibition"),     data=dat, random = ~ 1 | STUDY_ID / X)
res.f <- rma.mv(yi, vi, subset=(EXP_EF_DOMAIN=="flexibility"),     data=dat, random = ~ 1 | STUDY_ID / X)

### add summary polygons for the three subgroups
addpoly(res.wm, row=20.5, mlab=mlabfun("RE Model for Subgroup", res.wm))
addpoly(res.i, row= 9.5, mlab=mlabfun("RE Model for Subgroup", res.i))
addpoly(res.f, row= 1.5, mlab=mlabfun("RE Model for Subgroup", res.f))

### fit meta-regression model to test for subgroup differences
res_EF <- meta_mod(dat, dat$EXP_EF_DOMAIN)
# (res$sigma2[1] - res_EF$sigma2[1]) / res$sigma2[1] # 67%

### add text for the test of subgroup differences
text(-6, -1.8, pos=4, cex=1, bquote(paste("Test for Subgroup Differences: ",
                                            Q[M], " = ", .(formatC(res_EF$QM, digits=2, format="f")), ", df = ", .(res_EF$p - 1),
                                            ", p = ", .(formatC(res_EF$QMp, digits=2, format="f")))))
par(new=TRUE)
addpoly(res, row = -1, level = 95, annotate = TRUE, border = FALSE)
abline(h=0)
dev.off()



# Neurotypical vs clinical (Due to lack of subtypes of clinical sample, sensitivity analysis instead)
res_nt <- rma.mv(yi,
                 vi,
                 data = dat[dat$POPULATION_TYPE == "Neurotypical",],
                 random = ~ 1 | STUDY_ID / X) # g = 0.03 [-0.17, 0.24] (p-val = 0.76), 11 estimates

res_cl <- rma.mv(yi,
                 vi,
                 data = dat[dat$POPULATION_TYPE != "Neurotypical",],
                 random = ~ 1 | STUDY_ID / X) # g = 0.11 [-0.05, 0.28] (p-val = 0.18), 8 estimates

# 5. Sensitivity analysis ----
# Filtering data to exclude extracted/requested/calculated data
dat_excluded <-
  dat |>
  filter(MISSING_DATA != TRUE)

# Two-level random effects model with the filtered data
res_excluded <-
  rma.mv(
    yi = yi,
    V = vi,
    random =  ~ 1 | STUDY_ID / X,
    data = dat_excluded
  )
#print(res_excluded) # ES = 0.03 [-0.15, 0.21] (p = 0.72), 10 records

# 6. Publication bias ----
# 6-1. Publication bias assessment using PublicationBias package
pubBias <- svalue(
  yi = dat$yi,
  vi = dat$vi,
  q = 0,
  clustervar = dat$STUDY_ID,
  model = "robust",
  favor.positive = TRUE,
  small = TRUE,
  return.worst.meta = TRUE
) # Didn't yield svalue

corrected_meta(
  yi = dat$yi,
  vi = dat$vi,
  eta = 1.54,
  # Top psy-journals (Mathur & VanderWeele 2020)
  clustervar = dat$STUDY_ID,
  model = "robust",
  favor.positive = TRUE
) # ES = 0.06 (p-val = 0.41)

# 6-2. Publication bias assessment using Robust Bayesian meta-analysis
fit_RoBMA <- RoBMA(
  d = dat$yi,
  n = dat$N_OXT + dat$N_PLC,
  seed = 1,
  model = "PSMA"
  )
summary(fit_RoBMA)
summary(fit_RoBMA, type = "individual")

# BF_10 = 0.16 
# BF_rf = 1.41
# BF_PB = 0.81
# Model-averaged mean -0.02 [-0.25, 0.04]
# Model-averaged heterogeneity 0.10 [-0.25, 0.33]
# PET Mean = 0.02, Median = 1q = 0, 3q = 0.298
# PEESE Mean = 0.03, Median = 1q = 0, 3q = 0.61
# mu=0.033 (RE model)
fit_RoBMA_d <- RoBMA(
  d = dat_d$ES_d,
  n = dat$N_OXT + dat$N_PLC,
  seed = 1,
  model = "PSMA"
)
summary(fit_RoBMA_d)
summary(fit_RoBMA_d, type = "individual")
# BF_10 = 0.16
# BF_rf = 1.46
# BF_PB = 0.97
# Model-averaged mean -0.02 [-0.25, 0.03]
# Model-averaged heterogeneity 0.10 [0, 0.34]
# PET Mean = 0.02, Median = 1q = 0, 3q = 0.33
# PEESE Mean = 0.03, Median = 1q = 0, 3q = 0.57
# mu=0.03 (RE model)

# Funnel plot symmetry 
png(
  'figures/funnel.png',
  units = "in",
  height = 12,
  width = 15,
  res = 800
)
par(mar = c(5, 5, 1, 1))
funnel.rma(
  res,
  level = c(90, 95, 99),
  back = "white",
  shade = c("light grey", "#5C6BC0", "#9FA8DA"),
  #atransf = transf.ztor,
  legend = FALSE,
  xlab = "Hedge's g",
  # at=c(-1, -0.5, 0, 0.5, 1),
  col = case_when(
    dat$MISSING_DATA == TRUE ~ 'red',
    dat$MISSING_DATA == FALSE ~ 'black'
  ),
  cex = 2,
  digits = 2L,
  cex.lab = 2,
  cex.axis = 2
)
legend(
  "topright",
  inset = 0,
  bg = "white",
  legend = c(
    expression("0.10 < p" <= "1.00"),
    expression("0.05 < p" <= "0.10"),
    expression("0.01 < p" <= "0.05"),
    expression("0.00 < p" <= "0.01")
  ),
  bty = "o",
  fill = c("light grey", "#5C6BC0", "#9FA8DA", "white"),
  box.col = "black",
  cex = 2
)
dev.off()

# Egger's regression
metabias(
  x = dat$yi,
  seTE = dat$vi / sqrt(dat$N_OXT + dat$N_PLC),
  method.bias = "Egger"
)

# 7. Power calculation for individual study ----
# Make a data frame for power calculation
d <- data.frame(
  study = paste0(gsub(", .*$", "", dat$AUTHORS), " et al. (", dat$YEAR, ")"),
  yi = dat$yi,
  se = sqrt(dat$vi),
  upper = dat$yi + (1.96 * sqrt(dat$vi)),
  lower = dat$yi - (1.96 * sqrt(dat$vi))
)

# Power calculation
power <- mapower_ul(dat = d,
                    observed_es = 0.0954,
                    name = "Non-social executive function")
power_dat <- power$dat
power_med_dat <- power$power_median_dat
list_power <- list(power_med_dat)
#write.csv(power_dat, file = "tables/studyLevelPower.csv")


# Firepower plot
fp <- firepower(list_power, es = "Hedges' g")
fp_plot <- fp$fp_plot


# Export as .png
png(
  "figures/firepower.png",
  width = 15,
  height = 3,
  unit = "in",
  res = 800
)
fp_plot
dev.off()

# Extracting median power for a range of "true" effect size (observed, g = 0.1 to 1.0)
fpdata <- fp_plot$dat
fpdata <- as_tibble(fpdata)
fpdata

# 8. Visualizing risk of bias----
# install.packages("reactable")
# install.packages("reactablefmtr")
library(reactable)
library(dplyr)
library(tidyr)
library(reactablefmtr)

rob_data_combined <- read.csv("ROB2_combined.csv")

rob_types <- rob_data_combined %>% 
  mutate(
    D1_icons = dplyr::case_when(
      D1 == "Low" ~ "fas fa-circle-plus",
      D1 == "Some concerns" ~ "fas fa-circle-minus",
      D1 == "High" ~ "fas fa-circle-xmark",
      TRUE ~ "other"
    ),
    D1_colors = dplyr::case_when(
      D1 == "Low" ~ "#51b364",
      D1 == "Some concerns" ~ "#f0bd27",
      D1 == "High" ~ "#e03531",
      TRUE ~ "other"
    ),
    DS_icons = dplyr::case_when(
      DS == "Low" ~ "fas fa-circle-plus",
      DS == "Some concerns" ~ "fas fa-circle-minus",
      DS == "High" ~ "fas fa-circle-xmark",
      TRUE ~ "fas fa-hyphen"
    ),
    DS_colors = dplyr::case_when(
      DS == "Low" ~ "#51b364",
      DS == "Some concerns" ~ "#f0bd27",
      DS == "High" ~ "#e03531",
      TRUE ~ "#5A5A5A"
    ),
    D2_icons = dplyr::case_when(
      D2 == "Low" ~ "fas fa-circle-plus",
      D2 == "Some concerns" ~ "fas fa-circle-minus",
      D2 == "High" ~ "fas fa-circle-xmark",
      TRUE ~ "other"
    ),
    D2_colors = dplyr::case_when(
      D2 == "Low" ~ "#51b364",
      D2 == "Some concerns" ~ "#f0bd27",
      D2 == "High" ~ "#e03531",
      TRUE ~ "other"
    ),
    D3_icons = dplyr::case_when(
      D3 == "Low" ~ "fas fa-circle-plus",
      D3 == "Some concerns" ~ "fas fa-circle-minus",
      D3 == "High" ~ "fas fa-circle-xmark",
      TRUE ~ "other"
    ),
    D3_colors = dplyr::case_when(
      D3 == "Low" ~ "#51b364",
      D3 == "Some concerns" ~ "#f0bd27",
      D3 == "High" ~ "#e03531",
      TRUE ~ "other"
    ),
    D4_icons = dplyr::case_when(
      D4 == "Low" ~ "fas fa-circle-plus",
      D4 == "Some concerns" ~ "fas fa-circle-minus",
      D4 == "High" ~ "fas fa-circle-xmark",
      TRUE ~ "other"
    ),
    D4_colors = dplyr::case_when(
      D4 == "Low" ~ "#51b364",
      D4 == "Some concerns" ~ "#f0bd27",
      D4 == "High" ~ "#e03531",
      TRUE ~ "other"
    ),
    D5_icons = dplyr::case_when(
      D5 == "Low" ~ "fas fa-circle-plus",
      D5 == "Some concerns" ~ "fas fa-circle-minus",
      D5 == "High" ~ "fas fa-circle-xmark",
      TRUE ~ "other"
    ),
    D5_colors = dplyr::case_when(
      D5 == "Low" ~ "#51b364",
      D5 == "Some concerns" ~ "#f0bd27",
      D5 == "High" ~ "#e03531",
      TRUE ~ "other"
    ),
    overall_icons = dplyr::case_when(
      Overall == "Low" ~ "fas fa-circle-plus",
      Overall == "Some concerns" ~ "fas fa-circle-minus",
      Overall == "High" ~ "fas fa-circle-xmark",
      TRUE ~ "other"
    ),
    overall_colors = dplyr::case_when(
      Overall == "Low" ~ "#51b364",
      Overall == "Some concerns" ~ "#f0bd27",
      Overall == "High" ~ "#e03531",
      TRUE ~ "other"
    ))
remotes::install_github("kcuilla/reactablefmtr") 
library(reactablefmtr)

rob_table <- rob_types %>% 
  reactable(
    .,
      defaultPageSize = 12,
      columns = list(
      D1_colors = colDef(show = FALSE),
      D1_icons = colDef(show = FALSE),
      DS_colors = colDef(show = FALSE),
      DS_icons = colDef(show = FALSE),
      D2_colors = colDef(show = FALSE),
      D2_icons = colDef(show = FALSE),
      D3_colors = colDef(show = FALSE),
      D3_icons = colDef(show = FALSE),
      D4_colors = colDef(show = FALSE),
      D4_icons = colDef(show = FALSE),
      D5_colors = colDef(show = FALSE),
      D5_icons = colDef(show = FALSE),
      overall_colors = colDef(show = FALSE),
      overall_icons = colDef(show = FALSE),
      

      Design = colDef(
        name = "",
        style = group_merge_sort("Design"),
        align = "center"
      ),
      Study = colDef(
        align = "center", 
        style = list(borderLeft = "1px solid #eaeaea")
      ),
      D1 = colDef(
        align = "center",
        cell = icon_sets(., icon_color_ref = "D1_colors", icon_ref = "D1_icons", icon_position = "over", icon_size = 20)
      ),
      DS = colDef(
        na = "-",
        align = "center",
        cell = icon_sets(., icon_color_ref = "DS_colors", icon_ref = "DS_icons", icon_position = "over", icon_size = 20)
      ),
      D2 = colDef(
        align = "center",
        cell = icon_sets(., icon_color_ref = "D2_colors", icon_ref = "D2_icons", icon_position = "over", icon_size = 20)
      ),
      D3 = colDef(
        align = "center",
        cell = icon_sets(., icon_color_ref = "D3_colors", icon_ref = "D3_icons", icon_position = "over", icon_size = 20)
      ),
      D4 = colDef(
        align = "center",
        cell = icon_sets(., icon_color_ref = "D4_colors", icon_ref = "D4_icons", icon_position = "over", icon_size = 20)
      ),
      D5 = colDef(
        align = "center",
        cell = icon_sets(., icon_color_ref = "D5_colors", icon_ref = "D5_icons", icon_position = "over", icon_size = 20)
      ),
      Overall = colDef(
        align = "center",
        cell = icon_sets(., icon_color_ref = "overall_colors", icon_ref = "overall_icons", icon_position = "over", icon_size = 20)
      )
    )
  ) %>%
  google_font(font_family = "Arimo") %>%
  add_title("Risk of bias", align = "center") #%>%
  #save_reactable_test("rob_table.png")

library(webshot2)
library(htmlwidgets)

saveWidget(rob_table, html_file)
webshot(html_file, "rtableSnapshot.pdf")

library(ggplot2)
library(ggtext)
library(showtext)
# First argument = name in R
# Second argument = path to .otf-file
font_add('fa-reg', 'Font Awesome 6 Free-Regular-400.otf')
font_add('fa-solid', 'Font Awesome 6 Free-Solid-900.otf')

tib <- tibble(
  family = c('firasans', 'lora', 'lobster'),
  x = 0,
  y = seq(0.0, 1, length.out = 3),
  label = c("<span style='font-family:fa-solid'></span>",
            "<span style='font-family:fa-solid'></span>",
            "<span style='font-family:fa-solid'></span>"
            )
  )

p1 <- tib %>%
  ggplot(aes(x, y, label = label)) +
  geom_richtext(family = tib$family, size = 16, hjust = 0, col = c("#e03531", "#f0bd27", '#51b364'), label.colour = NA) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 1.1)) +
  theme_void()

tib <- tibble(
  family = c('firasans', 'lora', 'lobster'),
  x = 0,
  y = seq(0.0, 1, length.out = 3),
  label = c("High risk of bias",
            "Some concerns",
            "Low risk of bias"
  )
)

p2 <- tib %>%
  ggplot(aes(x, y, label = label)) +
  geom_richtext(family = tib$family, size = 12, hjust = 0, col = "black", label.colour = NA) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 1.1)) +
  theme_void()

p3 <- 

