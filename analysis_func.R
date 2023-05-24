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
  
  result <- bquote(paste(.(text),
                    " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                    ", df = ", .(res$k - res$p),
                    ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                    "Between-", I^2, " = ", .(formatC(I2.b, digits=2, format="f")), "%, ",
                    "Within-", I^2, " = ", .(formatC(I2.w[2], digits=2, format="f")), "%, ",
                    tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")"))
  
  return(result)
}

est_lab <- function(res) {
  result <- bquote(paste(.(formatC(
    res$b, digits = 2, format = "f"
  )),
  " [",
  .(formatC(
    res$ci.lb, digits = 2, format = "f"
  )),
  ",  ",
  .(formatC(
    res$ci.ub, digits = 2, format = "f"
  )),
  "]"))
  
  return(result)
}

specify_rows <- function(data, sublevel) {
  row_startnr <- c(3)
  row_endnr <- c()
  rows <- c()
  for (i in 1:length(c(levels(factor(data[[sublevel]]))))) {
    n = count(filter(data, data[[sublevel]] == c(levels(factor(data[[sublevel]])))[i]))
    if (length(row_startnr) < length(c(levels(factor(data[[sublevel]]))))) {
      row_startnr <- append(as.integer(row_startnr), n + as.integer(row_startnr[i]) + 4)
    }
    row_endnr <- append(row_endnr, n + as.integer(row_startnr[i]) - 1)
    rows <- append(rows, as.integer(row_startnr[i]):as.integer(row_endnr[i]))
  }
  df <- list(row_startnr = row_startnr, 
             row_endnr = row_endnr,
             rows = rows)
  return(df)
}

  

generate_mod_fp <- function(data, sublevel, res) {
  sublevels <- c(levels(factor(data[[sublevel]])))
  row_info <- specify_rows(data, sublevel)
  fp <- with(data, metafor::forest(yi, vi, ilab = EXP_OUTCOME_MEASURE, ilab.xpos = -3, xlim = c(-6, 4), 
                                  ylim = c(-3, 21 + 4 * length(sublevels)), at=c(-2, -1, 0, 1, 2),
                                  psize = 1, header = "Author(s) (Year)",
                                  rows = row_info$rows,
                                  order = data[[sublevel]],
                                  slab = paste0(gsub( ", .*$", "", AUTHORS), " et al. (", YEAR, ")"),
                                  xlab = substitute(paste("Hedges\'", italic("g")))))
  
  ### set font expansion factor (as in forest() above) and use a bold font
  op <- par(cex = 1, font=2)
  
  ### add additional column headings to the plot
  text(-3, 20 + 4 * length(sublevels), "Outcome measure")
  
  ### switch to bold italic font
  par(font = 4)
  
  ### add text for the subgroups
  for (i in 1:length(sublevels)) {
    text(-6, as.integer(row_info$row_endnr[i]) + 1, pos = 4, sublevels[i])
  }
  
  ### set par back to the original settings
  par(op)
  
  ### fit random-effects model in the subgroups
  ind <- 1
  for (i in sublevels) {
    res <- rma.mv(yi, 
                  vi, 
                  subset=(data[[sublevel]]==i), 
                  data = data, 
                  random = ~ 1 | STUDY_ID / X)
    addpoly(res, 
            row = as.integer(row_info$row_startnr[ind]) - 1.5, 
            mlab=list(mlabfun("RE Model for Subgroup", res)))
    ind = ind + 1
  }
  
  ### fit meta-regression model to test for subgroup differences
  mod_res <- meta_mod(data, data[[sublevel]])
  # (res$sigma2[1] - res_EF$sigma2[1]) / res$sigma2[1] # 67%
  
  ### add text for the test of subgroup differences
  text(-6, -2.5, pos=4, cex=1, bquote(paste("Test for Subgroup Differences: ", Q[M], 
                                            " = ", .(formatC(mod_res$QM, digits=2, format="f")), 
                                            ", df = ", .(mod_res$p - 1),
                                            ", p = ", .(formatC(mod_res$QMp, digits=2, format="f")))))
  par(new=TRUE)
  addpoly(res, row = -1, level = 95, annotate = TRUE, border = FALSE)
  abline(h=0)
  
  return(fp)
}  

# generate_mod_fp(dat, "EXP_DESIGN", res)
