library("opencesp")
library("np")
library("boostPM")
library("medicaldata")
library("causaldata")
library("genridge")
library("AppliedPredictiveModeling")
library("openintro")
library("janitor")
library("randomForestSRC")
library("rvinecopulib")
library("arf")
library("RGAN")
library("T4transport")
library("patchwork")
library("grid")
library("ggplot2")
library("dslabs")

wd <- ""
plot_dir <- ""

bor_cols <- function() {
  c("#534785", "#f55f5e", "#dec91a", "#51A3A3", "#82846D", "#A1C084", "#754F5B","#27476E", "#EB9486")
}

save_plot <- function(dir, width=7, height=5.5, p = NULL) {
  files <- list.files(path = dir, pattern = "plot[0-9]+\\.jpg", full.names = FALSE)
  n <- as.numeric(gsub("plot([0-9]+)\\.jpg", "\\1", files))
  if (length(n) > 0) max_n <- max(n, na.rm = TRUE)
  else max_n <- 0
  if(is.null(p)) p <- ggplot2::last_plot()
  filename <- paste(dir, "/plot", max_n + 1, ".jpg", sep="")
  ggplot2::ggsave(file = filename, plot = p,
         width = width, height = height, dpi = 320)
  return(filename)
}

fit_vine_cv <- function(data, max_trunc = 5) {
  max_trunc <- min(max_trunc, ncol(data) - 1)

  folds <- sample(rep(1:5, length.out = nrow(data)))

  scores <- rowMeans(sapply(1:5, function(fold) {
    fit <- vine(data[folds != fold, ], copula_controls = list(family_set = "nonparametric", nonpar_method = "quadratic", trunc_lvl = max_trunc))

    sapply(0:max_trunc, function(trunc_lvl) {
      mean(log(dvine(as.matrix(data[folds == fold, ]), truncate_model(fit, trunc_lvl))))
    })
  }))

  vine(data, copula_controls = list(family_set = "nonparametric", nonpar_method = "quadratic", trunc_lvl = which.max(scores) - 1))
}

sim_nsparse <- function(n, d) {
  r <- chol(0.5^abs(outer(seq_len(d), seq_len(d), "-")))
  as.data.frame(matrix(rnorm(n * d), n, d) %*% solve(t(r)))
}

sim_sparse <- function(n, d) {
  rho <- 0.25 / sqrt(d - 1)
  omega <- diag(d)
  omega[1, -1] <- rho
  omega[-1, 1] <- rho
  r <- chol(omega)
  as.data.frame(matrix(rnorm(n * d), n, d) %*% solve(t(r)))
}

sim_ind <- function(n, d) {
  as.data.frame(matrix(rnorm(n * d), n, d))
}

sim_manifold <- function(n, d) {
  u <- matrix(rnorm(n * 3), n, 3)
  u <- u / sqrt(rowSums(u^2))
  r <- exp(0.05 * rnorm(n))
  if(d > 3) {
    z <- matrix(rnorm(n * (d - 3), sd = 0.05), n, d - 3)
    x <- cbind(u * r, z)
  }
  else x <- u * r
  q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  list(data = as.data.frame(x %*% q), q = q)
}

dens_nsparse <- function(data) {
  r <- chol(0.5^abs(outer(seq_len(ncol(data)), seq_len(ncol(data)), "-")))
  -0.5 * (ncol(data) * log(2 * pi) - 2 * sum(log(diag(r))) + rowSums((as.matrix(data) %*% t(r))^2))
}

dens_sparse <- function(data) {
  d <- ncol(data)
  rho <- 0.25 / sqrt(d - 1)
  omega <- diag(d)
  omega[1, -1] <- rho
  omega[-1, 1] <- rho
  r <- chol(omega)
  -0.5 * (d * log(2 * pi) - 2 * sum(log(diag(r))) + rowSums((as.matrix(data) %*% t(r))^2))
}

dens_ind <- function(data) {
  rowSums(dnorm(as.matrix(data), log = TRUE))
}

dens_manifold <- function(data, q) {
  x <- as.matrix(data) %*% t(q)
  r <- sqrt(rowSums(x[, 1:3, drop = FALSE]^2))
  ldens <- dlnorm(r, 0, 0.05, log = TRUE) - log(4 * pi) - 2 * log(r)
  if(ncol(x) > 3) ldens <- ldens - (ncol(x) - 3) / 2 * log(2 * pi * 0.05^2) - rowSums(x[, 4:ncol(x), drop = FALSE]^2) / (2 * 0.05^2)
  ldens
}

sample_spiral <- function(n) {
  class <- rbinom(n, 1, 0.5)
  theta <- runif(n, 0, pi)

  data.frame(
    class = factor(class),
    x = (1 + theta) * cos(theta + class * pi) + rnorm(n, sd = 0.5),
    y = (1 + theta) * sin(theta + class * pi) + (1 - 2 * class) * 1.2 + rnorm(n, sd = 0.5)
  )
}


######################################


set.seed(1)
spiral <- sample_spiral(1000)

spiral_pgb <- pgb_seq_models(spiral)
spiral_arf <- adversarial_rf(spiral, num_trees = 200, verbose = FALSE)
spiral_arf_params <- forde(spiral_arf, spiral)

spiral_arf_synth <- forge(spiral_arf_params, n_synth = nrow(spiral))
spiral_pgb_synth <- pgb_jsample(spiral, nrow(spiral), models = spiral_pgb)

spiral_plot <- rbind(
  transform(spiral, method = "Original"),
  transform(spiral_arf_synth, method = "ARF"),
  transform(spiral_pgb_synth, method = "a-PGB")
)

spiral_plot$method <- factor(spiral_plot$method, levels = c("Original", "ARF", "a-PGB"))

spiral_plot[c("x", "y")] <- data.frame(x = -spiral_plot$y, y = spiral_plot$x)

ggplot(spiral_plot, aes(x, y, color = class)) +
  geom_point(alpha = 0.2) +
  scale_color_manual(values = bor_cols()) +
  facet_wrap(~ method, nrow = 1) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    strip.text = element_text(size = 14),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

save_plot(plot_dir, 7, 3)


######################################


sim_exp <- function(sim_func, dens_func, n_train, n_test, n_rep, d_vec) {
  res_pgb <- matrix(NA_real_, n_rep, length(d_vec))
  res_kde <- matrix(NA_real_, n_rep, length(d_vec))
  res_vines <- matrix(NA_real_, n_rep, length(d_vec))
  res_arf <- matrix(NA_real_, n_rep, length(d_vec))
  res_boostpm <- matrix(NA_real_, n_rep, length(d_vec))

  time_pgb <- matrix(NA_real_, n_rep, length(d_vec))
  time_kde <- matrix(NA_real_, n_rep, length(d_vec))
  time_vines <- matrix(NA_real_, n_rep, length(d_vec))
  time_arf <- matrix(NA_real_, n_rep, length(d_vec))
  time_boostpm <- matrix(NA_real_, n_rep, length(d_vec))

  for(j in seq_along(d_vec)) {
    print(j)
    for(i in seq_len(n_rep)) {
      print(i)
      sim <- sim_func(n_train + n_test, d_vec[j])
      if(is.list(sim) && !is.data.frame(sim)) {
        train <- sim$data[seq_len(n_train), , drop = FALSE]
        test <- sim$data[n_train + seq_len(n_test), , drop = FALSE]
        true_dens <- dens_func(test, sim$q)
      }
      else {
        train <- sim[seq_len(n_train), , drop = FALSE]
        test <- sim[n_train + seq_len(n_test), , drop = FALSE]
        true_dens <- dens_func(test)
      }
      time_pgb[i, j] <- system.time({
        pred_pgb <- pgb_jdens(train, test, visit.sequence = sample.int(ncol(train)))
      })[["elapsed"]]
      res_pgb[i, j] <- mean(abs(pred_pgb - true_dens))

      time_kde[i, j] <- system.time({
        bw <- npudensbw(dat = train, bwmethod = "normal-reference")
        pred_kde <- log(npudens(bws = bw, edat = test)$dens)
      })[["elapsed"]]
      res_kde[i, j] <- mean(abs(pred_kde - true_dens))

      time_vines[i, j] <- system.time({
        fit_vine <- fit_vine_cv(train)
        pred_vine <- log(dvine(test, fit_vine))
      })[["elapsed"]]
      res_vines[i, j] <- mean(abs(pred_vine - true_dens))

      time_arf[i, j] <- system.time({
        fit_arf <- adversarial_rf(train, num_trees = 200, verbose = FALSE, parallel = FALSE)
        params_arf <- forde(fit_arf, train, parallel = FALSE)
        pred_arf <- lik(params_arf, test, arf = fit_arf, log = TRUE, parallel = FALSE)
      })[["elapsed"]]
      res_arf[i, j] <- mean(abs(pred_arf - true_dens))

      time_boostpm[i, j] <- system.time({
        train_boostpm <- plogis(as.matrix(train))
        test_boostpm <- plogis(as.matrix(test))
        omega <- cbind(rep(0, ncol(train)), rep(1, ncol(train)))
        fit_boostpm <- fit_boostpm(train_boostpm, Omega = omega, progress = "none")
        pred_boostpm <- predict(fit_boostpm, test_boostpm, type = "log_density") + rowSums(log(test_boostpm) + log1p(-test_boostpm))
      })[["elapsed"]]
      res_boostpm[i, j] <- mean(abs(pred_boostpm - true_dens))
    }
  }
  list(error = list(pgb = res_pgb, kde = res_kde, vines = res_vines, arf = res_arf, boostpm = res_boostpm),
       time = list(pgb = time_pgb, kde = time_kde, vines = time_vines, arf = time_arf, boostpm = time_boostpm))
}

plot_sim_results <- function(res_nsparse, res_sparse, res_ind, res_manifold, d_vec, type = c("error", "time")) {
  type <- match.arg(type)

  plot_data <- do.call(rbind, Map(
    function(res, scenario) do.call(rbind, lapply(names(res[[type]]), function(method) {
      data.frame(
        scenario = factor(scenario, levels = c("Independent", "Sparse", "Dense", "Manifold")),
        method = factor(method, levels = c("pgb", "kde", "vines", "arf", "boostpm"), labels = c("a-PGB", "KDE", "Copulas", "ARF", "boostPM")),
        d = d_vec,
        mean = colMeans(res[[type]][[method]], na.rm = TRUE),
        sd = apply(res[[type]][[method]], 2, sd, na.rm = TRUE)
      )
    })),
    list(res_ind, res_sparse, res_nsparse, res_manifold),
    c("Independent", "Sparse", "Dense", "Manifold")
  ))

  ggplot(plot_data, aes(x = d, y = mean, color = method, group = method)) +
    lapply(c("KDE", "Copulas", "ARF", "boostPM", "a-PGB"), function(m) {
      d <- subset(plot_data, method == m)
      list(
        geom_errorbar(data = d, aes(ymin = mean - sd, ymax = mean + sd), width = 4.5, linewidth = 0.7),
        geom_line(data = d, linewidth = 0.8),
        geom_point(data = d, size = 2)
      )
    }) +
    scale_color_manual(
      values = c("a-PGB" = bor_cols()[2], "KDE" = bor_cols()[1], "Copulas" = bor_cols()[3], "ARF" = bor_cols()[4], "boostPM" = bor_cols()[5]),
      breaks = c("KDE", "Copulas", "ARF", "boostPM", "a-PGB")
    ) +
    facet_wrap(~scenario, nrow = 2, scales = "free_y") +
    labs(x = "Dimension", y = if(type == "error") "Mean absolute error" else "Execution time (seconds)", color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom", strip.text = element_text(size = 14))
}

n_train <- 100
n_test <- 100
n_rep <- 20
d_vec <- c(3, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

set.seed(2)
res_nsparse <- sim_exp(sim_nsparse, dens_nsparse, n_train, n_test, n_rep, d_vec)
res_sparse <- sim_exp(sim_sparse, dens_sparse, n_train, n_test, n_rep, d_vec)
res_ind <- sim_exp(sim_ind, dens_ind, n_train, n_test, n_rep, d_vec)
res_manifold <- sim_exp(sim_manifold, dens_manifold, n_train, n_test, n_rep, d_vec)

plot_sim_results(res_nsparse, res_sparse, res_ind, res_manifold, d_vec, "error")
save_plot(plot_dir, 9, 5.5)

plot_sim_results(res_nsparse, res_sparse, res_ind, res_manifold, d_vec, "time")
save_plot(plot_dir, 9, 5.5)


######################################


data("opt")
opt$PID <- NULL
opt$BL.Cortico <- NULL
opt[] <- lapply(opt, function(x) if(length(unique(x[!is.na(x)])) > 50) as.numeric(replace(as.character(x), as.character(x) == ".", NA)) else x)
opt[] <- lapply(opt, function(x) if(length(unique(x[!is.na(x)])) <= 4) factor(x) else x)
opt <- opt[, colSums(is.na(opt)) == 0]

data("peakVO2")
peakVO2[] <- lapply(peakVO2, function(x) if(length(unique(x)) == 2) factor(x) else x)

data("nhefs")
nhefs$id <- NULL
nhefs$seqn <- NULL
nhefs$censored <- NULL
nhefs[] <- lapply(nhefs, function(x) { attributes(x) <- NULL; x })
nhefs[] <- lapply(nhefs, function(x) if(length(unique(x)) <= 5) factor(x) else x)
nhefs$wt82 <- NULL
nhefs$wt71 <- NULL
nhefs$wtloss <- NULL
nhefs <- as.data.frame(nhefs)

data("prostate")
prostate$train <- NULL
prostate[] <- lapply(prostate, function(x) if(is.integer(x)) as.numeric(x) else x)
prostate$svi <- factor(prostate$svi)

data("supraclavicular")
supra <- supraclavicular
rm(supraclavicular)
supra$subject <- NULL
supra$nerve_block_censor <- NULL
supra$alfentanil <- NULL
supra$opioid_total <- NULL

data("blood_storage")
blood <- blood_storage
rm(blood_storage)
blood$AdjRadTherapy <- NULL
blood[] <- lapply(blood, function(x) if(length(unique(x[!is.na(x)])) <= 4) factor(x) else x)

data("AlzheimerDisease")
ad <- predictors
rm(predictors)
ad$diag <- diagnosis
rm(diagnosis)
ad$male <- factor(ad$male)

data("smartpill")
smartpill$Group <- NULL
smartpill$Gender <- NULL
smartpill$Race <- NULL

data("bdims")
bdims$sex <- NULL
bdims <- as.data.frame(bdims)

hepatitis <- readRDS(paste0(wd, "hepatitis.rds"))

df_l <- list(bdims = bdims, smartpill = smartpill, nhefs = nhefs, peakVO2 = peakVO2, opt = opt, supra = supra, prostate = prostate, blood = blood, hepatitis = hepatitis, ad = ad)

df_l <- lapply(df_l, function(df) {
  df <- df[, sapply(df, is.numeric)]
  df <- na.omit(df)
  df[] <- lapply(df, function(x) {
    if(all(x > 0)) x <- log(x)
    x
  })
  df[] <- lapply(df, function(x) {
    if(mean(x %% 1 == 0) > 0.9) {
      noise_sd <- sd(x) * sqrt(1.15^2 - 1)
      x <- x + rnorm(length(x), sd = noise_sd)
    }
    x
  })
  df[] <- lapply(df, function(x) {
    noise_sd <- sd(x) * sqrt(1.15^2 - 1)
    dup <- duplicated(x) | duplicated(x, fromLast = TRUE)
    if(any(dup)) x[dup] <- x[dup] + rnorm(sum(dup), sd = noise_sd)
    x
  })
  df <- clean_names(df)
  df
})

removed <- 0
df_l <- lapply(df_l, function(df) {
  keep <- Reduce(`&`, lapply(df, function(x) {
    q <- quantile(x, c(.25, .75))
    iqr <- q[2] - q[1]
    x >= q[1] - 5 * iqr & x <= q[2] + 5 * iqr
  }))
  removed <<- removed + sum(!keep)
  df[keep, , drop = FALSE]
})
print(removed)

labs <- c(
  smartpill = "SmartPill",
  nhefs = "Nutrition",
  supra = "Supraclavicular",
  prostate = "Prostate",
  peakVO2 = "Peak VO2",
  blood = "Blood storage",
  hepatitis = "Hepatitis",
  ad = "Alzheimer's disease",
  opt = "Obstetrics",
  bdims = "Measurements"
)

df_out <- data.frame(
  dataset = names(df_l),
  n_obs = sapply(df_l, nrow),
  n_vars = sapply(df_l, function(x) ncol(x))
)
df_out[order(df_out$n_obs), ]


real_exp <- function(data_list, n_rep) {
  res_pgb <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))
  res_kde <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))
  res_vines <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))
  res_arf <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))
  res_boostpm <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))

  time_pgb <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))
  time_kde <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))
  time_vines <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))
  time_arf <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))
  time_boostpm <- matrix(NA_real_, n_rep, length(data_list), dimnames = list(NULL, names(data_list)))

  for(j in seq_along(data_list)) {
    print(names(data_list)[j])
    df <- data_list[[j]]

    for(i in seq_len(n_rep)) {
      print(i)
      train_id <- sample.int(nrow(df), floor(nrow(df) / 2))
      train <- df[train_id, , drop = FALSE]
      test <- df[-train_id, , drop = FALSE]

      time_pgb[i, j] <- system.time({
        pred_pgb <- pgb_jdens(train, test, visit.sequence = sample.int(ncol(train)))
      })[["elapsed"]]
      res_pgb[i, j] <- mean(pred_pgb)

      time_kde[i, j] <- system.time({
        bw <- npudensbw(dat = train, bwmethod = "normal-reference")
        pred_kde <- log(npudens(bws = bw, edat = test)$dens)
      })[["elapsed"]]
      res_kde[i, j] <- mean(pred_kde)

      time_vines[i, j] <- system.time({
        fit_vine <- fit_vine_cv(train)
        pred_vine <- log(dvine(test, fit_vine))
      })[["elapsed"]]
      res_vines[i, j] <- mean(pred_vine)

      time_arf[i, j] <- system.time({
        fit_arf <- adversarial_rf(train, num_trees = 200, verbose = FALSE, parallel = FALSE)
        params_arf <- forde(fit_arf, train, parallel = FALSE)
        pred_arf <- lik(params_arf, test, arf = fit_arf, log = TRUE, parallel = FALSE)
      })[["elapsed"]]
      res_arf[i, j] <- mean(pred_arf)

      time_boostpm[i, j] <- system.time({
        train_mean <- colMeans(train)
        train_sd <- apply(train, 2, sd)
        train_boostpm <- plogis(sweep(sweep(as.matrix(train), 2, train_mean, "-"), 2, train_sd, "/"))
        test_boostpm <- plogis(sweep(sweep(as.matrix(test), 2, train_mean, "-"), 2, train_sd, "/"))
        omega <- cbind(rep(0, ncol(train)), rep(1, ncol(train)))
        boost_fit <- fit_boostpm(train_boostpm, Omega = omega, progress = "none")
        pred_boostpm <- predict(boost_fit, test_boostpm, type = "log_density") + rowSums(log(test_boostpm) + log1p(-test_boostpm)) - sum(log(train_sd))
      })[["elapsed"]]
      res_boostpm[i, j] <- mean(pred_boostpm)
    }
  }

  list(log_likelihood = list(pgb = res_pgb, kde = res_kde, vines = res_vines, arf = res_arf, boostpm = res_boostpm),
       time = list(pgb = time_pgb, kde = time_kde, vines = time_vines, arf = time_arf, boostpm = time_boostpm)
  )
}

latex_real_tables <- function(res_real, data_list) {
  slash <- "\\"
  method_labels <- c(pgb = "a-PGB", kde = "KDE", vines = "Copulas", arf = "ARF", boostpm = "boostPM")

  format_cell <- function(mean_value, sd_value, position) {
    mean_text <- sprintf("%.3f", mean_value)
    sd_text <- sprintf("%.3f", sd_value)

    if(position == 1) return(paste0(slash, "cellcolor{yellow!40}", slash, "colorbox{yellow!40}{", slash, "makecell{", slash, "textbf{", mean_text, "}", slash, slash, slash, "textbf{(", sd_text, ")}}}"))
    if(position == 2) return(paste0(slash, "cellcolor{yellow!15}", slash, "colorbox{yellow!15}{", slash, "makecell{", mean_text, slash, slash, "(", sd_text, ")}}"))
    paste0(slash, "makecell{", mean_text, slash, slash, "(", sd_text, ")}")
  }

  build_table <- function(results, maximize) {
    methods <- c(setdiff(names(results), "pgb"), "pgb")
    results <- results[methods]
    datasets <- colnames(results[[1]])
    datasets <- datasets[order(sapply(data_list[datasets], nrow))]

    columns <- paste0("|", paste(rep("c", length(methods) + 1), collapse = "|"), "|")
    header <- paste(c(paste0(slash, "textbf{Dataset}"), paste0(slash, "textbf{", unname(method_labels[methods]), "}")), collapse = " & ")

    lines <- c(
      paste0(slash, "begin{tabular}{", columns, "}"),
      paste0(slash, "hline"),
      paste0(header, " ", slash, slash, " ", slash, "hline")
    )

    for(dataset in datasets) {
      means <- sapply(results, function(x) mean(x[is.finite(x[, dataset]), dataset]))
      sds <- sapply(results, function(x) sd(x[is.finite(x[, dataset]), dataset]))
      ranking <- order(means, decreasing = maximize)
      positions <- match(seq_along(methods), ranking)
      cells <- mapply(format_cell, means, sds, positions, USE.NAMES = FALSE)
      row <- paste(c(unname(labs[dataset]), cells), collapse = " & ")
      lines <- c(lines, paste0(row, " ", slash, slash, " ", slash, "hline"))
    }

    lines <- c(lines, paste0(slash, "end{tabular}"))
    paste(lines, collapse = "\n")
  }

  log_likelihood <- build_table(res_real$log_likelihood, TRUE)
  time <- build_table(res_real$time, FALSE)

  cat(log_likelihood, "\n\n", time, "\n", sep = "")
}

n_rep <- 20

set.seed(3)
res_real <- real_exp(df_l, n_rep)

latex_real_tables(res_real, df_l)

######################################


order_exp <- function(data_list, n_rep, n_rand = 10) {
  order_names <- c(paste0("Random ", seq_len(n_rand)), "Treewidth")
  res <- array(NA_real_, c(n_rep, length(order_names), length(data_list)), dimnames = list(NULL, order_names, names(data_list)))

  for(k in seq_along(data_list)) {
    print(names(data_list)[k])
    df <- data_list[[k]]
    rand_orders <- replicate(n_rand, sample.int(ncol(df)), simplify = FALSE)

    for(i in seq_len(n_rep)) {
      print(i)
      train_id <- sample.int(nrow(df), floor(nrow(df) / 2))
      train <- df[train_id, , drop = FALSE]
      test <- df[-train_id, , drop = FALSE]

      for(j in seq_len(n_rand)) {
        visit_sequence <- rand_orders[[j]]
        pred <- pgb_jdens(train, test, visit.sequence = visit_sequence)
        res[i, j, k] <- mean(pred)
      }

      visit_sequence <- treewidth_order(train)
      pred <- pgb_jdens(train, test, visit.sequence = visit_sequence)
      res[i, n_rand + 1, k] <- mean(pred)
    }
  }

  res
}

plot_order_res <- function(res, res_real) {
  plot_df <- do.call(rbind, lapply(seq_len(dim(res)[3]), function(k) {
    do.call(rbind, lapply(seq_len(dim(res)[2]), function(j) {
      data.frame(
        dataset = dimnames(res)[[3]][k],
        order = dimnames(res)[[2]][j],
        mean = mean(res[, j, k], na.rm = TRUE),
        sd = sd(res[, j, k], na.rm = TRUE)
      )
    }))
  }))

  reference_means <- sapply(res_real$log_likelihood, function(x) apply(x, 2, function(y) mean(y[is.finite(y)])))

  reference_min <- apply(reference_means, 1, function(x) sort(x)[2])
  reference_range <- apply(reference_means, 1, max) - reference_min

  plot_df$mean <- (plot_df$mean - reference_min[plot_df$dataset]) / reference_range[plot_df$dataset]
  plot_df$sd <- plot_df$sd / reference_range[plot_df$dataset]

  dataset_order <- names(sort(tapply(plot_df$mean[plot_df$order != "Treewidth"], plot_df$dataset[plot_df$order != "Treewidth"], mean), decreasing = TRUE))
  plot_df$dataset <- factor(plot_df$dataset, levels = dataset_order)
  plot_df$order <- factor(plot_df$order, levels = c("Treewidth", setdiff(dimnames(res)[[2]], "Treewidth")))

  plot_df$order_type <- factor(ifelse(plot_df$order == "Treewidth", "Optimized", "Random"), levels = c("Optimized", "Random"))
  plot_df$plot_order <- NA

  for(ds in levels(plot_df$dataset)) {
    idx <- which(plot_df$dataset == ds & plot_df$order_type == "Random")
    plot_df$plot_order[idx] <- rank(-plot_df$mean[idx], ties.method = "first")
  }

  plot_df$plot_order[plot_df$order_type == "Optimized"] <- 0
  plot_df$plot_order <- factor(plot_df$plot_order, levels = 0:(dim(res)[2] - 1))

  mean_df <- aggregate(mean ~ dataset, plot_df[plot_df$order_type == "Random", ], mean)

  ggplot(plot_df, aes(x = dataset, y = mean, color = order_type, group = plot_order)) +
    geom_errorbar(data = mean_df, aes(x = dataset, ymin = mean, ymax = mean), inherit.aes = FALSE, width = 1.05, linewidth = 0.75, color = bor_cols()[8]) +
    geom_point(position = position_dodge(width = 0.9), size = 2) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9), width = 0.6, linewidth = 0.75) +
    scale_color_manual(values = bor_cols()[2:1], name = "Factorization order") +
    labs(x = NULL, y = "Normalized mean logarithmic score") +
    scale_x_discrete(labels = labs) +
    guides(x = guide_axis(angle = 45)) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 15, margin = margin(r = 10)),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
}

n_rep <- 20

set.seed(4)
res_order <- order_exp(df_l[sapply(df_l, ncol) < 50], n_rep)

plot_order_res(res_order, res_real)
save_plot(plot_dir, 9, 5.5)


######################################


plot_pairs <- function(synthetic, original) {
  size <- sqrt(ncol(synthetic))
  gap <- 2
  pair_width <- 2 * size
  cell_width <- pair_width + gap
  canvas_width <- 10 * pair_width + 9 * gap

  synth_norm <- rowSums(synthetic^2)
  original_norm <- rowSums(original^2)
  distances <- outer(synth_norm, original_norm, "+") - 2 * synthetic %*% t(original)
  nearest <- max.col(-distances, ties.method = "first")

  canvas <- matrix(0, 10 * size, canvas_width)

  bg <- do.call(rbind, lapply(seq_len(100), function(i) {
    row <- (i - 1) %/% 10
    col <- (i - 1) %% 10
    start <- col * cell_width

    data.frame(
      xmin = c(start + 0.5, start + size + 0.5),
      xmax = c(start + size + 0.5, start + 2 * size + 0.5),
      ymin = row * size + 0.5,
      ymax = (row + 1) * size + 0.5,
      type = factor(c("Synthetic samples", "Nearest original samples"), levels = c("Synthetic samples", "Nearest original samples"))
    )
  }))

  for(i in seq_len(100)) {
    row <- (i - 1) %/% 10
    col <- (i - 1) %% 10
    start <- col * cell_width
    rows <- row * size + seq_len(size)
    synth_cols <- start + seq_len(size)
    original_cols <- start + size + seq_len(size)

    canvas[rows, synth_cols] <- matrix(synthetic[i, ], size, size, byrow = TRUE)
    canvas[rows, original_cols] <- matrix(original[nearest[i], ], size, size, byrow = TRUE)
  }

  dat <- expand.grid(y = seq_len(nrow(canvas)), x = seq_len(ncol(canvas)))
  dat$value <- as.vector(canvas)
  dat$alpha <- dat$value^2

  ggplot() +
    geom_rect(data = bg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type)) +
    scale_fill_manual(
      values = c(
        "Synthetic samples" = colorspace::darken("#27476E", amount = 0.2),
        "Nearest original samples" = colorspace::darken("#754F5B", amount = 0.2)
      ),
      breaks = c("Synthetic samples", "Nearest original samples"),
      name = NULL
    ) +
    geom_raster(data = dat, aes(x, y, alpha = alpha), fill = "white") +
    scale_alpha(range = c(0, 1), guide = "none") +
    scale_y_reverse() +
    coord_equal(expand = FALSE) +
    theme_void() +
    theme(legend.position = "top")
}

set.seed(5)

mnist <- read_mnist()
ids <- which(mnist$train$labels == 3)

old <- seq(0, 1, length.out = 28)
new <- seq(0, 1, length.out = 14)
original <- t(vapply(seq_len(nrow(mnist$train$images[ids, , drop = FALSE])), function(i) {
  img <- matrix(mnist$train$images[ids, , drop = FALSE][i, ] / 255, 28, 28, byrow = TRUE)
  tmp <- t(vapply(seq_len(28), function(r) approx(old, img[r, ], new)$y, numeric(14)))
  img <- vapply(seq_len(14), function(c) approx(old, tmp[, c], new)$y, numeric(14))
  as.vector(t(img))
}, numeric(14^2)))

train <- as.data.frame(original)

active <- vapply(train, function(x) length(unique(x)) > 1, logical(1))
train_active <- train[, active, drop = FALSE]

coords <- expand.grid(col = 1:14, row = 1:14)
visit <- order((coords$row - 7.5)^2 + (coords$col - 7.5)^2, atan2(coords$row - 7.5, coords$col - 7.5))
visit <- match(visit[visit %in% which(active)], which(active))

fit <- pgb_seq_models(train_active, visit.sequence = visit)
synthetic_active <- pgb_jsample(train_active, 100, models = fit, visit.sequence = visit)

synthetic <- matrix(colMeans(train), 100, ncol(train), byrow = TRUE)
synthetic[, active] <- as.matrix(synthetic_active)
synthetic <- pmin(pmax(synthetic, 0), 1)

plot_pairs(synthetic, original)

save_plot(plot_dir, 7, 3.5)
