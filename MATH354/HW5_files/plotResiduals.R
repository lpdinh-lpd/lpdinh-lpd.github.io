histoqq <- function(data,
                    shadeDeviation = TRUE,
                    square = FALSE,
                    dist = "norm",
                    height.scale = 1,
                    cartesian = FALSE,
                    ...) {
  require(tidyverse)
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install 'MASS'")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Please install 'scales'")
  
  # ------------------------------------------------------------------
  # 1. Validate Input
  # ------------------------------------------------------------------
  d.name <- paste0("d", dist)
  q.name <- paste0("q", dist)
  
  if (!exists(q.name, mode = "function") || !exists(d.name, mode = "function")) {
    stop(paste0("Invalid distribution specified. Could not find functions '", 
                q.name, "()' or '", d.name, "()'." ))
  }
  
  qdist <- match.fun(q.name)
  ddist <- match.fun(d.name)
  
  countunique <- data |> na.omit() |> unique() |> length()
  if (countunique == 1) stop("There is no variation in your data object.")
  
  npass <- length(data)
  data <- data |> na.omit()
  if (npass > length(data)) {
    warning(paste0("There were ", (npass - length(data)), " missing observation(s). We removed these to make the plot."))
  }
  
  # ------------------------------------------------------------------
  # 2. Scale Data & Estimate Parameters
  # ------------------------------------------------------------------
  df <- tibble(x = as.numeric(data))
  user.params <- list(...)
  
  if (length(user.params) == 0) {
    if (dist == "unif") {
      x1 <- min(df$x, na.rm = TRUE)
      xn <- max(df$x, na.rm = TRUE)
      user.params <- list(
        min = x1 - (xn - x1) / (nrow(df) - 1),
        max = xn + (xn - x1) / (nrow(df) - 1)
      )
    } else if (dist == "binom") {
      x.bar <- mean(df$x, na.rm = TRUE)
      s2 <- var(df$x, na.rm = TRUE)
      max.x <- max(df$x, na.rm = TRUE)
      
      if (!is.na(s2) && s2 < x.bar && s2 > 0) {
        est.size <- round((x.bar^2) / (x.bar - s2))
        est.size <- max(est.size, max.x)
        est.prob <- x.bar / est.size
      } else {
        est.size <- max.x
        est.prob <- x.bar / est.size
      }
      user.params <- list(size = est.size, prob = est.prob)
    } else {
      fit.map <- c("norm" = "normal", "lnorm" = "log-normal", "exp" = "exponential",
                   "logis" = "logistic", "geom" = "geometric", "nbinom" = "negative binomial", 
                   "pois" = "Poisson")
      fit.name <- if (dist %in% names(fit.map)) fit.map[[dist]] else dist
      
      tryCatch({
        fit <- MASS::fitdistr(df$x, fit.name)
        user.params <- as.list(fit$estimate)
      }, error = function(e) {
        warning("MASS::fitdistr failed: ", e$message)
      })
    }
  }
  
  qdist.w <- function(p) do.call(qdist, c(list(p), user.params))
  ddist.w <- function(x) do.call(ddist, c(list(x), user.params))
  
  # ------------------------------------------------------------------
  # 3. QQ Line Geometry
  # ------------------------------------------------------------------
  y.quant <- quantile(df$x, c(0.25, 0.75), na.rm = TRUE)
  x.quant <- qdist.w(c(0.25, 0.75))
  slope <- ifelse(diff(x.quant) == 0, 1, diff(y.quant) / diff(x.quant))
  intercept <- y.quant[1] - slope * x.quant[1]
  
  len <- sqrt(slope^2 + 1)
  vx <- 1 / len
  vy <- slope / len
  nx <- -slope / len
  ny <- 1 / len
  
  # ------------------------------------------------------------------
  # 4. Bins and Density Data
  # ------------------------------------------------------------------
  test_probs <- seq(0.01, 0.99, length.out = 100)
  test_quants <- suppressWarnings(qdist.w(test_probs))
  test_quants <- test_quants[is.finite(test_quants)]
  is.discrete <- length(test_quants) > 0 && all(test_quants %% 1 == 0)
  
  range_x <- max(df$x, na.rm = TRUE) - min(df$x, na.rm = TRUE)
  bin.width <- range_x / nclass.Sturges(df$x)
  
  support.bounds <- qdist.w(c(0, 1))
  lower.bound <- support.bounds[1]
  upper.bound <- support.bounds[2]
  
  bin.min <- min(if (is.finite(lower.bound)) lower.bound else min(df$x, na.rm = TRUE), min(df$x, na.rm = TRUE))
  bin.max <- max(if (is.finite(upper.bound)) upper.bound else max(df$x, na.rm = TRUE), max(df$x, na.rm = TRUE))
  
  if (is.discrete) {
    bin.min <- floor(bin.min) - 0.5
    bin.max <- ceiling(bin.max) + 0.5
    breaks <- seq(bin.min, bin.max, by = 1)
    bin.width <- 1
    
    hist.data <- df |>
      mutate(
        bin.cut = cut(x, breaks = breaks, include.lowest = TRUE),
        bin.idx = as.numeric(bin.cut),
        bin.mid = bin.min + (bin.idx - 0.5) * bin.width
      ) |>
      group_by(bin.mid) |>
      summarise(count = n(), .groups = "drop")
  } else {
    hist_obj <- hist(df$x, plot = FALSE)
    breaks <- hist_obj$breaks
    bin.width <- diff(breaks)[1]
    hist.data <- tibble(bin.mid = hist_obj$mids, count = hist_obj$counts)
  }
  
  if (is.discrete) {
    kde.max <- max(table(df$x) / nrow(df))
  } else {
    dens.fit <- density(df$x)
    kde.max <- max(dens.fit$y, na.rm = TRUE)
  }
  
  hist.max.dens <- max(hist.data$count) / (nrow(df) * bin.width)
  max.empirical.dens <- max(kde.max, hist.max.dens)
  target.max.h <- range_x * height.scale
  density.multiplier <- target.max.h / max.empirical.dens
  bar.scale <- density.multiplier / (nrow(df) * bin.width)
  dens.cap <- max.empirical.dens * 2.5
  true.w <- sqrt((bin.width / slope)^2 + bin.width^2) / 2
  
  # Polygon creation
  hist.polygons <- hist.data |>
    mutate(
      expected.count = ddist.w(bin.mid) * nrow(df) * bin.width,
      deviation = abs(count - expected.count),
      x.c = (bin.mid - intercept) / slope,
      y.c = bin.mid,
      h = count * bar.scale
    ) |>
    slice(rep(1:n(), each = 4)) |>
    mutate(
      point_id = rep(1:4, times = nrow(hist.data)),
      dx = case_when(point_id %in% c(1, 4) ~ -true.w, TRUE ~ true.w),
      dy = case_when(point_id %in% c(1, 2) ~ 0, TRUE ~ h),
      final.x = x.c + dx * vx + dy * nx,
      final.y = y.c + dx * vy + dy * ny
    )
  
  # Densities
  plot.min <- max(min(c(df$x, qdist.w(0.001)), na.rm = TRUE) - bin.width, lower.bound, na.rm = TRUE)
  plot.max <- min(max(c(df$x, qdist.w(0.999)), na.rm = TRUE) + bin.width, upper.bound, na.rm = TRUE)
  
  if (is.discrete) {
    plot.range <- seq(ceiling(plot.min), floor(plot.max), by = 1)
  } else {
    base.seq <- seq(plot.min, plot.max, length.out = 1000)
    if (is.finite(lower.bound) && lower.bound >= plot.min && lower.bound <= plot.max) {
      base.seq <- c(base.seq, lower.bound - 1e-9, lower.bound, lower.bound + 1e-9)
    }
    if (is.finite(upper.bound) && upper.bound >= plot.min && upper.bound <= plot.max) {
      base.seq <- c(base.seq, upper.bound - 1e-9, upper.bound, upper.bound + 1e-9)
    }
    plot.range <- sort(unique(base.seq))
  }
  
  theoretical.density <- data.frame(val = plot.range) |>
    mutate(
      dens = pmin(replace_na(ddist.w(val), 0), dens.cap),
      h = dens * density.multiplier,
      x.base = (val - intercept) / slope - (if (is.discrete) 0.1 * vx else 0),
      y.base = val - (if (is.discrete) 0.1 * vy else 0),
      x = x.base + h * nx,
      y = y.base + h * ny
    )
  
  dens.fit <- density(df$x)
  if (is.discrete) {
    kde.curve <- as.data.frame(table(df$x)) |>
      mutate(val = as.numeric(as.character(Var1)), dens = Freq / nrow(df)) |>
      dplyr::select(val, dens) |>
      mutate(
        h = dens * density.multiplier,
        x.base = (val - intercept) / slope + 0.1 * vx,
        y.base = val + 0.1 * vy,
        x = x.base + h * nx,
        y = y.base + h * ny
      )
  } else {
    kde.curve <- data.frame(val = dens.fit$x, dens = dens.fit$y) |>
      mutate(
        h = dens * density.multiplier,
        x.base = (val - intercept) / slope,
        x = x.base + h * nx,
        y = val + h * ny
      )
  }
  
  # Bounding Box
  theoreticalquantiles <- qdist.w(ppoints(nrow(df)))
  all.coords.x <- c(theoreticalquantiles, hist.polygons$final.x, theoretical.density$x, kde.curve$x)
  all.coords.y <- c(df$x, hist.polygons$final.y, theoretical.density$y, kde.curve$y)
  limits <- range(c(all.coords.x, all.coords.y), na.rm = TRUE)
  
  # ------------------------------------------------------------------
  # 5. Build ggplot
  # ------------------------------------------------------------------
  if (!shadeDeviation) {
    histoqq.plot <- ggplot(df) +
      geom_polygon(data = hist.polygons, aes(x = final.x, y = final.y, group = bin.mid),
                   color = "black", alpha = 0.5)
  } else {
    max.abs.dev <- max(abs(hist.polygons$deviation), na.rm = TRUE)
    color.bound <- max(ceiling(nrow(df) / 5), max.abs.dev)
    histoqq.plot <- ggplot(df) +
      geom_polygon(data = hist.polygons, aes(x = final.x, y = final.y, group = bin.mid, fill = deviation),
                   color = "black", alpha = 0.5, show.legend = FALSE) +
      scale_fill_gradient(low = "grey90", high = "black", limits = c(0, color.bound), oob = scales::squish) +
      annotate("text", label = "(Shading intensity indicates deviations)", color = "grey10",
               x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5, size = 3, fontface = "italic")
  }
  
  if (is.discrete) {
    histoqq.plot <- histoqq.plot +
      geom_segment(data = theoretical.density, aes(x = x.base, y = y.base, xend = x, yend = y, color = "Gaussian-Assumed"), linewidth = 1) +
      geom_segment(data = kde.curve, aes(x = x.base, y = y.base, xend = x, yend = y, color = "Empirical"), linewidth = 1)
  } else {
    histoqq.plot <- histoqq.plot +
      geom_path(data = theoretical.density, aes(x = x, y = y, color = "Gaussian-Assumed"), linewidth = 1) +
      geom_path(data = kde.curve, aes(x = x, y = y, color = "Empirical"), linewidth = 1)
  }
  
  histoqq.plot <- histoqq.plot +
    stat_qq_line(aes(sample = x), distribution = qdist, dparams = user.params, color = "black", linewidth = 1, alpha = 0.75, show.legend = FALSE) +
    stat_qq(aes(sample = x), distribution = qdist, dparams = user.params, alpha = 0.75, show.legend = FALSE) +
    theme_bw() +
    labs(x = "Theoretical Quantiles", y = "Observed Quantiles", color = "") +
    scale_color_manual(values = c("Gaussian-Assumed" = "#0072B2", "Empirical" = "#E69F00"))
  
  # ------------------------------------------------------------------
  # 6. Set Coordinates (UN-SLANTED POLYGONS + PATCHWORK COMPATIBLE)
  # ------------------------------------------------------------------
  histoqq.plot <- histoqq.plot +
    coord_cartesian(xlim = limits, ylim = limits, expand = TRUE) +
    theme(aspect.ratio = 1)
  
  return(histoqq.plot)
}

plotResiduals <- function(mod,
                          shadeDeviation = TRUE,
                          square = FALSE,
                          std = FALSE,
                          plots = c("histoqq", "fittedVSresid")) {
  ##############################################################################
  # PRE PROCESSING
  ##############################################################################
  require(tidyverse)
  require(patchwork)
  all.plots <- c("histogram", "ecdf", "histoqq", "qq", "fittedVSresid")
  plots <- unique(plots)
  plots <- plots[order(match(plots, all.plots))]
  
  if (any(!(plots %in% all.plots))) {
    stop(
      "We currently support 5 plots. A histoqq plot ('histoqq'), a usual qq plot ('qq'), a histogram ('histogram'), a empirical cdf ('ecdf'), or a fitted vs residual plot ('fittedVSresid'). Please ensure your plots specification contains only these values."
    )
  }

  # Extract correct residual scale (sigma vs robust s)
  if (std) {
    mod.sigma <- 1
  } else if (inherits(mod, "rlm")) {
    mod.sigma <- mod$s
  } else {
    mod.sigma <- summary(mod)$sigma
  }
  
  if (std) {
    if (inherits(mod, "rlm")) {
      # Scaled robust residuals u_i = e_i / s
      ggdat <- data.frame(y.hat = fitted(mod), e = residuals(mod) / mod$s)
      Residual.Label <- "Scaled Robust Residual (u)"
    } else {
      ggdat <- data.frame(y.hat = fitted(mod), e = rstudent(mod))  
      Residual.Label <- "Studentized Residual"
    }
  } else if (inherits(mod, "rlm")) {
    # Weighted robust residuals sqrt(w_i) * e_i
    ggdat <- data.frame(y.hat = fitted(mod), e = sqrt(mod$w) * residuals(mod))
    Residual.Label <- "Weighted Robust Residual"
  } else if (!is.null(mod$weights)) {
    # WLS Pearson residuals
    ggdat <- data.frame(y.hat = fitted(mod), e = residuals(mod, type = "pearson"))
    Residual.Label <- "Weighted Residual"
  } else {
    ggdat <- data.frame(y.hat = fitted(mod), e = residuals(mod))
    Residual.Label <- "Residual"
  }
  
  mod.sigma = ifelse(std, 1, summary(mod)$sigma)
  plot.list <- vector(mode = "list", length = length(plots))
  
  ##############################################################################
  # COMPLETE NECESSARY PLOTS
  ##############################################################################
  for (i in 1:length(plots)) {
    curr.plot <- plots[i]
    ############################################################################
    # HISTO-QQ
    ############################################################################
    if (curr.plot == "histoqq") {
      ggdat.gaussian <- data.frame(x = seq(
        min(ggdat$e) - 0.25 * diff(range(ggdat$e)),
        max(ggdat$e) + 0.25 * diff(range(ggdat$e)),
        length.out = 1000
      )) %>%
        mutate(f = dnorm(
          x,
          #ei should have mean zero
          mean = 0,
          #ei should have common variance  (regression estimated)
          sd = mod.sigma
        ))
      
      plot.list[[i]] <- histoqq(
        ggdat$e,
        shadeDeviation = shadeDeviation,
        square = square,
        dist = "norm",
        mean = 0,
        sd = mod.sigma,
        cartesian = (length(plots)>1)
      )
      
      if(any(c("histogram", "ecdf") %in% plots)){
        plot.list[[i]] <- plot.list[[i]]+ guides(color = "none")
      }
      
      ##########################################################################
      # HISTOGRAM
      ##########################################################################
    } else if (curr.plot == "histogram") {
      ggdat.gaussian <- data.frame(x = seq(
        min(ggdat$e) - 0.25 * diff(range(ggdat$e)),
        max(ggdat$e) + 0.25 * diff(range(ggdat$e)),
        length.out = 1000
      )) %>%
        mutate(f = dnorm(
          x,
          #ei should have mean zero
          mean = 0,
          #ei should have common variance  (regression estimated)
          sd = mod.sigma
        ))
      ####################################
      # Create PDF of Residuals Plot"
      ####################################
      breaks <- hist(ggdat$e, plot = FALSE)$breaks
      plot.list[[i]] <- ggplot(ggdat, aes(x = e)) +
        geom_histogram(
          aes(y = after_stat(density)),
          breaks = breaks,
          fill = "lightgrey",
          color = "black"
        ) +
        geom_hline(yintercept = 0) +
        stat_density(aes(x = e, color = "Empirical"),
                     geom = "line", linewidth = 1,
                     position = "identity") +
        geom_line(data = ggdat.gaussian, 
                  aes(x = x, y = f, color = "Gaussian-Assumed"),
                  linewidth = 1) +
        theme_bw() +
        theme(aspect.ratio = 1) +
        xlab(Residual.Label) +
        ylab("Density") +
        labs(color = "") +
        theme(legend.position = "bottom")+
        scale_color_manual( #Okabe-Ito colorblind-friendly palette
          values = c(
            "Gaussian-Assumed" = "#0072B2", # Blue
            "Empirical"        = "#E69F00"  # Orange
          )
        )
      
      ##########################################################################
      # ECDF
      ##########################################################################
    } else if (curr.plot == "ecdf") {
      ####################################
      # Create Data for CDF of Residuals Plot"
      ####################################
      e.cdf.func <- ecdf(ggdat$e)
      e.cdf <- e.cdf.func(sort(ggdat$e))
      ggdat.ecdf <- data.frame(e = sort(ggdat$e), 
                               e.cdf = e.cdf)
      ggdat.gaussian <- data.frame(x = seq(
        min(ggdat$e) - 0.25 * diff(range(ggdat$e)),
        max(ggdat$e) + 0.25 * diff(range(ggdat$e)),
        length.out = 1000
      )) %>%
        mutate(cdf = pnorm(
          x,
          #ei should have mean zero
          mean = 0,
          #ei should have common variance  (regression estimated)
          sd = mod.sigma
        ))
      
      ####################################
      # Create ECDF of Residuals Plot"
      ####################################
      plot.list[[i]] <- ggplot(ggdat.ecdf, aes(x = e)) +
        geom_step(
          aes(y = e.cdf, color = "Empirical"),
          linewidth = 1
        ) +
        geom_line(data=ggdat.gaussian,
                  aes(x=x, y=cdf, color = "Gaussian-Assumed"),
                  linewidth = 1)+
        geom_hline(yintercept = 0) +
        theme_bw() +
        theme(aspect.ratio = 1) +
        xlab(Residual.Label) +
        ylab("Cumulative Probability") +
        labs(color = "") +
        theme(legend.position = "bottom")+
        scale_color_manual( #Okabe-Ito colorblind-friendly palette
          values = c(
            "Gaussian-Assumed" = "#0072B2", # Blue
            "Empirical"        = "#E69F00"  # Orange
          )
        )
      ##########################################################################
      # QQ
      ##########################################################################
    } else if (curr.plot == "qq") {
      plot.list[[i]] <- ggplot(data = ggdat, aes(sample = e)) + #standardize e
        geom_qq() +
        geom_qq_line() +
        theme_bw() +
        theme(aspect.ratio = 1) +
        xlab("Gaussian Quantiles") +
        ylab("Sample Quantiles")
      ##########################################################################
      # FITTED VERSUS RESIDUALS
      ##########################################################################
    } else if (curr.plot == "fittedVSresid"){
      ggdat.out3 <- ggdat %>%
        filter(abs(e) > 3 * mod.sigma)
      
      plot.list[[i]] <- ggplot(data = ggdat, aes(x = y.hat, y = e)) +
        geom_point(shape = 1) +
        geom_hline(yintercept = 0,
                   color = "#CC79A7",
                   linetype = "dashed",
                   linewidth = 1) +
        xlab(bquote("Fitted Values" ~ (hat(Y)))) +
        ylab(Residual.Label) +
        theme_bw() +
        theme(aspect.ratio = 1) +
        geom_hline(
          yintercept = c(-3, 3) * mod.sigma,
          color = "#CC79A7",
          linetype = "dotted",
          linewidth = 1
        ) +
        geom_point(data = ggdat.out3,
                   aes(x = y.hat, y = e),
                   fill = "#D55E00",
                   shape = 21)
    }
  }
  ####################################
  # Print Plots"
  ####################################
  if(length(plots)==1){
    plot.list[[1]]
  }else if(length(plots)==2){
    (plot.list[[1]] | plot.list[[2]]) + 
      plot_layout(guides = "collect") & 
      theme(legend.position = "bottom")
  }else if(length(plots)==3){
    (plot.list[[1]] | plot.list[[2]] | plot.list[[3]]) + 
      plot_layout(guides = "collect")& 
      theme(legend.position = "bottom")
  }else if(length(plots)==4){
    ((plot.list[[1]] | plot.list[[2]])/((plot.list[[3]])|plot.list[[4]])) + 
      plot_layout(guides = "collect")& 
      theme(legend.position = "bottom")
  }else if(length(plots)==5){
design <- "
EEAB
EECD
"
    wrap_plots(plot.list, design = design) + 
      plot_layout(guides = "collect") & 
      theme(legend.position = "bottom")        
  }
}


library(tidyverse)
beta0 <- 1.0
beta1 <- 0.5

mu.epsilon <- 0
sigma.epsilon <- sqrt(0.25)

n <- 50

set.seed(7272)
x <- runif(n = n, min = 0, max = 5)
log_y <- beta0 + beta1 * x + rnorm(n, mean = 0, sd = sigma.epsilon)
y <- exp(log_y)

dat <- tibble(x=x, y=y)

mod <- lm(y ~ x, data=dat)
plotResiduals(mod, plots=c("histoqq", "qq", "histogram", "ecdf", "fittedVSresid")) # 5
plotResiduals(mod, plots=c("histoqq", "qq", "ecdf", "fittedVSresid")) # 4
plotResiduals(mod, plots=c("histoqq", "ecdf", "fittedVSresid")) # 3
plotResiduals(mod)
plotResiduals(mod, plots="histoqq")