sig.region <- function(fit, predictor = NULL, moderator = NULL, mod.range = NULL, at.mod2.level = NULL, mod2.level.names = NULL, alpha = .05) {
  fit.check(fit)
  o <- fit.prep(fit, predictor, moderator)
  B <- fit.process(o, fit)
  
  for (i in names(o)) {
    assign(i, o[[i]])
  }
  for (i in names(B)) {
    assign(i, B[[i]])
  }
  
  if (length(moderator) == 2) {
    cz.list <- cz.prep(moderator[2], at.mod.level = at.mod2.level, mod.level.names = mod2.level.names, data = data, sig.region = NULL)
  }
  else {
    cz.list <- list(0)
  }
  
  if (!is.numeric(data[, moderator[1]]) || length(unique(data[, moderator[1]])) < 3) {
    stop("The moderator must be a continuous variable to generate confidence bands.", call. = FALSE)
  }
  
  tcrit <- abs(qt(alpha/2, df))
  #Using equations from Bauer & Curran (2005); need to solve a quadratic
  
  sig.region.list <- vector("list", length(cz.list[[1]]))
  

  if (length(moderator) == 1) {
    Q0 = function(cz) b[["pred"]] 
    Q1 = function(cz) b[["pred_mod1"]] 
    Q2 = function(cz) v["pred"] 
    Q3 = function(cz) 2*cov["pred.pred_mod1"] 
    Q4 = function(cz) v["pred_mod1"]
  }
  else if (length(moderator) == 2) {
    Q0 = function(cz) b[["pred"]] + b[["pred_mod2"]]*cz
    Q1 = function(cz) b[["pred_mod1"]] + b[["pred_mod1_mod2"]]*cz
    Q2 = function(cz) v["pred"] + v["pred_mod2"]*cz^2 + 2*cov["pred.pred_mod2"]*cz
    Q3 = function(cz) 2*cov["pred.pred_mod1"] + 2*cov["pred.pred_mod1_mod2"]*cz + 2*cov["pred_mod1.pred_mod2"]*cz + 2*cov["pred_mod2.pred_mod1_mod2"]*cz^2
    Q4 = function(cz) v["pred_mod1"] + v["pred_mod1_mod2"]*cz^2 + 2*cov["pred_mod1.pred_mod1_mod2"]*cz 
  }
  else {
    stop("Interactions of more than 2 variables are not yet supported.", call. = FALSE)
  }
  
  slope <- function(x, cz) {
    return(Q0(cz) + Q1(cz)*x)
  }
  lower.band <- function(x, cz) {
    return(Q0(cz) + Q1(cz)*x - tcrit*sqrt(Q2(cz) + Q3(cz)*x + Q4(cz)*x^2))
  }
  upper.band <- function(x, cz) {
    return(Q0(cz) + Q1(cz)*x + tcrit*sqrt(Q2(cz) + Q3(cz)*x + Q4(cz)*x^2))
  }
  
  for (i in seq_along(sig.region.list)) {
    roots <- c(NA, NA)
    
    cz <- cz.list[[1]][i]
    
    a0 <- Q1(cz)^2 - tcrit^2*Q4(cz)
    b0 <- 2*Q0(cz)*Q1(cz) - tcrit^2*Q3(cz)
    c0 <- Q0(cz)^2 - tcrit^2*Q2(cz)
    
    if (b0^2 - 4*a0*c0 > 0) {
      roots[1] <- (-b0 + sqrt(b0^2 - 4*a0*c0))/(2*a0)
      roots[2] <- (-b0 - sqrt(b0^2 - 4*a0*c0))/(2*a0)
    }
    else {
      roots <- rep(NA, 2)
    }
    
    low.bound <- min(roots)
    high.bound <- max(roots)
    
    mid.point <- mean(c(low.bound, high.bound))
    CI.mid.point <- c(lower.band(mid.point, cz.list[[1]][i]), upper.band(mid.point, cz.list[[1]][i]))
    if (!is.finite(mid.point)) sig <- NULL
    else if (between(0, CI.mid.point)) sig <- "outside"
    else sig <- "inside"
    
    sig.region.list[[i]] <- list(bounds = c(lower = low.bound, upper = high.bound),
                                 sig = sig)
    
  }
  
  
  if (length(mod.range) == 0 || !is.numeric(mod.range)) mod.range <- range(data[, moderator[1]])
  
  var.names <- list(outcome = outcome, predictor = predictor, moderator = moderator)
  
  out <- list(sig.region.list = sig.region.list,
              slope = slope,
              upper.band = upper.band,
              lower.band = lower.band,
              mod.range = mod.range,
              var.names = var.names,
              cz = cz.list)
  class(out) <- "sig.region"
  return(out)
}

plot.sig.region <- function(sr, mod.range = NULL, colors = NULL, ggplot2 = TRUE, ...) {
  args <- list(...)
  if (length(mod.range) > 0) {
    if (length(mod.range) == 2 && is.numeric(mod.range)) {
      sr$mod.range <- mod.range
    }
    else {
      warning("The argument to mod.range must a numeric vector of length 2. Ignoring mod.range", call. = FALSE)
    }
  }
  
  xlimits <- sort(sr$mod.range)
  n <- 101 #number of samples for plot
  d.list <- vector("list", length(sr$sig.region.list))
  for (i in seq_along(d.list)) {
    d <- data.frame(x = seq(xlimits[1], xlimits[2], length.out = n))
    d$slope <- sr$slope(d$x,  sr$cz[[1]][i])
    d$upper <- sr$upper.band(d$x, sr$cz[[1]][i])
    d$lower <- sr$lower.band(d$x, sr$cz[[1]][i])
    # if (sr$sig.region.list[[i]]$sig == "outside") d$sig <- ifelse(between(d$x, sr$sig.region.list[[i]]$bounds, inclusive = FALSE), "nonsig", "sig")
    # else d$sig <- ifelse(between(d$x, sr$sig.region.list[[i]]$bounds), "nonsig", "sig")
    d$level <- factor(rep(sr$cz[[1]][i], nrow(d)))
    d <- d[between(d$x, xlimits),]
    d.list[[i]] <- d
    
  }
  D <- do.call("rbind", d.list)
  
  #Color
  nlevels <- nlevels(D$level)
  if (length(args$colours) > 0) colors <- args$colours
  
  if (length(colors) == 0) {
    colors <- gg_color_hue(nlevels)
  }
  else {
    if (length(colors) == 1) colors <- rep(colors, nlevels)
    else if (length(colors) > nlevels) {
      colors <- colors[seq_len(nlevels)]
      warning(paste("Only using first", nlevels, "values in colors."), call. = FALSE)
    }
    else if (length(colors) < nlevels) {
      warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
      colors <- gg_color_hue(nlevels)
    }
    
    if (!all(sapply(colors, isColor))) {
      warning("The argument to colors contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
      colors <- gg_color_hue(nlevels)
    }
  }
  
  bounds <- numeric(2*nlevels)
  bounds.colors <- numeric(2*nlevels)
  cl <- 0
  for (i in seq_len(nlevels)) {
    cl <- cl + 1
    bounds[c(2*i-1,2*i)] <- sr$sig.region.list[[i]]$bounds
    bounds.colors[c(2*i-1,2*i)] <- colors[cl]
  }
  
  if (ggplot2) {
    sp <- ggplot(D, aes(x = x)) + 
      geom_line(aes(y = slope, color = level)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = level), alpha = .35) +
      geom_hline(yintercept = 0, color = "gray20") + 
      geom_vline(xintercept = bounds[between(bounds, xlimits)], color = bounds.colors[between(bounds, xlimits)], linetype = "dashed") +
      scale_linetype_identity() +
      scale_x_continuous(limits = xlimits, expand = c(0, 0.02)) +
      scale_fill_manual(values = colors) + 
      scale_color_manual(values = colors) +
      labs(x = sr$var.names[["moderator"]][1], y = paste("Slope of", sr$var.names[["predictor"]], "on", sr$var.names[["outcome"]]),
           title = "Confidence Bands")
      
    if (length(sr$var.names[["moderator"]]) > 1) {
      
        sp <- sp + labs(color = sr$var.names[["moderator"]][2],
             fill = sr$var.names[["moderator"]][2])
    }
    else {
      sp <- sp + theme(legend.position="none")
    }
  }
  return(sp)
}

print.sig.region <- function(sr, digits = 4,...) {
  cat("------- Johnson-Neyman Significance Bounds --------\n")
  cat(paste0("The significance bounds for the effect of ", sr$var.names[["predictor"]], " on ", sr$var.names[["outcome"]], " are:\n"))
  cat("---------------------------------------------------\n")
  for (i in seq_along(sr$sig.region.list)) {
    if (length(sr$var.names[["moderator"]]) == 2) cat(paste0("At ", sr$var.names[["moderator"]][2], " = ", sr$cz[[1]][i], ":\n\n"))
    if (all(!is.na(sr$sig.region.list[[i]]$bounds))) {
      cat(paste0(sr$var.names[["moderator"]][1], " = ", round(sr$sig.region.list[[i]]$bounds["lower"], digits), " and ", round(sr$sig.region.list[[i]]$bounds["upper"], digits), "\n\n"))
      cat(paste0("The simple slope is significant *", sr$sig.region.list[[i]]$sig, "* these bounds.\n"))
      if (any(!between(sr$sig.region.list[[i]]$bounds, sr$mod.range))) {
        which <- c("lower", "upper")[!between(sr$sig.region.list[[i]]$bounds, sr$mod.range)]
        cat(paste("Note: the", word.list(which), ifelse(length(which) == 1, "bound is", "bounds are"), "outside the range of the moderator in this sample.\n"))
      }
    }
    else {
      cat("The simple slope is never significant.\n")
    }
    cat("---------------------------------------------------\n")
  }
  
}