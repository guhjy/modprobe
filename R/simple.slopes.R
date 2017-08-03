simple.slopes <- function(fit, predictor = NULL, moderator = NULL, at.mod.level = NULL, 
                          mod.level.names = NULL, alpha = .05, sig.region = NULL, ...) {
  fit.check(fit)
  o <- fit.prep(fit, predictor, moderator)
  B <- fit.process2(o, fit)
  
  for (i in names(o)) {
    assign(i, o[[i]])
  }
  for (i in names(B)) {
    assign(i, B[[i]])
  }
  
  pred.range <- range(data[, predictor])
  
  cz.list <- cz.prep(o, B, at.mod.level = at.mod.level, mod.level.names = mod.level.names, data = data, sig.region = sig.region)
  cz <- list(val = setNames(expand.grid(cz.list[["cz"]]), names(o$vars)[o$vars %in% names(cz.list[["cz"]])]),
             names = setNames(expand.grid(lapply(cz.list[["cz"]], function(x) if (length(names(x)) > 0) names(x) else rep("", length(x))), stringsAsFactors = FALSE), 
                              names(o$vars)[o$vars %in% names(cz.list[["cz"]])]))
  cz$val.exp <- cbind(seq_len(nrow(cz$val)), cz$val) #Create ID variable to preserve order after merge()
  if (length(cz.list$mod.exp.list) > 0) {
    for (i in names(cz$val)) {
      if (i %in% names(cz.list$mod.exp.list)) {
        cz$val.exp <- merge(cz$val.exp, cz.list$mod.exp.list[[i]], by = i, sort = FALSE)
        cz$val.exp <- cz$val.exp[, names(cz$val.exp) != i]
      }
    }
  }
  cz$val.exp <- cz$val.exp[order(cz$val.exp[,1]),-1, drop = FALSE] #To preserve correct order disrupted by merge()
  simple.lines <- vector("list", nrow(cz[["val.exp"]])) 
  
  for (i in seq_along(simple.lines)) {
    a <- matrix(1, ncol = 2, nrow = length(coefs), dimnames = list(names(coefs), c("0", "1")))
    
    #For w0
    a[,"0"][indices.with[["pred"]]] <- 0
    a[,"0"][!(indices.with[["intercept"]] | indices.with[["pred"]] | apply(simplify2array(indices.with[["mod"]]), 1, any))] <- 0 #Not intercept, pred, or mods
    for (m in names(indices.with[["mod"]])) {
      a[,"0"][!indices.with[["pred"]] & indices.with[["mod"]][[m]]] <- a[,"0"][!indices.with[["pred"]] & indices.with[["mod"]][[m]]]*cz$val.exp[i, m]
    }
    
    #For w1
    a[,"1"][!indices.with[["pred"]]] <- 0
    a[,"1"][!(indices.with[["intercept"]] | indices.with[["pred"]] | apply(simplify2array(indices.with[["mod"]]), 1, any))] <- 0 #Not intercept, pred, or mods
    for (m in names(indices.with[["mod"]])) {
      a[,"1"][indices.with[["pred"]] & indices.with[["mod"]][[m]]] <- a[,"1"][indices.with[["pred"]] & indices.with[["mod"]][[m]]]*cz$val.exp[i, m]
    }
    
    w <- drop(t(a) %*% coefs)
    se <- sqrt(diag(t(a) %*% vcov %*% a))
    
    t <- w/se
    p <- 2*pt(abs(t), df, lower.tail = FALSE)
    ci.levels <- c(alpha/2, 1-alpha/2)
    ci <- matrix(w + se %o% qt(ci.levels, df), nrow = length(w), ncol = 2, 
                 dimnames = list(names(w), paste0(100*ci.levels, "%")))
    simple.lines[[i]] <- list(w=w, se=se, t=t, p=p, ci=ci)
  }
  
  var.names <- list(outcome = outcome, predictor = predictor, moderator = moderator)
  
  out <- list(simple.lines = simple.lines,
              var.names = var.names,
              pred.range = pred.range,
              link = link,
              cz = cz)
  
  class(out) <- "simple.slopes"
  return(out)
}

plot.simple.slopes <- function(s, pred.range = NULL, colors = NULL, ...) {
  args <- list(...)
  
  if (length(pred.range) > 0) {
    if (length(pred.range) == 2 && is.numeric(pred.range)) {
      s$pred.range <- pred.range
    }
    else {
      warning("The argument to pred.range must a numeric vector of length 2. Ignoring pred.range", call. = FALSE)
    }
  }
  xlimits <- s$pred.range
  nmods <- length(s$var.names$moderator)
  nlines <- nrow(s$cz[["names"]])
  
  #Color
  nlevels.mod1 <- length(unique(s$cz[["names"]][, 1]))
  if (length(args$colours) > 0) colors <- args$colours
  
  if (length(colors) == 0) {
    colors <- gg_color_hue(nlevels.mod1)
  }
  else {
    if (length(colors) == 1) colors <- rep(colors, nlevels.mod1)
    else if (length(colors) > nlevels.mod1) {
      colors <- colors[seq_len(nlevels.mod1)]
      warning(paste("Only using first", nlevels.mod1, "values in colors."), call. = FALSE)
    }
    else if (length(colors) < nlevels.mod1) {
      warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
      colors <- gg_color_hue(nlevels.mod1)
    }
    
    if (!all(sapply(colors, isColor))) {
      warning("The argument to colors contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
      colors <- gg_color_hue(nlevels.mod1)
    }
  }
  
  plot.data.list <- vector("list", nlines)
  len <- 100
  x <- matrix(c(rep(1, len), seq(min(xlimits), max(xlimits), length.out = len)), ncol = 2)
  
  for (j in seq_along(plot.data.list)) {
    Y <- with(s$simple.lines[[j]], x %*% w)
    # if (s$link == "logit") Y <- plogis(Y)
    # else if (s$link == "probit") Y <- pnorm(Y)
    Y <- make.link(s$link)$linkinv(Y)
    
    d <- as.data.frame(matrix(c(x[,-1], Y), nrow = len, byrow = FALSE))
    plot.data.list[[j]] <- setNames(cbind(d, as.data.frame(simplify2array(lapply(seq_len(nmods), function(q) rep(s$cz[["names"]][j, q], length(x)))))),
                                    c("x", "y", paste0("mod", seq_len(ncol(s$cz[["names"]])), ".level")))
  }
  plot.data <- do.call("rbind", plot.data.list)
  ss <- ggplot(data = plot.data, aes(x=x, y=y, color = mod1.level)) + 
    labs(x=s$var.names[["predictor"]], 
         y=s$var.names[["outcome"]], 
         color = s$var.names[["moderator"]][1],
         title = "Simple Slopes") +
    scale_color_manual(values = colors) 
  if (nmods == 2) {
    ss <- ss + geom_line(aes(linetype = mod2.level), size = 1) +
      labs(linetype = s$var.names[["moderator"]][2])
  }
  else {
    ss <- ss + geom_line(size = 1)
  }
  
  return(ss)
  
}

print.simple.slopes <- function(s, digits = 3, ...) {
  cat("Simple Slopes\n\n")
  cat(paste0("The coefficient for ", s$var.names$predictor, " on ", s$var.names$outcome, ":\n"))
  df <- setNames(cbind(s$cz$names, 
                       do.call("rbind", lapply(s$simple.lines, function(j) c(sapply(j[-5], function(x) x["1"]), sapply(j[5], function(x) x["1",])))),
                       ""),
                 c(s$var.names$moderator, c("Estimate", "Std.Err", "t value", "Pr(>|t|)", "CI LB", "CI UB", "")))
  df[,ncol(df)] <- ifelse(df[,"Pr(>|t|)"] < 0.001, "***", ifelse(df[,"Pr(>|t|)"] < 0.01, "** ", ifelse(df[,"Pr(>|t|)"] < 0.05, "*  ", ifelse(df[,"Pr(>|t|)"] < 0.1, ".  ", "   "))))
  print.data.frame(round_df(df, digits = digits))
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  
}

summary.simple.slopes <- function(s, ...) {
  cat("Simple Slopes\n")
  for (j in seq_along(s$simple.lines)) {
    cat(paste0("\nAt ", word.list(sapply(seq_along(s$var.names[["moderator"]]), 
                                         function(i) paste0(s$var.names[["moderator"]][i], " = ", if (is.numeric(s$cz[["val"]][j, i])) round(s$cz[["val"]][j, i], digits = 4) 
                                                            else s$cz[["val"]][j, i], 
                                                            ifelse(s$cz[["names"]][j, i] != s$cz[["val"]][j, i], paste0(" (", s$cz[["names"]][j, i], ")"), "")))), "...\n"))
    cat("Coefficients:\n")
    d <- with(s$simple.lines[[j]], matrix(c(w, se, t, p),
                                          nrow = 2))
    dimnames(d) <- list(c("(Intercept)", s$var.names[["predictor"]]),
                        c("Estimate", "Std.Err", "t value", "Pr(>|t|)"))
    printCoefmat(d, ...)
  }
}