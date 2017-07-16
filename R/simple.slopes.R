simple.slopes <- function(fit, predictor = NULL, moderator = NULL, at.mod.level = NULL, 
                          mod.level.names = NULL, pred.range = NULL, alpha = .05, sig.region = NULL, ...) {
  fit.check(fit)
  o <- fit.prep(fit, predictor, moderator)
  B <- fit.process(o, fit)
  
  for (i in names(o)) {
    assign(i, o[[i]])
  }
  for (i in names(B)) {
    assign(i, B[[i]])
  }
  
  if (length(pred.range) == 0 || !is.numeric(pred.range)) pred.range <- range(data[, predictor])
  
  cz.list <- cz.prep(moderator, at.mod.level = at.mod.level, mod.level.names = mod.level.names, data = data, sig.region = sig.region)
  cz.mat <- setNames(expand.grid(cz.list), paste0("mod", seq_along(moderator)))
  cz.names.mat <- setNames(expand.grid(lapply(cz.list, function(x) if (length(names(x)) > 0) names(x) else rep("", length(x))), stringsAsFactors = FALSE), 
                           paste0("mod", seq_along(moderator)))
  cz <- list(val = cz.mat,
             names = cz.names.mat)
  simple.lines <- vector("list", nrow(cz.mat))
  for (j in seq_along(simple.lines)) {
    #Only correct with one moderator
    # mod <- unlist(lapply(seq_along(moderator), function(x) sapply(combn(paste0("mod", seq_along(moderator)), x, simplify = FALSE), paste, collapse = "_")))
    # w0 <- b["intercept"] + sum(sapply(mod, function(x) b[x]*cz.mat[j, x]))
    if (length(moderator) == 1) {
      if (all(attr(o, "cat.mod") == FALSE)) { #continuous moderator
        w0 <- b[["intercept"]] + b[["mod1"]]*cz.mat[j, "mod1"]
        se.w0 <- sqrt(v["intercept"] + v["mod1"]*cz.mat[j, "mod1"]^2 + 2*cov["intercept.mod1"]*cz.mat[j, "mod1"])
       
        w1 <- b[["pred"]] + b[["pred_mod1"]]*cz.mat[j, "mod1"]
        se.w1 <- sqrt(v["pred"] + 2*cz.mat[j, "mod1"]*cov["pred.pred_mod1"] + cz.mat[j, "mod1"]^2*v["pred_mod1"])

      }
      else { #categorical predictor
        print(cz$names[j,1])
        w0 <- b[["intercept"]][cz$names[j,1]]
        se.w0 <- sqrt(v[["intercept"]][cz$names[j,1]])
        
        w1 <- b[["mod1"]][cz$names[j,1]]
        se.w1 <- sqrt(v[["mod1"]][cz$names[j,1]])
        
      }

    }
    else if (length(moderator) == 2) {
      w0 <- b[["intercept"]] + b[["mod1"]]*cz.mat[j, "mod1"] + b[["mod2"]]*cz.mat[j, "mod2"] + b[["mod1_mod2"]]*cz.mat[j, "mod1"]*cz.mat[j, "mod2"]
      se.w0 <- sqrt(v["intercept"] + v["mod1"]*cz.mat[j, "mod1"]^2 + v["mod2"]*cz.mat[j, "mod2"]^2  + v["mod1_mod2"]*cz.mat[j, "mod1"]^2*cz.mat[j, "mod2"]^2 +
                      2*(cov["intercept.mod1"]*cz.mat[j, "mod1"] + cov["intercept.mod2"]*cz.mat[j, "mod2"] +
                           cov["intercept.mod1_mod2"]*cz.mat[j, "mod1"]*cz.mat[j, "mod2"] + 
                           cov["mod1.mod2"]*cz.mat[j, "mod1"]*cz.mat[j, "mod2"] +
                           cov["mod1.mod1_mod2"]*cz.mat[j, "mod1"]^2*cz.mat[j, "mod2"] +
                           cov["mod2.mod1_mod2"]*cz.mat[j, "mod1"]*cz.mat[j, "mod2"]^2))
      
      w1 <- b[["pred"]] + b[["pred_mod1"]]*cz.mat[j, "mod1"] + b[["pred_mod2"]]*cz.mat[j, "mod2"] + b[["pred_mod1_mod2"]]*cz.mat[j, "mod1"]*cz.mat[j, "mod2"]
      se.w1 <- sqrt(v["pred"] + v["pred_mod1"]*cz.mat[j, "mod1"]^2 + v["pred_mod2"]*cz.mat[j, "mod2"]^2 + v["pred_mod1_mod2"]*cz.mat[j, "mod1"]^2*cz.mat[j, "mod2"]^2 +
                      2*(cov["pred.pred_mod1"]*cz.mat[j, "mod1"] + cov["pred.pred_mod2"]*cz.mat[j, "mod2"] + cov["pred.pred_mod1_mod2"]*cz.mat[j, "mod1"]*cz.mat[j, "mod2"] +
                           cov["pred_mod1.pred_mod2"]*cz.mat[j, "mod1"]*cz.mat[j, "mod2"] + cov["pred_mod1.pred_mod1_mod2"]*cz.mat[j, "mod1"]^2*cz.mat[j, "mod2"] +
                           cov["pred_mod2.pred_mod1_mod2"]*cz.mat[j, "mod1"]*cz.mat[j, "mod2"]^2))
      
    }
    
    for (k in c("w0", "w1", "se.w0", "se.w1")) {
      assign(k, unname(get(k)))
    }
    
    t.w0 <- w0/se.w0
    p.w0 <- 2*pt(abs(t.w0), df, lower.tail = FALSE)
    ci.w0 <- c(w0 - qt(alpha/2, df)*se.w0, w0 + qt(alpha/2, df)*se.w0)
    
    t.w1 <- w1/se.w1
    p.w1 <- 2*pt(abs(t.w1), df, lower.tail = FALSE)
    ci.w1 <- c(w1 - qt(alpha/2, df)*se.w1, w1 + qt(alpha/2, df)*se.w1)
    
    simple.lines[[j]] <- list(w0=w0, se.w0=se.w0, t.w0=t.w0, p.w0=p.w0, ci.w0=ci.w0,
                              w1=w1, se.w1=se.w1, t.w1=t.w1, p.w1=p.w1, ci.w1=ci.w1)
    
  }
  
  var.names <- list(outcome = outcome, predictor = predictor, moderator = moderator)
  
  out <- list(simple.lines = simple.lines,
              var.names = var.names,
              pred.range = pred.range,
              cz = cz)
  class(out) <- "simple.slopes"
  return(out)
  
}

plot.simple.slopes <- function(s, pred.range = NULL, ...) {
  
  if (length(pred.range) > 0) {
    if (length(pred.range) == 2 && is.numeric(pred.range)) {
      s$pred.range <- pred.range
    }
    else {
      warning("The argument to pred.range must a numeric vector of length 2. Ignoring pred.range", call. = FALSE)
    }
  }
  xlimits <- s$pred.range
  
  plot.data.list <- vector("list", length(s$simple.lines))
  if (length(s$var.names$moderator) == 1) {
    x <- seq(min(xlimits), max(xlimits), length.out = 5)
    
    for (j in seq_along(plot.data.list)) {
      plot.data.list[[j]] <- data.frame(x = x, 
                                        y = with(s$simple.lines[[j]], w0 + w1*x),
                                        mod1.level = rep(s$cz[["names"]][j, 1], length(x)))
    }
    plot.data <- do.call("rbind", plot.data.list)
    ggplot(data = plot.data, aes(x=x, y=y, color = mod1.level)) + 
      labs(x=s$var.names[["predictor"]], 
           y=s$var.names[["outcome"]], 
           color = s$var.names[["moderator"]][1],
           linetype = s$var.names[["moderator"]][2],
           title = "Simple Slopes") +
      geom_line(size = 1)
    
  }
  else if (length(s$var.names$moderator) == 2) {
    x <- seq(min(xlimits), max(xlimits), length.out = 5)
    
    for (j in seq_along(plot.data.list)) {
      
      plot.data.list[[j]] <- data.frame(x = x,
                                        y = with(s$simple.lines[[j]], w0 + w1*x),
                                        mod1.level = rep(s$cz[["names"]][j, 1], length(x)),
                                        mod2.level = rep(s$cz[["names"]][j, 2], length(x)))
    }
    plot.data <- do.call("rbind", plot.data.list)
    ggplot(data = plot.data, aes(x=x, y=y, color = mod1.level, linetype = mod2.level)) + 
      labs(x=s$var.names[["predictor"]], 
           y=s$var.names[["outcome"]], 
           color = s$var.names[["moderator"]][1],
           linetype = s$var.names[["moderator"]][2],
           title = "Simple Slopes") +
      geom_line(size = 1)
  }
  
}

print.simple.slopes <- function(s, ...) {
  cat("Simple Slopes\n")
  for (j in seq_along(s$simple.lines)) {
    cat(paste0("\nAt ", word.list(sapply(seq_along(s$var.names[["moderator"]]), 
                                         function(i) paste0(s$var.names[["moderator"]][i], " = ", if (is.numeric(s$cz[["val"]][j, i])) round(s$cz[["val"]][j, i], digits = 4) 
                                                                                                  else s$cz[["val"]][j, i], 
                                                            ifelse(s$cz[["names"]][j, i] != s$cz[["val"]][j, i], paste0(" (", s$cz[["names"]][j, i], ")"), "")))), "...\n"))
    cat("Coefficients:\n")
    d <- with(s$simple.lines[[j]], matrix(c(w0, se.w0, t.w0, p.w0,
                                            w1, se.w1, t.w1, p.w1),
                                          nrow = 2, byrow = TRUE))
    dimnames(d) <- list(c("(Intercept)", s$var.names[["predictor"]]),
                        c("Estimate", "Std.Err", "t value", "Pr(>|t|)"))
    printCoefmat(d, ...)
  }
  
}