
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
    cz.list <- list(NA)
  }
  
  #####
  # cz <- 1
  # tcrit <- abs(qt(alpha, df))
  # a0 <- tcrit^2*v["pred_mod1"] + tcrit^2*v["pred_mod1_mod2"]*cz^2 + 2*tcrit^2*cov["pred_mod2.pred_mod1_mod2"]*cz -
  #   b["pred_mod1"]^2 - 2*b["pred_mod1"]*b["pred_mod1_mod2"]*cz - b["pred_mod1_mod2"]^2*cz^2
  # 
  # b0 <- 2*tcrit^2*cov["pred.pred_mod1"] + 2*tcrit^2*cov["pred.pred_mod1_mod2"]*cz + 2*tcrit^2*cov["pred_mod1.pred_mod2"]*cz +
  #   2*tcrit^2*cov["pred_mod2.pred_mod1_mod2"]*cz^2 - 2*b["pred"]*b["pred_mod1"] - 2*b["pred"]*b["mod1_mod2"]*cz -
  #   2*b["pred_mod1"]*b["pred_mod2"]*cz - 2*b["pred_mod2"]*b["pred_mod1_mod2"]*cz^2
  # #print(c(a0=unname(a0), b0=unname(b0)))
  # print(cov); stop()
  #####
  
  if (!is.numeric(data[, moderator[1]]) || length(unique(data[, moderator[1]])) < 3) {
    stop("The moderator must be a continuous variable to generate confidence bands.", call. = FALSE)
  }
  
  tcrit <- abs(qt(alpha, df))
  #Using equations from Bauer & Curran (2005); need to solve a quadratic
  
  sig.region.list <- vector("list", length(cz.list[[1]]))
  
    if (length(moderator) == 1) {
      #Computing functions for the bands
      slope <- function(x, cz) {
        return(b["pred"] + b["pred_mod1"]*x)
      }
      upper.band <- function(x, cz) {
        return(b["pred"] + b["pred_mod1"]*x + tcrit * sqrt(v["pred"] + 2 * x * cov["pred.pred_mod1"] + x^2 * v["pred_mod1"]))
      }
      
      lower.band <- function(x, cz) {
        return(b["pred"] + b["pred_mod1"]*x - tcrit * sqrt(v["pred"] + 2 * x * cov["pred.pred_mod1"] + x^2 * v["pred_mod1"]))
      }
      
      for (i in seq_along(sig.region.list)) {
        roots <- c(NA, NA)
        a0 <- (tcrit^2 * v["pred_mod1"]) - b["pred_mod1"]^2
        b0 <- 2 * (tcrit^2*cov["pred.pred_mod1"] - b["pred"] * b["pred_mod1"])
        c0 <- tcrit^2*v["pred"] - b["pred"]^2
        roots[1] <- (-b0 + sqrt(b0^2 - 4*a0*c0))/(2*a0)
        roots[2] <- (-b0 - sqrt(b0^2 - 4*a0*c0))/(2*a0)
        
        if (all(is.finite(roots))) {
          low.bound <- min(roots)
          high.bound <- max(roots)
          
          mid.point <- mean(c(low.bound, high.bound))
          CI.mid.point <- c(lower.band(mid.point), upper.band(mid.point))
          if (between(0, CI.mid.point)) sig <- "outside"
          else sig <- "inside"
        }
        else {
          low.bound <- high.bound <- NA
          
          sig <- "nowhere"
        }

        sig.region.list[[i]] <- list(bounds = c(lower = low.bound, upper = high.bound),
                                     sig = sig)
      }

      

    }
    else if (length(moderator) == 2) {
      slope <- function(x, cz) {
        return(b["pred"] + b["pred_mod1"]*x + b["pred_mod2"]*cz + b["pred_mod1_mod2"]*cz*x)
      }
      upper.band <- function(x, cz) {
        return((b["pred"]+b["pred_mod1"]*x+b["pred_mod2"]*cz+b["pred_mod1_mod2"]*x*cz) + 
                 tcrit*sqrt(v["pred"] + v["pred_mod1"]*x^2 + v["pred_mod2"]*cz^2 + 
                              v["pred_mod1_mod2"]*x^2*cz^2 + 2*cov["pred.pred_mod1"]*x +
                              2*cov["pred.pred_mod2"]*cz + 2*cov["pred.pred_mod1_mod2"]*x*cz +
                              2*cov["pred_mod1.pred_mod2"]*x*cz + 2*cov["pred_mod1.pred_mod1_mod2"]*x^2*cz + 
                              2*cov["pred_mod2.pred_mod1_mod2"]*x*cz^2))
                              
      }
      
      lower.band <- function(x, cz) {
        return((b["pred"]+b["pred_mod1"]*x+b["pred_mod2"]*cz+b["pred_mod1_mod2"]*x*cz) - 
                 tcrit*sqrt(v["pred"] + v["pred_mod1"]*x^2 + v["pred_mod2"]*cz^2 + 
                              v["pred_mod1_mod2"]*x^2*cz^2 + 2*cov["pred.pred_mod1"]*x +
                              2*cov["pred.pred_mod2"]*cz + 2*cov["pred.pred_mod1_mod2"]*x*cz +
                              2*cov["pred_mod1.pred_mod2"]*x*cz + 2*cov["pred_mod1.pred_mod1_mod2"]*x^2*cz + 
                              2*cov["pred_mod2.pred_mod1_mod2"]*x*cz^2))      
        }
      
      for (i in seq_along(sig.region.list)) {
        #Possibly incorrect; doesn't line up with bands; may be problem with bands
        roots <- c(NA, NA)
        
        cz <- cz.list[[1]][i]
        # a0 <- (tcrit^2 * (v["pred_mod1"] + (2 * cz.list[[1]][i] * cov["pred_mod1.pred_mod1_mod2"]) + (cz.list[[1]][i]^2 * v["pred_mod1_mod2"])) - (b["pred_mod1"] + (b["pred_mod1_mod2"] * cz.list[[1]][i]))^2)
        # b0 <- 2 * ((tcrit^2 * (cov["pred.pred_mod1"] + (cov["pred.pred_mod1_mod2"] * cz.list[[1]][i]) + (cov["pred_mod1.pred_mod2"] * cz.list[[1]][i]) + (cov["pred_mod2.pred_mod1_mod2"] * cz.list[[1]][i]^2))) - ((b["pred"] + b["pred_mod2"] * cz.list[[1]][i]) * (b["pred_mod1"] + (b["pred_mod1_mod2"] * cz.list[[1]][i]))))
        # c0 <- (tcrit^2 * (v["pred"] + (2 * cz.list[[1]][i] * cov["pred.pred_mod1"]) + (cz.list[[1]][i]^2 * v["pred_mod1"])) - (b["pred"] + (b["pred_mod2"] * cz.list[[1]][i]))^2)
        a0 <- tcrit^2*v["pred_mod1"] + tcrit^2*v["pred_mod1_mod2"]*cz^2 + 2*tcrit^2*cov["pred_mod2.pred_mod1_mod2"]*cz -
          b["pred_mod1"]^2 - 2*b["pred_mod1"]*b["pred_mod1_mod2"]*cz - b["pred_mod1_mod2"]^2*cz^2
        
        b0 <- 2*tcrit^2*cov["pred.pred_mod1"] + 2*tcrit^2*cov["pred.pred_mod1_mod2"]*cz + 2*tcrit^2*cov["pred_mod1.pred_mod2"]*cz +
          2*tcrit^2*cov["pred_mod2.pred_mod1_mod2"]*cz^2 - 2*b["pred"]*b["pred_mod1"] - 2*b["pred"]*b["mod1_mod2"]*cz -
          2*b["pred_mod1"]*b["pred_mod2"]*cz - 2*b["pred_mod2"]*b["pred_mod1_mod2"]*cz^2
        
        c0 <- tcrit^2*v["pred"] + tcrit^2*v["pred_mod2"]*cz^2 + 2*tcrit^2*cov["pred.pred_mod2"]*cz - b["pred"]^2 - 
          2*b["pred"]*b["pred_mod2"]*cz - b["pred_mod2"]^2*cz^2
        roots[1] <- (-b0 + sqrt(b0^2 - 4*a0*c0))/(2*a0)
        roots[2] <- (-b0 - sqrt(b0^2 - 4*a0*c0))/(2*a0)
        
        low.bound <- min(roots)
        high.bound <- max(roots)
        
        mid.point <- mean(c(low.bound, high.bound))
        CI.mid.point <- c(lower.band(mid.point, cz.list[[1]][i]), upper.band(mid.point, cz.list[[1]][i]))
        if (!is.finite(mid.point)) sig <- "nowhere"
        else if (between(0, CI.mid.point)) sig <- "outside"
        else sig <- "inside"

        sig.region.list[[i]] <- list(bounds = c(lower = low.bound, upper = high.bound),
                                     sig = sig)

      }
      


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

plot.sig.region <- function(sr, mod.range = NULL, ...) {
  
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
    #d$level <- factor(rep(1, nrow(d)))
    d <- d[between(d$x, xlimits),]
    d.list[[i]] <- d
   
  }
  #d <- data.frame(x = sort(c(sr$bounds, seq(min(c(xlimits[1], sr$bounds)), max(c(xlimits[2],sr$bounds)), length.out = n))))
  D <- do.call("rbind", d.list)
  
  bounds <- numeric(2*length(sr$sig.region.list))
  bounds.colors <- numeric(2*length(sr$sig.region.list))
  cl <- 0
  for (i in seq_along(sr$sig.region.list)) {
    cl <- cl + 1
    bounds[c(2*i-1,2*i)] <- sr$sig.region.list[[i]]$bounds
    bounds.colors[c(2*i-1,2*i)] <- gg_color_hue(cl)
  }
  sp <- ggplot(D, aes(x = x)) + 
    scale_linetype_identity() +
    geom_hline(yintercept = 0, color = "gray20") + 
    geom_vline(xintercept = bounds[between(bounds, xlimits)], color = bounds.colors[between(bounds, xlimits)], linetype = "dashed") +
    labs(x = sr$var.names[["moderator"]][1], y = paste("Slope of", sr$var.names[["predictor"]], "on", sr$var.names[["outcome"]]),
         title = "Confidence Bands") +
    # scale_y_continuous(limits = c(min(unlist(c(0, d[between(d$x, xlimits), !names(d) %in% c("x", "sig")])), na.rm = TRUE),
    #                               max(unlist(c(0, d[between(d$x, xlimits), !names(d) %in% c("x", "sig")])), na.rm = TRUE))) +
    scale_x_continuous(limits = xlimits, expand = c(0, 0.02)) 
  
  if (length(sr$var.names[["moderator"]]) > 1) {
    sp <- sp + geom_line(aes(y = slope, color = level)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = level), alpha = .35) 
  }
  else {
    sp <- sp + geom_line(aes(y = slope)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .35)  
  }
  return(sp)
  
}

print.sig.region <- function(sr, digits = 4,...) {
  cat("Johnson-Neyman Significance Bounds\n-----------------------------\n")
  for (i in seq_along(sr$sig.region.list)) {
    cat(paste0("The significance bounds for the effect of ", sr$var.names[["predictor"]], " on ", sr$var.names[["outcome"]], " are:\n"))
    if (all(!is.na(sr$cz[[1]]))) cat(paste0("At ", sr$var.names[["moderator"]][2], " = ", sr$cz[[1]][i], ":\n\n"))
    cat(paste0(sr$var.names[["moderator"]][1], " = ", round(sr$sig.region.list[[i]]$bounds["lower"], digits), " and ", round(sr$sig.region.list[[i]]$bounds["upper"], digits), "\n\n"))
    cat(paste0("The simple slope is significant *", sr$sig.region.list[[i]]$sig, "* these bounds."))
    if (any(!between(sr$sig.region.list[[i]]$bounds, sr$mod.range))) {
      which <- c("lower", "upper")[!between(sr$sig.region.list[[i]]$bounds, sr$mod.range)]
      cat(paste("\nNote: the", word.list(which), ifelse(length(which) == 1, "bound is", "bounds are"), "outside the range of the moderator in this sample."))
    }
    cat("\n-----------------------------\n")
  }

}