fit.prep <- function(fit, predictor = NULL, moderator = NULL) {
  outcome <- all.vars(fit$terms)[1]
  classes <- attr(fit$terms,"dataClasses")
  x.vars <- intersect(attr(fit$terms, "term.labels"), names(classes))
  intercept <- "(Intercept)"
  
  #Check predictor
  blank.predictor <- length(predictor) == 0
  if (blank.predictor) {
    predictor <- x.vars[x.vars %in% names(classes)[classes == "numeric"] & !x.vars %in% moderator][1]
  }
  else {
    if (length(predictor) > 1) {
      warning("The argument to predictor may only be of length 1. Using the first entry instead.", call. = FALSE)
    }
    if (!predictor[1] %in% x.vars) {
      stop(paste(word.list(predictor[1], quote = TRUE), "is not the name of a variable in the model."), call. = FALSE)
    }
  }
  
  #Check moderator
  blank.moderator <- length(moderator) == 0
  if (blank.moderator) {
    moderator <- x.vars[x.vars != predictor][1]
  }
  else {
    if (length(moderator) > 2) {
      stop("Only up to two moderators are supported.", call. = FALSE)
    }
    if (any(!moderator %in% x.vars)) {
      if (any(moderator %in% x.vars)) {
        wl <- word.list(moderator[!moderator %in% x.vars], is.are = TRUE, quotes = TRUE)
        warning(paste0(wl, " not the name", ifelse(attr(wl, "plural"), "s", ""), " of any variable",
                       ifelse(attr(wl, "plural"), "s", ""), " in  the model and will be ignored."), call. = FALSE)
        moderator <- intersect(moderator, x.vars)
      }
      else stop("No names supplied to moderator are the names of variables in the model.", call. = FALSE)
    }
  }
  pred.and.mod <- c(predictor, moderator)
  #Interactions
  interaction <- vector("list", length(moderator))
  for (i in seq_along(interaction)) {
    possible.interactions <- apply(gtools::permutations(length(pred.and.mod), i+1, pred.and.mod), 1, paste, collapse = ":")
    interaction[[i]] <- intersect(possible.interactions, attr(fit$terms, "term.labels"))
    counts <- setNames(sapply(pred.and.mod, function(x) sum(grepl(x, interaction[[i]], fixed = TRUE))),
                       pred.and.mod)
    if (any(counts < choose(length(pred.and.mod) - 1, i))) stop("Not all interactions are present.", call. = FALSE)
  }

  if (classes[predictor] != "numeric") {
    stop(paste(predictor, "must be a numeric variable."), call. = FALSE)
  }
  
  xv <- c("pred", paste0("mod", seq_along(moderator)))
  bnames <- unlist(lapply(seq_along(xv), 
                                 function(x) sapply(combn(xv, 
                                                          x, simplify = F), 
                                                    paste, collapse="_")))
  bnames <- c("intercept", bnames)
  
  ordered.ints <- vector("list", length(interaction))
  for (I in seq_along(ordered.ints)) {
    nway <- I + 1
    combs <- combn(c(predictor, moderator), nway, simp = F)
    ordered.ints[[I]] <- character(choose(length(xv), nway))
    for (i in seq_along(ordered.ints[[I]])) {
      ordered.ints[[I]][i] <- interaction[[I]][sapply(interaction[[I]], function(x) all(sapply(combs[[i]], grepl, x, fixed = TRUE)))]
    }
  }
  
  vars <- setNames(c(intercept, predictor, moderator, unlist(ordered.ints)),
                   bnames)
  
  out <- list(outcome = outcome,
              intercept = intercept,
              predictor = predictor,
              moderator = moderator,
              interaction = interaction,
              vars = vars)
  attr(out, "cat.mod") <- setNames(classes[moderator] != "numeric",
                                   moderator)
  return(out)
}

fit.process <- function(o, fit) {
  
  coefs <- coef(fit)
  vc <- vcov(fit)
  
  if (all(attr(o, "cat.mod") == FALSE)) {
    b <- setNames(as.list(coefs[o$vars]),
                  names(o$vars))
    
    v <- setNames(diag(vc)[o$vars], names(o$vars))
    
    covcomb <- combn(names(o$vars), 2, simplify = FALSE)
    cov <- setNames(sapply(seq_along(covcomb), function(i) vc[o$vars[covcomb[[i]][1]], 
                                                              o$vars[covcomb[[i]][2]]]),
                    sapply(covcomb, paste, collapse = "."))

  }
  else {
    coefnames <- paste0(o$moderator, fit$xlevels[[o$moderator]])
    missing.coef <- coefnames[!coefnames %in% names(coefs)]
    
    b <- setNames(vector("list", 1 + length(o$moderator)),
                  c("intercept", paste0("mod", seq_along(o$moderator))))
    
    b[["intercept"]] <- setNames(c(coefs[o$intercept], coefs[coefnames[coefnames != missing.coef]] + coefs[o$intercept]),
                   c(missing.coef, coefnames[coefnames != missing.coef]))
    b[["mod1"]] <- setNames(c(coefs[o$predictor], coefs[paste(o$predictor, coefnames[coefnames != missing.coef], sep = ":")] + coefs[o$predictor]),
                   c(missing.coef, coefnames[coefnames != missing.coef]))
    
    variances <- diag(vc)
    v <- setNames(vector("list", 1+ length(o$moderator)),
                  c("intercept", paste0("mod", seq_along(o$moderator))))
    v[["intercept"]] <- setNames(c(variances[o$intercept], variances[coefnames[coefnames != missing.coef]] + variances[o$intercept] + 2*vc[o$intercept, coefnames[coefnames != missing.coef]]),
                  c(missing.coef, coefnames[coefnames != missing.coef]))
    v[["mod1"]] <- setNames(c(variances[o$predictor], variances[paste(o$predictor, coefnames[coefnames != missing.coef], sep = ":")] + variances[o$predictor] + 2*vc[o$predictor, paste(o$predictor, coefnames[coefnames != missing.coef], sep = ":")]),
                           c(missing.coef, coefnames[coefnames != missing.coef]))
    
    cov <- NULL
  }
  

  df <- fit$df.residual
  
  data <- fit$model
  
  out <- list(data = data,
              df = df,
              b = b,
              v = v, 
              cov = cov
  )
  return(out)
}

cz.prep <- function(moderator, at.mod.level = NULL, mod.level.names = NULL, data = NULL, sig.region = NULL) {
  cz <- vector("list", length(moderator))
  if (!is.list(at.mod.level)) {
    at.mod.level <- list(at.mod.level)
  }
  for (i in seq_along(cz)) {
    if (length(at.mod.level[[i]]) > 0) {
      if (is.numeric(at.mod.level[[i]])) {
        if (is.numeric(data[, moderator[i]])) {
          cz[[i]] <- at.mod.level[[i]]
        }
        else if (is.logical(data[, moderator[i]])) {
          if (!all(at.mod.level[[i]] %in% c(0, 1))) {
            warning(paste(moderator[i], "is logical, but the argument to at.mod.level contains values other than 0 or 1. Setting at.mod.level to c(FALSE, TRUE)."))
            cz[[i]] <- c(FALSE, TRUE)
          }
        }
        else if (is.factor(data[, moderator[i]])) {
          cz[[i]] <- levels(data[, moderator[i]])[at.mod.level[[i]]]
        }
        else {
          stop("The argument to at.mod.level is numeric but the moderator is not.", call. = FALSE)
        }
      }
      else if (is.character(at.mod.level[[i]])) {
        if (is.character(data[, moderator[i]]) || is.factor(data[, moderator[i]])) {
          cz[[i]] <- at.mod.level[[i]]
        }
        else {
          stop("The argument to at.mod.level is a character but the moderator is not a character or factor.", call. = FALSE)
        }
      }
      else if (is.logical(at.mod.level[[i]])) {
        if (is.logical(data[, moderator[i]])) {
          cz[[i]] <- at.mod.level[[i]]
        }
        else {
          stop("The argument to at.mod.level is a logical but the moderator is not logical.", call. = FALSE)
        }
      }
      else {
        stop(paste("The argument to at.mod.level must be the same type as", moderator[i]), call. = FALSE)
      }
      if (length(mod.level.names) > 0) {
        if (!is.list(mod.level.names)) {
          mod.level.names <- list(mod.level.names)
        }
        if (length(mod.level.names[[i]]) == length(cz[[i]])) {
          if (length(unique(mod.level.names[[i]])) != length(mod.level.names[[i]])) {
            warning("mod.level.names contains duplicate values. Ignoring mod.level.names.", call. = FALSE)
          }
          else names(cz[[i]]) <- mod.level.names[[i]]
        }
        else {
          warning("mod.level.names is not the same length as at.mod.level. Ignoring mod.level.names.", cal. = FALSE)
        }
      }
    }
    if (length(sig.region) > 0) {
      if (length(cz[[i]]) > 0) {
        warning("An argument to sig.region was specified, but will be ignored. Using at.model.level instead.", call. = FALSE)
      }
      else {
        if (length(sig.region$bounds) == 2 && is.numeric(sig.region$bounds)) {
          cz[[i]] <- sig.region$bounds
          names(cz[[i]]) <- c("Lower Bound", "Upper Bound")
        }
        else {
          warning("An argument to sig.region was specified, but it doesn't contain significance bounds, and will be ignored.", call. = FALSE)
        }
      }
    }
    if (length(cz[[i]]) == 0) {
      if (is.numeric(data[, moderator[i]]) && length(unique(data[, moderator[i]])) > 3) {
        mod.mean <- mean(data[, moderator[i]])
        mod.sd <- sd(data[, moderator[i]])
        cz[[i]] <- setNames(c(mod.mean - mod.sd, mod.mean, mod.mean + mod.sd),
                            c("Mean - SD", "Mean", "Mean + SD"))
      }
      else if (is.factor(data[, moderator[i]])) {
        cz[[i]] <- levels(data[, moderator[i]])
      }
      else {
        cz[[i]] <- sort(unique(data[, moderator[i]]))
      }
    }
    if (length(names(cz[[i]])) == 0) {
      names(cz[[i]]) <- cz[[i]]
    }
  }
  return(cz)
}

between <- function(x, range, inclusive = TRUE) {
  if (length(range) != 2) stop("range must be of length 2.")
  if (any(is.na(c(x, range)))) return(FALSE)
  range <- sort(range)
  if (inclusive) return(x >= range[1] & x <= range[2])
  else return(x > range[1] & x < range[2])
}

word.list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  if (quotes) word.list <- sapply(word.list, function(x) paste0("\"", x, "\""))
  if (L == 0) {
    out <- ""
    attr(out, "plural") = FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") = FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") = FALSE
    }
    else {
      and.or <- match.arg(and.or)
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or," "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "), 
                     word.list[L], sep = paste0(", ", and.or," "))
        
      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") = TRUE
    }
    
    
  }
  return(out)
}

fit.check <- function(fit) {
  classes <- attr(fit$terms,"dataClasses")
  x.vars <- intersect(attr(fit$terms, "term.labels"), names(classes))
  
  if (!all(class(fit) == "lm")) {
    if (any(class(fit) == "glm")) {
      if (fit$family != "gaussian" || fit$link != "identity") {
        stop("The input to fit must be an lm object or a glm object with family = \"guassian\".", call. = FALSE)
      }
    }
    else stop("The input to fit must be a lm object or a glm object with family = \"guassian\".", call. = FALSE)
  }
  if (length(all.vars(fit$terms)) < 3) {
    stop("The regresison needs to contain at least 2 predictors and an outcome.", call. = FALSE)
  }
  if (all(classes[names(classes) %in% x.vars] != "numeric")) {
    stop("No variables in the model are numeric. The predictor must be numeric, but the moderator may not be.", call. = FALSE)
  }
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 50, c = 100)[1:n]
}

# Make coeffs for factors 