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
    stop(paste("The predictor must be a numeric variable."), call. = FALSE)
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
    coefnames <- setNames(paste0(o$moderator, fit$xlevels[[o$moderator]]),
                          fit$xlevels[[o$moderator]])
    missing.coef <- coefnames[!coefnames %in% names(coefs)]
    
    b <- setNames(vector("list", 1 + length(o$moderator)),
                  c("intercept", paste0("mod", seq_along(o$moderator))))
    
    b[["intercept"]] <- setNames(c(coefs[o$intercept], coefs[coefnames[coefnames != missing.coef]] + coefs[o$intercept]),
                   c(names(missing.coef), names(coefnames[coefnames != missing.coef])))
    b[["mod1"]] <- setNames(c(coefs[o$predictor], coefs[paste(o$predictor, coefnames[coefnames != missing.coef], sep = ":")] + coefs[o$predictor]),
                   c(names(missing.coef), names(coefnames[coefnames != missing.coef])))
    
    variances <- diag(vc)
    v <- setNames(vector("list", 1+ length(o$moderator)),
                  c("intercept", paste0("mod", seq_along(o$moderator))))
    v[["intercept"]] <- setNames(c(variances[o$intercept], variances[coefnames[coefnames != missing.coef]] + variances[o$intercept] + 2*vc[o$intercept, coefnames[coefnames != missing.coef]]),
                  c(names(missing.coef), names(coefnames[coefnames != missing.coef])))
    v[["mod1"]] <- setNames(c(variances[o$predictor], variances[paste(o$predictor, coefnames[coefnames != missing.coef], sep = ":")] + variances[o$predictor] + 2*vc[o$predictor, paste(o$predictor, coefnames[coefnames != missing.coef], sep = ":")]),
                           c(names(missing.coef), names(coefnames[coefnames != missing.coef])))
    
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
  cz <- setNames(vector("list", length(moderator)), moderator)
  if (length(at.mod.level) > 0) {
    if (!is.list(at.mod.level)) {
      at.mod.level <- setNames(list(at.mod.level), moderator[seq_along(list(at.mod.level))])
    }
    if (length(at.mod.level) > length(moderator)) {
      plural <- length(moderator) > 1
      stop(paste(length(at.mod.level), "moderator levels were specified in at.mod.level, but there", ifelse(plural, "are", "is"),
           "only", length(moderator), ifelse(plural, "moderators", "moderator"), "specified."), call. = FALSE)
    }
    else if (length(intersect(names(at.mod.level), moderator)) < length(moderator)){
      new.at.mod.level <- at.mod.level[names(at.mod.level) %in% moderator]
      if (length(new.at.mod.level) < length(moderator)) {
        new.at.mod.level <- c(new.at.mod.level, setNames(lapply(seq_len(length(moderator) - length(new.at.mod.level)),
                                                       function(x) {if (x <= sum(!names(at.mod.level) %in% moderator)) at.mod.level[!names(at.mod.level) %in% moderator][[x]]
                                                       else NULL}),
                                                       moderator[!moderator%in%names(at.mod.level)]))
      }
      at.mod.level <- new.at.mod.level
    }
  
    if (length(mod.level.names) > 0) {
      if (!is.list(mod.level.names)) {
        mod.level.names <- setNames(list(mod.level.names), moderator[seq_along(list(mod.level.names))])
      }
      if (length(mod.level.names) > length(moderator)) {
        plural <- length(moderator) > 1
        stop(paste(length(mod.level.names), "moderator levels were specified in mod.level.names, but there", ifelse(plural, "are", "is"),
                   "only", length(moderator), ifelse(plural, "moderators", "moderator"), "specified."), call. = FALSE)
      }
      else if (length(intersect(names(mod.level.names), moderator)) < length(moderator)){
        new.mod.level.names <- mod.level.names[names(mod.level.names) %in% moderator]
        if (length(new.mod.level.names) < length(moderator)) {
          new.mod.level.names <- c(new.mod.level.names, setNames(lapply(seq_len(length(moderator) - length(new.mod.level.names)),
                                                                  function(x) {if (x <= sum(!names(mod.level.names) %in% moderator)) mod.level.names[!names(mod.level.names) %in% moderator][[x]]
                                                                    else NULL}),
                                                           moderator[!moderator %in% names(mod.level.names)]))
        }
        mod.level.names <- new.mod.level.names
      }
    } 
  } 
    
  for (i in moderator) {
    if (length(at.mod.level[[i]]) > 0) {
      if (is.numeric(at.mod.level[[i]])) {
        if (is.numeric(data[, i])) {
          cz[[i]] <- at.mod.level[[i]]
        }
        else if (is.logical(data[, i])) {
          if (!all(at.mod.level[[i]] %in% c(0, 1))) {
            warning(paste(i, "is logical, but the argument to at.mod.level contains values other than 0 or 1. Setting at.mod.level to c(FALSE, TRUE)."))
            cz[[i]] <- c(FALSE, TRUE)
          }
        }
        else if (is.factor(data[, i])) {
          cz[[i]] <- levels(data[, i])[at.mod.level[[i]]]
        }
        else {
          stop("The argument to at.mod.level is numeric but the moderator is not.", call. = FALSE)
        }
      }
      else if (is.character(at.mod.level[[i]])) {
        if (is.character(data[, i]) || is.factor(data[, i])) {
          cz[[i]] <- at.mod.level[[i]]
        }
        else {
          stop("The argument to at.mod.level is a character but the moderator is not a character or factor.", call. = FALSE)
        }
      }
      else if (is.logical(at.mod.level[[i]])) {
        if (is.logical(data[, i])) {
          cz[[i]] <- at.mod.level[[i]]
        }
        else {
          stop("The argument to at.mod.level is a logical but the moderator is not logical.", call. = FALSE)
        }
      }
      else {
        stop(paste("The argument to at.mod.level must be the same type as", i), call. = FALSE)
      }
      if (length(mod.level.names) > 0) {
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
        warning("An argument to sig.region was specified, but will be ignored. Using at.mod.level instead.", call. = FALSE)
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
      if (is.numeric(data[, i]) && length(unique(data[, i])) > 3) {
        mod.mean <- mean(data[, i])
        mod.sd <- sd(data[, i])
        cz[[i]] <- setNames(c(mod.mean - mod.sd, mod.mean, mod.mean + mod.sd),
                            c("Mean - SD", "Mean", "Mean + SD"))
      }
      else if (is.factor(data[, i])) {
        cz[[i]] <- levels(data[, i])
      }
      else {
        cz[[i]] <- sort(unique(data[, i]))
      }
    }
    if (length(names(cz[[i]])) == 0) {
      names(cz[[i]]) <- cz[[i]]
    }
  }
  return(cz)
}

between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
  if (!all(is.numeric(x))) stop("x must be a numeric vector.", call. = FALSE)
  if (length(range) != 2) stop("range must be of length 2.", call. = FALSE)
  if (any(is.na(range) | !is.numeric(range))) stop("range must contain numeric entries only.", call. = FALSE)
  range <- sort(range)
  
  if (any(is.na(x))) {
    if (length(na.action) != 1 || !is.atomic(na.action)) stop("na.action must be an atomic vector of length 1.", call. = FALSE)
  }
  if (inclusive) out <- ifelse(is.na(x), na.action, x >= range[1] & x <= range[2])
  else out <- ifelse(is.na(x), na.action, x > range[1] & x < range[2])
  
  return(out)
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
        stop("The input to fit must be an lm object or a glm object with family = \"gaussian\".", call. = FALSE)
      }
    }
    else stop("The input to fit must be an lm object or a glm object with family = \"gaussian\".", call. = FALSE)
  }
  if (length(all.vars(fit$terms)) < 3) {
    stop("The regression needs to contain at least 2 predictors and an outcome.", call. = FALSE)
  }
  if (all(classes[names(classes) %in% x.vars] != "numeric")) {
    stop("No variables in the model are numeric. At least the predictor must be numeric.", call. = FALSE)
  }
}

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  return(hcl(h = hues, l = 50, c = 100)[1:n])
}

isColor <- function(x) {
  tryCatch(is.matrix(col2rgb(x)), 
           error = function(e) FALSE)
}