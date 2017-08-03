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
  names(moderator) <- moderator
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
  
  link <- fit$family$link
  if (link %in% c("probit", "logit")) outcome <- paste0("P(", outcome, ")")
  
  vars <- setNames(c(intercept, predictor, moderator, unlist(ordered.ints)),
                   bnames)
  
  out <- list(outcome = outcome,
              intercept = intercept,
              predictor = predictor,
              moderator = moderator,
              interaction = interaction,
              vars = vars,
              link = link)
  attr(out, "cat.mod") <- setNames(classes[moderator] != "numeric",
                                   moderator)
  return(out)
}

fit.process2 <- function(o, fit) {
  data <- fit$model
  df <- fit$df.residual
  
  coefs <- coef(fit)
  vcov <- vcov(fit)
  
  new.vars <- o$vars
  has.int <- setNames(grepl(":", o$vars), o$vars)
  for (i in names(attr(o, "cat.mod"))[attr(o, "cat.mod") == TRUE]) {
    coefnames <- setNames(paste0(i, fit$xlevels[[i]]),
                          fit$xlevels[[i]])
    missing.coef <- coefnames[!coefnames %in% names(coefs)]
    
    o$moderator <- append(o$moderator, unname(coefnames[coefnames!=missing.coef]), which(o$moderator == i))[-which(o$moderator == i)]
    names(o$moderator)[names(o$moderator) == ""] <- i
    has.i <- setNames(grepl(i, o$vars), o$vars)
    for (j in o$vars) {
      if (has.i[j]) {
        if (has.int[j]) {
          split <- unlist(strsplit(j, ":", fixed = TRUE))
          new.v <- sapply(coefnames[coefnames!=missing.coef], 
                          function(x) {
                            s <- split; s[s==i] <- x
                            paste(s, collapse = ":")
                          })
          split.names <- unlist(strsplit(names(o$vars)[o$vars==j], "_", fixed = TRUE))
          new.names <- sapply(seq_along(coefnames[coefnames!=missing.coef]), 
                              function(x) {
                                s <- split.names; s[s==names(o$vars)[o$vars==i]] <- paste(s[s==names(o$vars)[o$vars==i]], x, sep = ".")
                                paste(s, collapse = "_")
                              })
        }
        else {
          new.v <- coefnames[coefnames!=missing.coef]
          new.names <- paste(names(new.vars)[new.vars==j],
                             seq_along(coefnames[coefnames!=missing.coef]),
                             sep = ".")
        }
        new.vars <- setNames(c(new.vars[new.vars!=j], new.v),
                             c(names(new.vars)[new.vars!=j], new.names))
      }
    }
  }
  o$vars <- new.vars
  
  in.vars <- setNames(names(coefs) %in% o$vars, names(coefs))
  names(coefs)[in.vars] <- sapply(names(coefs)[in.vars], 
                                  function(x) names(o$vars)[o$vars == x])
  dimnames(vcov) <- lapply(dimnames(vcov), function(x) names(coefs))
  indices.with.pred <- setNames(in.vars & grepl("pred", names(coefs), fixed = TRUE), names(coefs))
  indices.with.mod <- setNames(lapply(names(o$vars)[o$vars %in% o$moderator], 
                                      function(x) setNames(in.vars & grepl(x, names(coefs), fixed = TRUE), 
                                                           names(coefs))),
                               names(o$vars)[o$vars %in% o$moderator])
  indices.with.intercept <- setNames(in.vars & grepl("intercept", names(coefs), fixed = TRUE), names(coefs))
  
  out <- list(data = data,
              df = df,
              mod.exp = o$moderator,
              coefs = coefs,
              vcov = vcov, 
              vars = o$vars,
              indices.with = list(intercept = indices.with.intercept,
                                  pred = indices.with.pred,
                                  mod = indices.with.mod))
  
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

cz.prep <- function(o, B, at.mod.level = NULL, mod.level.names = NULL, data = NULL, sig.region = NULL) {
  moderator <- o$moderator
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
  
  cz <- setNames(vector("list", length(moderator)), moderator)
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
  
  if (FALSE) {
  cz.list <- setNames(vector("list", length(B$mod.exp)), B$mod.exp)
  for (i in moderator) {
    if (attr(o, "cat.mod")[i] == TRUE) {
      if (length(at.mod.level[[i]]) == 0) {
        at.mod.level[[i]] <- levels(data[, i])
      }
      if (length(at.mod.level[[i]]) > 0) {
        if (is.numeric(at.mod.level[[i]])) {
          at.mod.level[[i]] <- levels(data[, i])[at.mod.level[[i]]]
        }
        if (is.character(at.mod.level[[i]])) {
          cz.list[names(B$mod.exp) == i] <- lapply(names(cz.list[names(B$mod.exp) == i]),
                                                   function(x) as.numeric(paste0(i, at.mod.level[[i]]) == x))
        }
        else {
          stop(paste("The argument to at.mod.level must be the same type as", i), call. = FALSE)
        }
      }
      if (length(sig.region) > 0) {
        if (length(at.mod.level[[i]]) > 0) {
          warning("An argument to sig.region was specified, but will be ignored. Using at.mod.level instead.", call. = FALSE)
        }
        else {
          if (length(sig.region$bounds) == 2 && is.numeric(sig.region$bounds)) {
            at.mod.level[[i]] <- sig.region$bounds
            mod.level.names[[i]] <- c("Lower Bound", "Upper Bound")
          }
          else {
            warning("An argument to sig.region was specified, but it doesn't contain significance bounds, and will be ignored.", call. = FALSE)
          }
        }
      }
    }
    else {
      if (length(at.mod.level[[i]]) == 0) {
        if (length(unique(data[, i])) > 3) {
          mod.mean <- mean(data[, i])
          mod.sd <- sd(data[, i])
          cz.list[[i]] <- setNames(c(mod.mean - mod.sd, mod.mean, mod.mean + mod.sd),
                                   c("Mean - SD", "Mean", "Mean + SD"))
        }
        else {
          cz.list[[i]] <- sort(unique(data[, i]))
        }
      }
      else {
        if (is.numeric(at.mod.level[[i]])) {
          if (is.numeric(data[, i])) {
            cz.list[[i]] <- at.mod.level[[i]]
          }
          else if (is.logical(data[, i])) {
            if (!all(at.mod.level[[i]] %in% c(0, 1))) {
              warning(paste(i, "is logical, but the argument to at.mod.level contains values other than 0 or 1. Setting at.mod.level to c(FALSE, TRUE)."))
              cz.list[[i]] <- c(FALSE, TRUE)
            }
          }
          else if (is.factor(data[, i])) {
            cz.list[[i]] <- levels(data[, i])[at.mod.level[[i]]]
          }
          else {
            stop("The argument to at.mod.level is numeric but the moderator is not.", call. = FALSE)
          }
        }
        else if (is.logical(at.mod.level[[i]])) {
          if (is.logical(data[, i])) {
            cz.list[[i]] <- at.mod.level[[i]]
          }
          else {
            stop("The argument to at.mod.level is a logical but the moderator is not logical.", call. = FALSE)
          }
        }
        else {
          stop(paste("The argument to at.mod.level must be the same type as", i), call. = FALSE)
        }
        if (length(mod.level.names) > 0) {
          if (length(mod.level.names[[i]]) == length(cz.list[[i]])) {
            if (length(unique(mod.level.names[[i]])) != length(mod.level.names[[i]])) {
              warning("mod.level.names contains duplicate values. Ignoring mod.level.names.", call. = FALSE)
            }
            else names(cz.list[[i]]) <- mod.level.names[[i]]
          }
          else {
            warning("mod.level.names is not the same length as at.mod.level. Ignoring mod.level.names.", cal. = FALSE)
          }
        }
      }
      if (length(sig.region) > 0) {
        if (length(cz.list[[i]]) > 0) {
          warning("An argument to sig.region was specified, but will be ignored. Using at.mod.level instead.", call. = FALSE)
        }
        else {
          if (length(sig.region$bounds) == 2 && is.numeric(sig.region$bounds)) {
            cz.list[[i]] <- sig.region$bounds
            names(cz.list[[i]]) <- c("Lower Bound", "Upper Bound")
          }
          else {
            warning("An argument to sig.region was specified, but it doesn't contain significance bounds, and will be ignored.", call. = FALSE)
          }
        }
      }
    }
  }
  
  cz.list2 <- setNames(vector("list", length(moderator)), moderator)
  for (j in moderator) {
    if (attr(o, "cat.mod")[j] == TRUE) {
      cz.list.df <- as.data.frame(cz.list[B$mod.exp[names(B$mod.exp) == j]])
      cz.list2[[j]] <- as.character(unsplitfactor(cz.list.df, j, replace = TRUE, sep = "", 
                              dropped.level = levels(data[,j])[!paste0(j, levels(data[,j])) %in% B$mod.exp[names(B$mod.exp) == j]])[,])
      names(cz.list2[[j]]) <- cz.list2[[j]]
    }
    else cz.list2[[j]] <- cz.list[[j]]
  }

  return(list(cz = cz.list2,
              mod.exp.list = mod.exp.list))
  }
  mod.exp.list <- setNames(vector("list", sum(attr(o, "cat.mod"))), names(attr(o, "cat.mod"))[attr(o, "cat.mod")])
  for (i in names(mod.exp.list)) {
    mod.exp.list[[i]] <- setNames(data.frame(levels(data[, i]), model.matrix(~f, data = data.frame(f = factor(levels(data[, i]))))[,-1]),
                                 c(names(o$vars)[o$vars == i], names(B$vars)[sapply(B$vars, function(x) x %in% B$mod.exp[names(B$mod.exp) == i])]))
  }
  names(mod.exp.list) <- sapply(names(mod.exp.list), function(x) names(o$vars)[o$vars == x])
  
  return(list(cz = cz,
         mod.exp.list = mod.exp.list))
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
      if (fit$family$family != "gaussian" || fit$family$link != "identity") {
        #stop("The input to fit must be an lm object or a glm object with family = \"gaussian\".", call. = FALSE)
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

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[, nums] <- round(df[, nums], digits = digits)
  return(df)
}

unsplitfactor <- function(data, var.name, replace = TRUE, sep = "_", dropped.level = NULL, dropped.na = TRUE) {
  
  if (!is.data.frame(data)) stop("data must be a data.frame containing the variables to unsplit.", call = FALSE)
  if (!is.character(var.name)) stop("var.name must be a string containing the name of the variables to unsplit.", call. = FALSE)
  if (length(dropped.level) > 0 && length(var.name) > 1) {
    warning("dropped.level cannot be used with multiple var.names and will be ignored.", call. = FALSE, immediate. = TRUE)
    dropped.level <- NULL
  }
  
  if (!is.character(var.name)) stop("var.name must be a character vector containing the name of the variable to unsplit.", call. = FALSE)
  if (length(sep) > 1 || !is.character(sep)) stop("sep must be a character vector of length 1 containing the seperating character in the names of the split variables.", call. = FALSE)
  if (length(dropped.level) > 1 && !is.atomic(dropped.level)) {
    warning("dopped.level must be an atomic vector of length 1 containing the value of the dropped category of the split variable. It will be ignored.", call. = FALSE, immediate. = TRUE)
    dropped.level <- NULL
  }
  not.the.stem <- character(0)
  
  for (v in var.name) {
    dropped.level0 <- dropped.level
    var.to.combine <- data[, startsWith(names(data), paste0(v, sep)), drop = FALSE]
    if (length(var.to.combine) == 0) {
      not.the.stem <- c(not.the.stem, paste0(v, sep))
      next
    }
    
    if (!all(rowSums(apply(var.to.combine, 2, is.na)) %in% c(0, ncol(var.to.combine)))) {
      stop("The variables in data selected based on var.name and sep do not seem to form a split variable based on the <NA> pattern.", call. = FALSE)
    }
    NA.column <- character(0)
    
    if (!isTRUE(dropped.na)) {
      NA.column <- paste0(v, sep, ifelse(dropped.na == FALSE, "NA", dropped.na))
      if (NA.column %in% names(var.to.combine)) {
        var.to.combine[var.to.combine[, NA.column] == 1,] <- NA
        var.to.combine <- var.to.combine[, names(var.to.combine) != NA.column, drop = FALSE]
      }
      else {
        stop(paste("There is no variable called", word.list(NA.column, quotes = TRUE), "to generate the NA values."), call. = FALSE)
      }
    }
    var.sum <- rowSums(var.to.combine)
    if (isTRUE(all.equal(unique(var.sum), 1))) {
      #Already unsplit
    }
    else if (isTRUE(all.equal(sort(unique(var.sum)), c(0, 1)))) {
      #Missing category
      
      if (length(dropped.level) == 0) {
        k.levels0 <- sapply(names(var.to.combine), function(x) strsplit(x, paste0(v, sep))[[1]][2])
        
        if (suppressWarnings(all(!is.na(as.numeric(k.levels0))))) {
          dropped.level0 <- as.character(min(as.numeric(k.levels0)) - 1)
          dropped.name <- paste0(v, sep, dropped.level0)
        }
        else {
          message("The dropped category will be set to NA.")
          dropped.name <- dropped.level0 <- NA_character_
        }
        
      }
      else dropped.name <- paste0(v, sep, dropped.level)
      var.to.combine <- setNames(data.frame(1-var.sum, var.to.combine),
                                 c(dropped.name, names(var.to.combine)))
      
    }
    else {
      stop("The variables in data selected based on var.name and sep do not seem to form a split variable based on the row sums.", call. = FALSE)
    }
    
    k.levels <- sapply(names(var.to.combine), function(x) strsplit(x, paste0(v, sep))[[1]][2])
    
    k <- rep(NA_character_, nrow(data))
    for (i in seq_along(k.levels)) {
      k <- ifelse(var.to.combine[, i] == 1, k.levels[i], k)
    }
    
    k <- factor(k, levels = k.levels)
    
    
    if (replace) {
      where <- which(names(data) %in% c(names(var.to.combine), NA.column))
      
      data[, min(where)] <- k
      remove.cols <- where[where!=min(where)]
      if (length(remove.cols) > 0) data <- data[, -remove.cols, drop = FALSE]
      names(data)[min(where)] <- v
    }
    else {
      data <- cbind(data, setNames(data.frame(k), v))
    }
  }
  
  if (length(not.the.stem) > 0) warning(paste0(word.list(not.the.stem, is.are = TRUE, quotes = TRUE), " not the stem of any variables in data and will be ignored. Ensure var.name and sep are correct."), call. = FALSE)
  
  return(data)
}

splitfactor <- function(data, var.name, replace = TRUE, sep = "_", drop.level = NULL, drop.first = c(TRUE, FALSE, "if2"), drop.singleton = FALSE, drop.na = TRUE, check = TRUE) {
  #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
  #Retains all categories unless only 2 levels, in which case only the second level is retained.
  #If variable only has one level, will delete.
  #var.name= the name of the variable to split when data is specified
  #data=data set to be changed
  
  if (is.data.frame(data)) {
    if (check) {
      factor.names <- names(data)[sapply(data, function(x) is.factor(x) || is.character(x))]
      if (missing(var.name)) {
        var.name <- factor.names
      }
      else if (is.character(var.name)) {
        if (any(var.name %in% factor.names)) {
          if (any(!var.name %in% factor.names)) {
            not.in.factor.names <- var.name[!var.name %in% factor.names]
            warning(paste(word.list(not.in.factor.names, "and", is.are = TRUE), 
                          "not the name(s) of factor variable(s) in data and will not be split."), 
                    call. = FALSE)
          }
          var.name <- var.name[var.name %in% factor.names]
        }
        else {
          stop("No names in var.name are names of factor variables in data.", call. = FALSE)
        }
      }
      else {
        stop("var.name must be a character vector of the name(s) of factor variable(s) in data.", call. = FALSE)
      }
      if (length(factor.names) == 0) {
        stop("There are no factor variables to split in data.", call. = FALSE)
      }
    }
    else {
      if (missing(var.name) || !is.character(var.name)) {
        stop("var.name must be a character vector of the names of variables in data.", call. = FALSE)
      }
      else {
        if (any(var.name %in% names(data))) {
          if (any(!var.name %in% names(data))) {
            not.in.data.names <- var.name[!var.name %in% names(data)]
            warning(paste(word.list(not.in.data.names, "and", is.are = TRUE), 
                          "not the name(s) of variable(s) in data and will not be split."), 
                    call. = FALSE)
          }
          var.name <- var.name[var.name %in% names(data)]
        }
        else {
          stop("No names in var.name are names of variables in data.", call. = FALSE)
        }
      }
    }
    
  }
  else if (is.atomic(data)) {
    dep <- deparse(substitute(data))
    data <- data.frame(data)
    if (missing(var.name)) {
      names(data) <- dep
    }
    else if (is.vector(var.name) && (is.atomic(var.name) || is.factor(var.name))) {
      if (length(var.name) == 0) {
        names(data) <- dep
      }
      else if (length(var.name) == 1) {
        names(data) <- var.name
      }
      else {
        warning("Only using the first item of var.name.", call. = FALSE)
        names(data) <- var.name[1]
      }
    }
    else {
      stop("var.name must be an atomic or factor vector of length 1 with the stem of the new variable.", call. = FALSE)
    }
    var.name <- names(data)
  }
  else {
    stop("data must a be a data.frame or an atomic vector.", call. = FALSE)
  }
  
  if (length(drop.level) > 0 && length(var.name) > 1) {
    warning("drop.level cannot be used with multiple entries to var.name. Ignoring drop.level.", call. = FALSE)
    drop.level <- NULL
  }
  drop.na <- setNames(rep(drop.na, length(var.name)), var.name)
  for (v in var.name) {
    drop <- character(0)
    x <- data[, names(data) == v] <- factor(data[, names(data) == v], exclude = NULL)
    
    skip <- FALSE
    if (nlevels(x) > 1) {
      k <- model.matrix(as.formula(paste0("~", v, "- 1")), data = data)
      
      if (any(is.na(levels(x)))) {
        
        if (drop.na[v]) {
          k[k[,is.na(levels(x))] == 1,] <- NA
          #k <- k[, !is.na(levels(x)), drop = FALSE]
        }
      }
      else drop.na[v] <- FALSE
      
    }
    else {
      if (drop.singleton) {
        data <- data[, names(data)!=v, drop = FALSE]
        skip <- TRUE
      }
      else {
        k <- matrix(1, ncol = 1, nrow = length(x))
        colnames(k) <- paste0(v, levels(x)[1])
      }
    }
    
    if (!skip) {
      colnames(k) <- paste(v, sapply(strsplit(colnames(k), v, fixed = TRUE), function(n) paste(n, collapse = "")), sep = sep)
      
      if (length(drop.level) > 0) {
        if (is.character(drop.level) && length(drop.level) == 1 && drop.level %in% levels(x)) {
          drop <- drop.level
        }
        else {
          stop(paste("drop must be the name of a level of", v, "which is to be dropped."), call. = FALSE)
        }
      }
      else {
        if ((ncol(k) == 2 && drop.first %in% c("if2", TRUE)) ||
            (ncol(k) > 2 && drop.first == TRUE)) {
          drop <- levels(x)[1]
        }
      }
      
      dropl <- rep(FALSE, ncol(k))
      if (length(drop) > 0) {
        dropl[!is.na(levels(x)) & levels(x) %in% drop] <- TRUE
      }
      if (drop.na[v]) dropl[is.na(levels(x))] <- TRUE
      k <- k[, !dropl, drop = FALSE]
      
      if (ncol(data) == 1) {
        data <- data.frame(k)
      }
      else if (replace) {
        if (match(v, names(data)) == 1){
          data <- data.frame(k, data[, names(data)!=v, drop = FALSE])
        }
        else if (match(v, names(data)) == ncol(data)) {
          data <- data.frame(data[, names(data)!=v, drop = FALSE], k)
        }
        else {
          where <- match(v, names(data))
          data <- data.frame(data[, 1:(where-1), drop = FALSE], k, data[, (where+1):ncol(data), drop = FALSE])
        }
      }
      else {
        data <- data.frame(data, k)
      }
      
    }
    
  }
  
  return(data)
}
