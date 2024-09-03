AObootBetween <- function(
    var.between,
    var.dv,
    var.id,
    levels.b1,
    levels.b2 = NULL,
    eff.si = c("pes", "ges"),
    data,
    silence = FALSE,
    n.sim = 1000,
    alpha = .05,
    seed = 1234,
    n.round = 2){

  set.seed(seed)

  if (is.null(levels.b2)){
    no.test = 0
    no.aov = 0

    ## Prepare dataset
    names(data)[names(data) == var.id] <- "id"

    names(data)[names(data) == var.dv] <- "value"

    names(data)[names(data) == var.between] <- "factor1"

    all.var <- c("id", "factor1", "value")

    boots.data <- data[, colnames(data) %in% all.var]

    ## Prepare arrays
    names <- c(NA)

    k = 1

    for(i in 1:length(levels.b1)){
      for(j in 2:length(levels.b1)){
        if(i < j){
          names[k] <- paste0(levels.b1[i], " - ", levels.b1[j])
          k = k + 1
        }
      }
    }

    array.F <- array(NA, dim = c(1, 1, n.sim))
    dimnames(array.F) <- list("group",
                              "F",
                              1:n.sim)

    array.es <- array(NA, dim = c(1, 1, n.sim))
    dimnames(array.es) <- list("group",
                               "eta",
                               1:n.sim)

    array.em <- array(NA, dim = c(length(levels.b1), 2, n.sim))
    dimnames(array.em) <- list(levels.b1,
                               c("est", "SE"),
                               1:n.sim)

    array.ph <- array(NA, dim = c(length(names), 3, n.sim))
    dimnames(array.ph) <- list(names,
                               c("est", "SE", "d"),
                               1:n.sim)

    ## original aov to compare
    boots.data$factor1 <- factor(boots.data$factor1, levels = levels.b1)

    suppressMessages(orig.aov <- aov_car(value ~ factor1 + Error(id),
                                         data = boots.data,
                                         fun_aggregate = mean,
                                         anova_table = list(es = eff.si)))

      for(i in 1:n.sim){ # bootstrapping loop

        if(!silence){

          message("Progress:", i, "/", n.sim)

        }

      ## Sample
      order <- sample(1:nrow(boots.data),
                      replace = TRUE)
      data.BS <- boots.data[order, ]

      data.BS$factor1 <- factor(data.BS$factor1, levels = levels.b1)


      ## Control for all cells
      boots.con <- table(data.BS$factor1)

      if (any(boots.con == 0) == TRUE) {
        next
      }

      ## Test
      suppressMessages(boots.aov <- aov_car(value ~ factor1 + Error(id),
                                            data = data.BS,
                                            fun_aggregate = mean,
                                            anova_table = list(es = eff.si)))

      no.aov = no.aov + 1

      ## Array aov
      df.F <- boots.aov$anova_table$F

      array.F[, , i] <- df.F

      df.es <- eval(parse(text=paste0("boots.aov$anova_table$", eff.si)))

      array.es[, , i] <- df.es

      ## Post hoc
      boots.em <- suppressMessages(emmeans(boots.aov, specs = "factor1"))

      boot.SE <- as.data.frame(summary(boots.em))

      # check for unique SEs to avoid NaNs
      vec.SE <- round(boot.SE$SE, 6)

      if(length(unique(vec.SE)) < length(boot.SE$SE)) {
        next
      }

      df.em <- c(summary(boots.em)$emmean,
                 summary(boots.em)$SE)

      array.em[, , i] <- df.em

      ## Array post hoc
      no.test = no.test + 1

      boots.pairs <- pairs(boots.em,
                           adjust = "bonferroni")

      d.ph <- as.numeric(summary(
        eff_size(boots.em,
                 sigma = summary.lm(boots.aov$lm)$sigma,
                 edf = boots.aov$lm$df.residual)
      )$effect.size)

      df.ph <- c(summary(boots.pairs)$estimate,
                 summary(boots.pairs)$SE,
                 d.ph)

      array.ph[, , i] <- df.ph

    } # bootstrapping loop

    out.F <- list(
      length(which(array.F[1, , ] > orig.aov$anova_table$F[1]))/n.sim
    )

    orig.aov$anova_table$`Pr(>F)` = as.numeric(out.F)

    out.es <- list(mean = apply(array.es, c(1, 2),
                                function(x){
                                  mean(x, na.rm = TRUE)}),
                   LL = apply(array.es, c(1, 2),
                              function(x){
                                quantile(x, probs = alpha / 2, na.rm = TRUE)}),
                   UL = apply(array.es, c(1, 2),
                              function(x){
                                quantile(x, probs = 1 - (alpha / 2),
                                         na.rm = TRUE)}))

    orig.aov$anova_table$`mean effect size` <- out.es$mean
    orig.aov$anova_table$`LL effect size` <- out.es$LL
    orig.aov$anova_table$`UL effect size` <- out.es$UL
    orig.aov$anova_table$`no tests` <- no.aov

    out.em <- list(mean = apply(array.em, c(1, 2),
                                function(x){
                                  mean(x, na.rm = TRUE)}),
                   LL = apply(array.em, c(1, 2),
                              function(x){
                                quantile(x, probs = alpha / 2, na.rm = TRUE)}),
                   UL = apply(array.em, c(1, 2),
                              function(x){
                                quantile(x, probs = 1 - (alpha / 2),
                                         na.rm = TRUE)}))

    dat.em <- data.frame("group" = levels.b1)
    dat.em$`est mean` <- round(out.em$mean[1:length(levels.b1)], n.round)
    dat.em$`est LL` <- round(out.em$LL[1:length(levels.b1)], n.round)
    dat.em$`est UL` <- round(out.em$UL[1:length(levels.b1)], n.round)
    dat.em$`SE mean` <- round(out.em$mean[(length(levels.b1) + 1):
                                            (2 * length(levels.b1))], n.round)
    dat.em$`SE LL` <- round(out.em$LL[(length(levels.b1) + 1):
                                        (2 * length(levels.b1))], n.round)
    dat.em$`SE UL` <- round(out.em$UL[(length(levels.b1) + 1):
                                        (2 * length(levels.b1))], n.round)

    out.ph <- list(mean = apply(array.ph, c(1, 2),
                                function(x){
                                  mean(x, na.rm = TRUE)}),
                   LL = apply(array.ph, c(1, 2),
                              function(x){
                                quantile(x, probs = alpha / 2, na.rm = TRUE)}),
                   UL = apply(array.ph, c(1, 2),
                              function(x){
                                quantile(x, probs = 1 - (alpha / 2),
                                         na.rm = TRUE)}))

    dat.ph <- data.frame("comparison" = names)
    dat.ph$`d mean` <- round(out.ph$mean[(2 * length(names) + 1):
                                           (3 * length(names))], n.round)
    dat.ph$`d LL` <- round(out.ph$LL[(2 * length(names) + 1):
                                       (3 * length(names))], n.round)
    dat.ph$`d UL` <- round(out.ph$UL[(2 * length(names) + 1):
                                       (3 * length(names))], n.round)
    dat.ph$`est mean` <- round(out.ph$mean[1:length(names)], n.round)
    dat.ph$`est LL` <- round(out.ph$LL[1:length(names)], n.round)
    dat.ph$`est UL` <- round(out.ph$UL[1:length(names)], n.round)
    dat.ph$`SE mean` <- round(out.ph$mean[(length(names) + 1):
                                            (2 * length(names))], n.round)
    dat.ph$`SE LL` <- round(out.ph$LL[(length(names) + 1):
                                        (2 * length(names))], n.round)
    dat.ph$`SE UL` <- round(out.ph$UL[(length(names) + 1):
                                        (2 * length(names))], n.round)

    output <- list(type.aov = "One-way between ANOVA",
                   factor = levels.b1,
                   anova = round(orig.aov$anova_table, n.round),
                   em = dat.em,
                   no.test = no.test,
                   ph = dat.ph)

    class(output) <- "AOboot_one"

    return(output)

  } else if (!is.null(levels.b2)) {
    no.test1 = 0
    no.test2 = 0
    no.test3 = 0
    no.test4 = 0
    no.aov = 0

    ## Prepare data
    names(data)[names(data) == var.id] <- "id"

    names(data)[names(data) == var.dv] <- "value"

    for (i in 1:length(var.between)) {
      names(data)[names(data) == var.between[i]] <- paste0("factor", i)
    }

    all.var <- c("id", "factor1", "factor2", "value")

    boots.data <- data[, colnames(data) %in% all.var]

    ## Prepare arrays
    names.b1 <- c(NA)

    k = 1

    for(i in 1:length(levels.b1)){
      for(j in 2:length(levels.b1)){
        if(i < j){
          names.b1[k] <- paste0(levels.b1[i], " - ", levels.b1[j])
          k = k + 1
        }
      }
    }

    names.b2 <- c(NA)

    k = 1

    for(i in 1:length(levels.b2)){
      for(j in 2:length(levels.b2)){
        if(i < j){
          names.b2[k] <- paste0(levels.b2[i], " - ", levels.b2[j])
          k = k + 1
        }
      }
    }

    array.F <- array(NA, dim = c(3, 1, n.sim))
    dimnames(array.F) <- list(c("factor 1", "factor 2", "interaction"),
                              "F",
                              1:n.sim)

    array.es <- array(NA, dim = c(3, 1, n.sim))
    dimnames(array.es) <- list(c("factor 1", "factor 2", "interaction"),
                               "eta",
                               1:n.sim)

    array.em1 <- array(NA, dim = c(length(levels.b1), 2, n.sim))
    dimnames(array.em1) <- list(levels.b1,
                                c("est", "SE"),
                                1:n.sim)

    array.ph1 <- array(NA, dim = c(length(names.b1), 3, n.sim))
    dimnames(array.ph1) <- list(names.b1,
                                c("est", "SE", "d"),
                                1:n.sim)

    array.em2 <- array(NA, dim = c(length(levels.b2), 2, n.sim))
    dimnames(array.em2) <- list(levels.b2,
                                c("est", "SE"),
                                1:n.sim)

    array.ph2 <- array(NA, dim = c(length(names.b2), 3, n.sim))
    dimnames(array.ph2) <- list(names.b2,
                                c("est", "SE", "d"),
                                1:n.sim)

    array.em3 <- array(NA, dim = c(length(levels.b1), length(levels.b2), 2,
                                  n.sim))
    dimnames(array.em3) <- list(levels.b1,
                                levels.b2,
                                c("est", "SE"),
                                1:n.sim)

    array.ph3 <- array(NA, dim = c(length(names.b1), length(levels.b2), 3,
                                  n.sim))
    dimnames(array.ph3) <- list(names.b1,
                                levels.b2,
                                c("est", "SE", "d"),
                                1:n.sim)

    array.em4 <- array(NA, dim = c(length(levels.b2), length(levels.b1), 2,
                                  n.sim))
    dimnames(array.em4) <- list(levels.b2,
                                levels.b1,
                                c("est", "SE"),
                                1:n.sim)

    array.ph4 <- array(NA, dim = c(length(names.b2), length(levels.b1), 3,
                                  n.sim))
    dimnames(array.ph4) <- list(names.b2,
                                levels.b1,
                                c("est", "SE", "d"),
                                1:n.sim)

    array.em.list <- list(array.em1, array.em2, array.em3, array.em4)
    array.ph.list <- list(array.ph1, array.ph2, array.ph3, array.ph4)
    no.list <- list(no.test1, no.test2, no.test3, no.test4)

    ## original aov to compare
    boots.data$factor1 <- factor(boots.data$factor1, levels = levels.b1)
    boots.data$factor2 <- factor(boots.data$factor2, levels = levels.b2)

    suppressMessages(orig.aov <- aov_car(value ~ factor1 * factor2 + Error(id),
                                         data = boots.data,
                                         fun_aggregate = mean,
                                         anova_table = list(es = eff.si)))


      for(i in 1:n.sim){ # bootstrapping loop

        if(!silence){

          message("Progress:", i, "/", n.sim)

        }

      ## Sample
      order <- sample(1:nrow(boots.data),
                      replace = TRUE)
      data.BS <- boots.data[order, ]

      data.BS$factor1 <- factor(data.BS$factor1, levels = levels.b1)
      data.BS$factor2 <- factor(data.BS$factor2, levels = levels.b2)

      ## Control for all cells
      boots.con <- table(data.BS$factor1, data.BS$factor2)

      if (any(boots.con == 0) == TRUE) {
        next
      }

      ## Test
      suppressMessages(boots.aov <- aov_car(
        value ~ factor1 * factor2 + Error(id),
        data = data.BS,
        fun_aggregate = mean,
        anova_table = list(es = eff.si)))

      no.aov = no.aov + 1

      ## Array aov
      df.F <- boots.aov$anova_table$F

      array.F[, , i] <- df.F

      df.es <- eval(parse(text=paste0("boots.aov$anova_table$", eff.si)))

      array.es[, , i] <- df.es

      ## Post hoc
      m = 1

      for (j in 1:nrow(boots.aov$anova_table)) { # post hoc loop
        ronam <- rownames(boots.aov$anova_table)[j]

        m = j

        if(length(strsplit(ronam, ":")[[1]]) == 1) { # effects if
          boots.em <- suppressMessages(emmeans(boots.aov, specs = ronam))

          boot.SE <- as.data.frame(summary(boots.em))

          # check for unique SEs to avoid NaNs
          br = 0

          vec.SE <- round(boot.SE$SE, 6)

          if(length(unique(vec.SE)) < length(boot.SE$SE)) {
            br = 1
          }

          if(br == 0) {
            no.list[m] = no.list[m][[1]] + 1

            df.em <- c(summary(boots.em)$emmean,
                       summary(boots.em)$SE)

            array.em.list[[m]][, , i] <- df.em


            boots.pairs <- pairs(boots.em, adjust = "bonferroni")

            d.ph <- as.numeric(summary(
              eff_size(boots.em,
                       sigma = summary.lm(boots.aov$lm)$sigma,
                       edf = boots.aov$lm$df.residual)
            )$effect.size)

            df.ph <- c(summary(boots.pairs)$estimate,
                       summary(boots.pairs)$SE,
                       d.ph)

            array.ph.list[[m]][, , i] <- df.ph
          }

        } else if (length(strsplit(ronam, ":")[[1]]) == 2) { # effects if
          boots.em <- suppressMessages(
            emmeans(boots.aov, specs = strsplit(ronam, ":")[[1]][1],
                    by =  strsplit(ronam, ":")[[1]][2]))

          boot.SE <- as.data.frame(summary(boots.em))

          # check for unique SEs to avoid NaNs
          for (k in 1:length(levels(summary(boots.em)$factor2))) {
            br = 0

            vec.SE <- round(boot.SE$SE[which(boot.SE$factor2 ==
                                               levels(summary(
                                                 boots.em)$factor2)[k])],
                            6)


            if(length(unique(vec.SE)) < length(levels(
              summary(boots.em)$factor1))) {
              br = 1
              break
            }
          }

          if (br == 0) {
            no.list[m] = no.list[m][[1]] + 1

            df.em <- c(summary(boots.em)$emmean,
                       summary(boots.em)$SE)

            array.em.list[[m]][, , ,i] <- df.em


            boots.pairs <- pairs(boots.em, adjust = "bonferroni")

            d.ph <- as.numeric(summary(
              eff_size(boots.em,
                       sigma = summary.lm(boots.aov$lm)$sigma,
                       edf = boots.aov$lm$df.residual)
            )$effect.size)

            df.ph <- c(summary(boots.pairs)$estimate,
                       summary(boots.pairs)$SE,
                       d.ph)

            array.ph.list[[m]][, , , i] <- df.ph
          }


          m = m + 1


          boots.em <- suppressMessages(
            emmeans(boots.aov, specs = strsplit(ronam, ":")[[1]][2],
                    by =  strsplit(ronam, ":")[[1]][1]))

          boot.SE <- as.data.frame(summary(boots.em))

          # check for unique SEs to avoid NaNs
          for (k in 1:length(levels(summary(boots.em)$factor1))) {
            br = 0

            vec.SE <- round(boot.SE$SE[which(boot.SE$factor1 ==
                                               levels(summary(
                                                 boots.em)$factor1)[k])],
                            6)

            if(length(unique(vec.SE)) < length(levels(
              summary(boots.em)$factor2))) {
              br = 1
              break
            }
          }

          if (br == 0) {
            no.list[m] = no.list[m][[1]] + 1

            df.em <- c(summary(boots.em)$emmean,
                       summary(boots.em)$SE)

            array.em.list[[m]][, , , i] <- df.em


            boots.pairs <- pairs(boots.em, adjust = "bonferroni")

            d.ph <- as.numeric(summary(
              eff_size(boots.em, sigma = summary.lm(boots.aov$lm)$sigma,
                       edf = boots.aov$lm$df.residual))$effect.size)

            df.ph <- c(summary(boots.pairs)$estimate,
                       summary(boots.pairs)$SE,
                       d.ph)

            array.ph.list[[m]][, , , i] <- df.ph
          }

        } #effects if

      } #post hoc loop

    } # bootstrapping loop

    no.test1 <- no.list[1][[1]]
    no.test2 <- no.list[2][[1]]
    no.test3 <- no.list[3][[1]]
    no.test4 <- no.list[4][[1]]

    out.F <- list(
      length(which(array.F[1, , ] > orig.aov$anova_table$F[1]))/n.sim,
      length(which(array.F[2, , ] > orig.aov$anova_table$F[2]))/n.sim,
      length(which(array.F[3, , ] > orig.aov$anova_table$F[3]))/n.sim
    )

    orig.aov$anova_table$`Pr(>F)` = as.numeric(out.F)

    out.es <- list(mean = apply(array.es, c(1, 2),
                                function(x){
                                  mean(x, na.rm = TRUE)}),
                   LL = apply(array.es, c(1, 2),
                              function(x){
                                quantile(x, probs = alpha / 2, na.rm = TRUE)}),
                   UL = apply(array.es, c(1, 2),
                              function(x){
                                quantile(x, probs = 1 - (alpha / 2),
                                         na.rm = TRUE)}))

    orig.aov$anova_table$`mean effect size` <- out.es$mean
    orig.aov$anova_table$`LL effect size` <- out.es$LL
    orig.aov$anova_table$`UL effect size` <- out.es$UL
    orig.aov$anova_table$`no tests` <- no.aov

    array.em1 <- array.em.list[[1]]
    array.em2 <- array.em.list[[2]]
    array.em3 <- array.em.list[[3]]
    array.em4 <- array.em.list[[4]]

    array.ph1 <- array.ph.list[[1]]
    array.ph2 <- array.ph.list[[2]]
    array.ph3 <- array.ph.list[[3]]
    array.ph4 <- array.ph.list[[4]]

    out.em1 <- list(mean = apply(array.em1, c(1, 2),
                                 function(x){mean(x, na.rm = TRUE)}),
                    LL = apply(array.em1, c(1, 2),
                               function(x){quantile(x, probs = alpha / 2,
                                                    na.rm = TRUE)}),
                    UL = apply(array.em1, c(1, 2),
                               function(x){quantile(x,
                                                    probs = 1 - (alpha / 2),
                                                    na.rm = TRUE)}))

    out.em2 <- list(mean = apply(array.em2, c(1, 2),
                                 function(x){mean(x, na.rm = TRUE)}),
                    LL = apply(array.em2, c(1, 2),
                               function(x){quantile(x,
                                                    probs = alpha / 2,
                                                    na.rm = TRUE)}),
                    UL = apply(array.em2, c(1, 2),
                               function(x){quantile(x,
                                                    probs = 1 - (alpha / 2),
                                                    na.rm = TRUE)}))

    out.em3 <- list(mean = apply(array.em3, c(1, 2, 3),
                                 function(x){mean(x, na.rm = TRUE)}),
                    LL = apply(array.em3, c(1, 2, 3),
                               function(x){quantile(x, probs = alpha / 2,
                                                    na.rm = TRUE)}),
                    UL = apply(array.em3, c(1, 2, 3),
                               function(x){quantile(x,
                                                    probs = 1 - (alpha / 2),
                                                    na.rm = TRUE)}))

    out.em4 <- list(mean = apply(array.em4, c(1, 2, 3),
                                 function(x){mean(x, na.rm = TRUE)}),
                    LL = apply(array.em4, c(1, 2, 3),
                               function(x){quantile(x, probs = alpha / 2,
                                                    na.rm = TRUE)}),
                    UL = apply(array.em4, c(1, 2, 3),
                               function(x){quantile(x,
                                                    probs = 1 - (alpha / 2),
                                                    na.rm = TRUE)}))

    out.ph1 <- list(mean = apply(array.ph1, c(1, 2),
                                 function(x){mean(x, na.rm = TRUE)}),
                    LL = apply(array.ph1, c(1, 2),
                               function(x){quantile(x, probs = alpha / 2,
                                                    na.rm = TRUE)}),
                    UL = apply(array.ph1, c(1, 2),
                               function(x){quantile(x,
                                                    probs = 1 - (alpha / 2),
                                                    na.rm = TRUE)}))

    out.ph2 <- list(mean = apply(array.ph2, c(1, 2),
                                 function(x){mean(x, na.rm = TRUE)}),
                    LL = apply(array.ph2, c(1, 2),
                               function(x){quantile(x, probs = alpha / 2,
                                                    na.rm = TRUE)}),
                    UL = apply(array.ph2, c(1, 2),
                               function(x){quantile(x,
                                                    probs = 1 - (alpha / 2),
                                                    na.rm = TRUE)}))

    out.ph3 <- list(mean = apply(array.ph3, c(1, 2, 3),
                                 function(x){mean(x, na.rm = TRUE)}),
                    LL = apply(array.ph3, c(1, 2, 3),
                               function(x){quantile(x, probs = alpha / 2,
                                                    na.rm = TRUE)}),
                    UL = apply(array.ph3, c(1, 2, 3),
                               function(x){quantile(x,
                                                    probs = 1 - (alpha / 2),
                                                    na.rm = TRUE)}))

    out.ph4 <- list(mean = apply(array.ph4, c(1, 2, 3),
                                 function(x){mean(x, na.rm = TRUE)}),
                    LL = apply(array.ph4, c(1, 2, 3),
                               function(x){quantile(x, probs = alpha / 2,
                                                    na.rm = TRUE)}),
                    UL = apply(array.ph4, c(1, 2, 3),
                               function(x){quantile(x,
                                                    probs = 1 - (alpha / 2),
                                                    na.rm = TRUE)}))

    dat.em1 <- data.frame("group" = levels.b1)
    dat.em1$`est mean` <- round(out.em1$mean[1:length(levels.b1)], n.round)
    dat.em1$`est LL` <- round(out.em1$LL[1:length(levels.b1)], n.round)
    dat.em1$`est UL` <- round(out.em1$UL[1:length(levels.b1)], n.round)
    dat.em1$`SE mean` <- round(out.em1$mean[(length(levels.b1) + 1):
                                              (2 * length(levels.b1))], n.round)
    dat.em1$`SE LL` <- round(out.em1$LL[(length(levels.b1) + 1):
                                          (2 * length(levels.b1))], n.round)
    dat.em1$`SE UL` <- round(out.em1$UL[(length(levels.b1) + 1):
                                          (2 * length(levels.b1))], n.round)

    dat.em2 <- data.frame("group" = levels.b2)
    dat.em2$`est mean` <- round(out.em2$mean[1:length(levels.b2)], n.round)
    dat.em2$`est LL` <- round(out.em2$LL[1:length(levels.b2)], n.round)
    dat.em2$`est UL` <- round(out.em2$UL[1:length(levels.b2)], n.round)
    dat.em2$`SE mean` <- round(out.em2$mean[(length(levels.b2) + 1):
                                              (2 * length(levels.b2))], n.round)
    dat.em2$`SE LL` <- round(out.em2$LL[(length(levels.b2) + 1):
                                          (2 * length(levels.b2))], n.round)
    dat.em2$`SE UL` <- round(out.em2$UL[(length(levels.b2) + 1):
                                          (2 * length(levels.b2))], n.round)

    dat.em3 <- expand.grid("group factor 1" = levels.b1,
                           "group factor 2" = levels.b2)
    dat.em3$`est mean` <- round(out.em3$mean[1:(length(levels.b2) *
                                                  length(levels.b1))], n.round)
    dat.em3$`est LL` <- round(out.em3$LL[1:(length(levels.b2) *
                                              length(levels.b1))], n.round)
    dat.em3$`est UL` <- round(out.em3$UL[1:(length(levels.b2) *
                                              length(levels.b1))], n.round)
    dat.em3$`SE mean` <- round(out.em3$mean[(length(levels.b2) *
                                               length(levels.b1) + 1):
                                              ((2 * length(levels.b2)) *
                                                 length(levels.b1))], n.round)
    dat.em3$`SE LL` <- round(out.em3$LL[(length(levels.b2) * length(levels.b1)
                                         + 1):((2 * length(levels.b2)) *
                                                 length(levels.b1))], n.round)
    dat.em3$`SE UL` <- round(out.em3$UL[(length(levels.b2) * length(levels.b1)
                                         + 1):((2 * length(levels.b2)) *
                                                 length(levels.b1))], n.round)

    dat.em3 <- dat.em3[, c(2, 1, 3:8)]
    dat.em3 <- dat.em3[order(dat.em3$`group factor 2`), ]

    dat.em4 <- expand.grid("group factor 2" = levels.b2,
                           "group factor 1" = levels.b1)
    dat.em4$`est mean` <- round(out.em4$mean[1:(length(levels.b1) *
                                                  length(levels.b2))], n.round)
    dat.em4$`est LL` <- round(out.em4$LL[1:(length(levels.b1) *
                                              length(levels.b2))], n.round)
    dat.em4$`est UL` <- round(out.em4$UL[1:(length(levels.b1) *
                                              length(levels.b2))], n.round)
    dat.em4$`SE mean` <- round(out.em4$mean[(length(levels.b1) *
                                               length(levels.b2) + 1):
                                              ((2 * length(levels.b1)) *
                                                 length(levels.b2))], n.round)
    dat.em4$`SE LL` <- round(out.em4$LL[(length(levels.b1) * length(levels.b2)
                                         + 1):((2 * length(levels.b1)) *
                                                 length(levels.b2))], n.round)
    dat.em4$`SE UL` <- round(out.em4$UL[(length(levels.b1) * length(levels.b2)
                                         + 1):((2 * length(levels.b1)) *
                                                 length(levels.b2))], n.round)

    dat.em4 <- dat.em4[, c(2, 1, 3:8)]
    dat.em4 <- dat.em4[order(dat.em4$`group factor 1`), ]

    dat.ph1 <- data.frame("comparison" = names.b1)
    dat.ph1$`d mean` <- round(out.ph1$mean[(2 * length(names.b1) + 1):
                                             (3 * length(names.b1))], n.round)
    dat.ph1$`d LL` <- round(out.ph1$LL[(2 * length(names.b1) + 1):
                                         (3 * length(names.b1))], n.round)
    dat.ph1$`d UL` <- round(out.ph1$UL[(2 * length(names.b1) + 1):
                                         (3 * length(names.b1))], n.round)
    dat.ph1$`est mean` <- round(out.ph1$mean[1:length(names.b1)], n.round)
    dat.ph1$`est LL` <- round(out.ph1$LL[1:length(names.b1)], n.round)
    dat.ph1$`est UL` <- round(out.ph1$UL[1:length(names.b1)], n.round)
    dat.ph1$`SE mean` <- round(out.ph1$mean[(length(names.b1) + 1):
                                              (2 * length(names.b1))], n.round)
    dat.ph1$`SE LL` <- round(out.ph1$LL[(length(names.b1) + 1):
                                          (2 * length(names.b1))], n.round)
    dat.ph1$`SE UL` <- round(out.ph1$UL[(length(names.b1) + 1):
                                          (2 * length(names.b1))], n.round)

    dat.ph2 <- data.frame("comparison" = names.b2)
    dat.ph2$`d mean` <- round(out.ph2$mean[(2 * length(names.b2) + 1):
                                             (3 * length(names.b2))], n.round)
    dat.ph2$`d LL` <- round(out.ph2$LL[(2 * length(names.b2) + 1):
                                         (3 * length(names.b2))], n.round)
    dat.ph2$`d UL` <- round(out.ph2$UL[(2 * length(names.b2) + 1):
                                         (3 * length(names.b2))], n.round)
    dat.ph2$`est mean` <- round(out.ph2$mean[1:length(names.b2)], n.round)
    dat.ph2$`est LL` <- round(out.ph2$LL[1:length(names.b2)], n.round)
    dat.ph2$`est UL` <- round(out.ph2$UL[1:length(names.b2)], n.round)
    dat.ph2$`SE mean` <- round(out.ph2$mean[(length(names.b2) + 1):
                                              (2 * length(names.b2))], n.round)
    dat.ph2$`SE LL` <- round(out.ph2$LL[(length(names.b2) + 1):
                                          (2 * length(names.b2))], n.round)
    dat.ph2$`SE UL` <- round(out.ph2$UL[(length(names.b2) + 1):
                                          (2 * length(names.b2))], n.round)

    dat.ph3 <- expand.grid("comparison" = names.b1,
                           "group factor 2" = levels.b2)
    dat.ph3$`d mean` <- round(out.ph3$mean[((2 * length(levels.b2)) *
                                              length(names.b1) + 1):
                                             ((3 * length(levels.b2)) *
                                                length(names.b1))], n.round)
    dat.ph3$`d LL` <- round(out.ph3$LL[((2 * length(levels.b2)) *
                                          length(names.b2) + 1):
                                         ((3 * length(levels.b2)) *
                                            length(names.b2))], n.round)
    dat.ph3$`d UL` <- round(out.ph3$UL[((2 * length(levels.b2)) *
                                          length(names.b2) + 1):
                                         ((3 * length(levels.b2)) *
                                            length(names.b2))], n.round)
    dat.ph3$`est mean` <- round(out.ph3$mean[1:(length(levels.b2) *
                                                  length(names.b1))], n.round)
    dat.ph3$`est LL` <- round(out.ph3$LL[1:(length(levels.b2) *
                                              length(names.b1))], n.round)
    dat.ph3$`est UL` <- round(out.ph3$UL[1:(length(levels.b2) *
                                              length(names.b1))], n.round)
    dat.ph3$`SE mean` <- round(out.ph3$mean[(length(levels.b2) *
                                               length(names.b1) + 1):
                                              ((2 * length(levels.b2)) *
                                                 length(names.b1))], n.round)
    dat.ph3$`SE LL` <- round(out.ph3$LL[(length(levels.b2) * length(names.b1)
                                         + 1):((2 * length(levels.b2)) *
                                                 length(names.b1))], n.round)
    dat.ph3$`SE UL` <- round(out.ph3$UL[(length(levels.b2) * length(names.b1)
                                         + 1):((2 * length(levels.b2)) *
                                                 length(names.b1))], n.round)

    dat.ph3 <- dat.ph3[order(dat.ph3$`group factor 2`), c("group factor 2",
                                                          "comparison",
                                                          "d mean", "d LL",
                                                          "d UL", "est mean",
                                                          "est LL", "est UL",
                                                          "SE mean", "SE LL",
                                                          "SE UL")]

    dat.ph4 <- expand.grid("comparison" = names.b2,
                           "group factor 1" = levels.b1)
    dat.ph4$`d mean` <- round(out.ph4$mean[((2 * length(levels.b1)) *
                                              length(names.b2) + 1):
                                             ((3 * length(levels.b1)) *
                                                length(names.b2))], n.round)
    dat.ph4$`d LL` <- round(out.ph4$LL[((2 * length(levels.b1)) *
                                          length(names.b2)  + 1):
                                         ((3 * length(levels.b1)) *
                                            length(names.b2))], n.round)
    dat.ph4$`d UL` <- round(out.ph4$UL[((2 * length(levels.b1)) *
                                          length(names.b2) + 1):
                                         ((3 * length(levels.b1)) *
                                            length(names.b2))], n.round)
    dat.ph4$`est mean` <- round(out.ph4$mean[1:(length(levels.b1) *
                                                  length(names.b2))], n.round)
    dat.ph4$`est LL` <- round(out.ph4$LL[1:(length(levels.b1) *
                                              length(names.b2))], n.round)
    dat.ph4$`est UL` <- round(out.ph4$UL[1:(length(levels.b1) *
                                              length(names.b2))], n.round)
    dat.ph4$`SE mean` <- round(out.ph4$mean[(length(levels.b1) *
                                               length(names.b2) + 1):
                                              ((2 * length(levels.b1)) *
                                                 length(names.b2))], n.round)
    dat.ph4$`SE LL` <- round(out.ph4$LL[(length(levels.b1) * length(names.b2)
                                         + 1):((2 * length(levels.b1)) *
                                                 length(names.b2))], n.round)
    dat.ph4$`SE UL` <- round(out.ph4$UL[(length(levels.b1) * length(names.b2)
                                         + 1):((2 * length(levels.b1)) *
                                                 length(names.b2))], n.round)

    dat.ph4 <- dat.ph4[order(dat.ph4$`group factor 1`), c("group factor 1",
                                                          "comparison",
                                                          "d mean", "d LL",
                                                          "d UL", "est mean",
                                                          "est LL", "est UL",
                                                          "SE mean", "SE LL",
                                                          "SE UL")]

    output <- list(type.aov = "Two-way between ANOVA",
                   factor1 = levels.b1,
                   factor2 = levels.b2,
                   anova = round(orig.aov$anova_table, n.round),
                   em.1 = dat.em1,
                   no.test1 = no.test1,
                   ph.1 = dat.ph1,
                   em.2 = dat.em2,
                   no.test2 = no.test2,
                   ph.2 = dat.ph2,
                   em.3 = dat.em3,
                   no.test3 = no.test3,
                   ph.3 = dat.ph3,
                   em.4 = dat.em4,
                   no.test4 = no.test4,
                   ph.4 = dat.ph4)

    class(output) <- "AOboot_two"
    output

    }
}
