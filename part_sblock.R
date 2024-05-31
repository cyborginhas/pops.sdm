## %######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
## %######################################################%##


#' pre_tr_te
#'
#' @noRd
pre_tr_te <- function(data, p_names, h) {
    train <- list()
    test <- list()

    if (any(c("train", "train-test", "test")
    %in%
        unique(data[, p_names[h]]))) {
        np2 <- 1

        filt <- grepl("train", data[, p_names[h]])
        train[[1]] <- data[filt, ] %>%
            dplyr::select(-p_names[!p_names == p_names[h]])

        filt <- grepl("test", data[, p_names[h]])
        test[[1]] <- data[filt, ] %>%
            dplyr::select(-p_names[!p_names == p_names[h]])
    } else {
        np2 <- max(data[p_names[h]])

        for (i in 1:np2) {
            train[[i]] <- data[data[p_names[h]] != i, ] %>%
                dplyr::select(-p_names[!p_names == p_names[h]])

            test[[i]] <- data[data[p_names[h]] == i, ] %>%
                dplyr::select(-p_names[!p_names == p_names[h]])
        }
    }
    return(list(train = train, test = test, np2 = np2))
}



# Inverse bioclim
#'
#' @noRd
#'
bio <- function(data, env_layer) {
    . <- NULL
    if (class(data)[1] != "data.frame") {
        data <- data.frame(data)
    }
    if (!methods::is(env_layer, "SpatRaster")) {
        env_layer <- terra::rast(env_layer)
    }

    data <- na.omit(data)

    result <- env_layer[[1]]
    result[] <- NA

    minv <- apply(data, 2, min)
    maxv <- apply(data, 2, max)
    vnames <- names(data)

    data_2 <- data %>%
        na.omit() %>%
        apply(., 2, sort) %>%
        data.frame()

    rnk <- function(x, y) {
        b <- apply(y, 1, FUN = function(z) sum(x < z))
        t <- apply(y, 1, FUN = function(z) sum(x == z))
        r <- (b + 0.5 * t) / length(x)
        i <- which(r > 0.5)
        r[i] <- 1 - r[i]
        r * 2
    }

    var_df <- terra::as.data.frame(env_layer)
    var_df <- na.omit(var_df)

    k <- (apply(t(var_df) >= minv, 2, all) &
        apply(t(var_df) <= maxv, 2, all))

    for (j in vnames) {
        var_df[k, j] <- rnk(
            data_2[, j],
            var_df[k, j, drop = FALSE]
        )
    }
    var_df[!k, ] <- 0
    res <- apply(var_df, 1, min)
    result[as.numeric(names(res))] <- res
    return(result)
}

inv_bio <- function(e, p) {
    if (!methods::is(e, "SpatRaster")) {
        e <- terra::rast(e)
    }
    r <- bio(data = terra::extract(e, p)[-1], env_layer = e)
    r <- (r - terra::minmax(r)[1]) /
        (terra::minmax(r)[2] - terra::minmax(r)[1])
    r <- r <= 0.01 # environmental constrain
    r[which(r[, ] == FALSE)] <- NA
    return(r)
}


#' Inverse geo
#'
#' @noRd
#'
inv_geo <- function(e, p, d) {
    colnames(p) <- c("x", "y")
    p <- terra::vect(p, geom = c("x", "y"), crs = terra::crs(e))
    b <- terra::buffer(p, width = d)
    b <- terra::rasterize(b, e, background = 0)
    e <- terra::mask(e, b, maskvalues = 1)
    return(e)
}

#' Boyce
#'
#' @description This function calculate Boyce index performance metric. Codes were adapted from
#' enmSdm package.
#'
#' @noRd
boyce <- function(pres,
                  contrast,
                  n_bins = 101,
                  n_width = 0.1) {
    lowest <- min(c(pres, contrast), na.rm = TRUE)
    highest <- max(c(pres, contrast), na.rm = TRUE) + .Machine$double.eps
    window_width <- n_width * (highest - lowest)

    lows <- seq(lowest, highest - window_width, length.out = n_bins)
    highs <- seq(lowest + window_width + .Machine$double.eps, highest, length.out = n_bins)

    ## initiate variables to store predicted/expected (P/E) values
    freq_pres <- NA
    freq_contrast <- NA

    # tally proportion of test presences/background in each class
    for (i in 1:n_bins) {
        # number of presence predictions in a class
        freq_pres[i] <-
            sum(pres >= lows[i] & pres < highs[i], na.rm = TRUE)

        # number of background predictions in this class
        freq_contrast[i] <-
            sum(contrast >= lows[i] & contrast < highs[i], na.rm = TRUE)
    }

    # mean bin prediction
    mean_pred <- rowMeans(cbind(lows, highs))

    # add small number to each bin that has 0 background frequency but does have a presence frequency > 0
    if (any(freq_pres > 0 & freq_contrast == 0)) {
        small_value <- 0.5
        freq_contrast[freq_pres > 0 & freq_contrast == 0] <- small_value
    }

    # remove classes with 0 presence frequency
    if (any(freq_pres == 0)) {
        zeros <- which(freq_pres == 0)
        mean_pred[zeros] <- NA
        freq_pres[zeros] <- NA
        freq_contrast[zeros] <- NA
    }

    # remove classes with 0 background frequency
    if (any(0 %in% freq_contrast)) {
        zeros <- which(freq_pres == 0)
        mean_pred[zeros] <- NA
        freq_pres[zeros] <- NA
        freq_contrast[zeros] <- NA
    }

    P <- freq_pres / length(pres)
    E <- freq_contrast / length(contrast)
    PE <- P / E

    # remove NAs
    rm_nas <- stats::complete.cases(data.frame(mean_pred, PE))
    # mean_pred <- mean_pred[rm_nas]
    # PE <- PE[rm_nas]

    # calculate Boyce index
    result <- stats::cor(
        x = ifelse(is.na(mean_pred), 0, mean_pred),
        y = ifelse(is.na(PE), 0, PE), method = "spearman"
    )
    return(result)
}

# maxnet:::predict.maxnet()
#' Predict maxnet
#' @importFrom stats model.matrix
#' @noRd
predict_maxnet <- function(object, newdata, clamp = TRUE, type = c("link", "exponential", "cloglog", "logistic"), ...) {
    categoricalval <- function(x, category) {
        ifelse(x == category, 1, 0)
    }
    thresholdval <- function(x, knot) {
        ifelse(x >= knot, 1, 0)
    }
    hingeval <- function(x, min, max) {
        pmin(1, pmax(0, (x - min) / (max - min)))
    }

    if (clamp) {
        for (v in intersect(names(object$varmax), names(newdata))) {
            newdata[, v] <- pmin(
                pmax(newdata[, v], object$varmin[v]),
                object$varmax[v]
            )
        }
    }
    terms <- sub(
        "hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
        names(object$betas)
    )
    terms <- sub(
        "categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\"\\2\")",
        terms
    )
    terms <- sub(
        "thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
        terms
    )
    f <- formula(paste("~", paste(terms, collapse = " + "), "-1"))
    mm <- model.matrix(f, data.frame(newdata))
    if (clamp) {
        mm <- t(pmin(
            pmax(t(mm), object$featuremins[names(object$betas)]),
            object$featuremaxs[names(object$betas)]
        ))
    }
    link <- (mm %*% object$betas) + object$alpha
    type <- match.arg(type)

    if (type == "link") {
        return(link)
    }
    if (type == "exponential") {
        return(exp(link))
    }
    if (type == "cloglog") {
        return(1 - exp(0 - exp(object$entropy + link)))
    }
    if (type == "logistic") {
        return(1 / (1 + exp(-object$entropy - link)))
    }
}

#' Outliers with Reverse Jackknife
#'
#' @noRd
#'
rev_jack <- function(v) {
    v2 <- v
    v <- unique(v)
    lgh <- length(v) - 1
    t1 <- (0.95 * sqrt(length(v))) + 0.2
    x <- sort(v)
    y <- rep(0, lgh)
    for (i in seq_len(lgh)) {
        x1 <- x[i + 1]
        if (x[i] < mean(v)) {
            y[i] <- (x1 - x[i]) * (mean(v) - x[i])
        } else {
            y[i] <- (x1 - x[i]) * (x1 - mean(v))
        }
    }
    my <- mean(y)
    z <- y / (sqrt(sum((y - my)^2) / lgh))
    out <- rep(0, length(v2))
    if (any(z > t1)) {
        f <- which(z > t1)
        v <- x[f]
        if (v < median(x)) {
            xa <- (v2 <= v) * 1
            out <- out + xa
        }
        if (v > median(x)) {
            xb <- (v2 >= v) * 1
            out <- out + xb
        }
    } else {
        out <- out
    }
    return(which(out == 1))
}

#' Calculate amount of data for each training dataset in a given partition
#'
#' @noRd
#'
n_training <- function(data, partition) {
    . <- partt <- NULL
    if (any(c("train", "train-test", "test")
    %in%
        (data %>%
            dplyr::select(dplyr::starts_with({{ partition }})) %>%
            dplyr::pull() %>%
            unique()))) {
        nn_part <- data %>%
            dplyr::select(dplyr::starts_with({{ partition }})) %>%
            apply(., 2, table) %>%
            data.frame()
        nn_part <- nn_part %>% dplyr::mutate(partt = rownames(nn_part))
        nn_part$partt[grepl("train", nn_part$partt)] <- "train"
        nn_part <- nn_part %>%
            dplyr::filter(partt == "train") %>%
            dplyr::select(-partt)
        nn_part <- colSums(nn_part)
    } else {
        data <- data %>%
            dplyr::select(dplyr::starts_with({{ partition }}))

        nn_part <- list()
        for (ppp in 1:ncol(data)) {
            nn_part[[ppp]] <- data %>%
                dplyr::pull(ppp) %>%
                table() %>%
                c()
            sm <- nn_part[[ppp]] %>% sum()
            nn_part[[ppp]] <- sapply(nn_part[[ppp]], function(x) sum(sm - x))
        }
        nn_part <- unlist(nn_part)
    }
    return(nn_part)
}

#' Calculate number of coefficient for gam models
#'
#' @noRd
#'
n_coefficients <- function(data, predictors, predictors_f = NULL, k = 10) {
    data <- data.frame(data)
    if (k < 0) {
        k <- 10
    }
    if (!is.null(predictors_f)) {
        n_levels <- rep(NA, length(predictors_f))
        for (fff in 1:length(predictors_f)) {
            n_levels[fff] <- unique(data[, predictors_f]) %>%
                na.omit() %>%
                length()
        }
        n_levels <- sum(n_levels)
    } else {
        n_levels <- 0
    }
    n <- (k - 1) * length(predictors) + n_levels
    return(n)
}

#' Euclidean distance for extrapolation
#'
#' @noRd
#'
euc_dist <- function(x, y) {
    if (!methods::is(x, "matrix")) {
        x <- as.matrix(x)
    }
    if (!methods::is(y, "matrix")) {
        y <- as.matrix(y)
    }
    result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
    for (ii in 1:nrow(y)) {
        result[, ii] <- sqrt(colSums((t(x) - y[ii, ])^2))
    }
    rownames(result) <- rownames(x)
    colnames(result) <- rownames(y)
    return(result)
}

#' Moran I, based on ape package
#'
#' @noRd
#'
morani <- function(x, weight, na.rm = FALSE, scaled = TRUE) {
    if (dim(weight)[1] != dim(weight)[2]) {
        stop("'weight' must be a square matrix")
    }
    n <- length(x)
    if (dim(weight)[1] != n) {
        stop("'weight' must have as many rows as observations in 'x'")
    }
    ei <- -1 / (n - 1)
    nas <- is.na(x)
    if (any(nas)) {
        if (na.rm) {
            x <- x[!nas]
            n <- length(x)
            weight <- weight[!nas, !nas]
        } else {
            warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
            return(NA)
        }
    }
    rs <- rowSums(weight)
    rs[rs == 0] <- 1
    weight <- weight / rs
    s <- sum(weight)
    m <- mean(x)
    y <- x - m
    cv <- sum(weight * y %o% y)
    v <- sum(y^2)
    res <- (n / s) * (cv / v)
    if (scaled) {
        imax <- (n / s) * (sd(rowSums(weight) * y) / sqrt(v / (n - 1)))
        res <- res / imax
    }

    return(res)
}

#' Mahalanobis distance
#'
#' @noRd
#'
mah_dist <- function(x, y, cov) {
    if (!methods::is(x, "matrix")) {
        x <- as.matrix(x)
    }
    if (!methods::is(y, "matrix")) {
        y <- as.matrix(y)
    }
    result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
    for (ii in 1:nrow(y)) {
        # root square of squared Mahalanobis distance
        result[, ii] <- sqrt(mahalanobis(x = x, center = y[ii, ], cov = cov))
    }
    rownames(result) <- rownames(x)
    colnames(result) <- rownames(y)
    return(result)
}

#' @description Debugged version of the function 'spart_block' from the package
#' flexsdm. This function is used to partition the data into blocks based
#' on the spatial autocorrelation, and environmental similarity.

part_sblock_upd <- function(
    env_layer, data, x, y, pr_ab, n_part = 3, min_res_mult = 3,
    max_res_mult = 200, num_grids = 30, min_occ = 10, prop = 0.5) {
    data <- dplyr::tibble(data)
    data <- data[, c(pr_ab, x, y)]
    colnames(data) <- c("pr_ab", "x", "y")
    if (any(!unique(data[, "pr_ab"][[1]]) %in% c(0, 1))) {
        stop(
            "values in pr_ab column did not match with 0 and 1:\nunique list values in pr_ab column are: ",
            paste(unique(data[, "pr_ab"]), collapse = " ")
        )
    }
    data <- dplyr::tibble(data, terra::extract(env_layer, data[, 2:3])[-1])
    filt <- stats::complete.cases(data)
    if (sum(!filt) > 0) {
        data <- data[filt, ]
        message(sum(!filt), " rows were excluded from database because NAs were found")
    }
    rm(filt)
    pa <- data %>% dplyr::pull(pr_ab)
    cell_size <- seq(terra::res(env_layer[[1]])[1] * min_res_mult,
        terra::res(env_layer[[1]])[1] * max_res_mult,
        length.out = num_grids
    )
    message(
        "The following grid cell sizes will be tested:\n",
        paste(round(cell_size, 2), collapse = " | "), "\n"
    )
    message("Creating basic raster mask...\n")
    mask <- env_layer[[1]]
    names(mask) <- "group"
    mask[!is.na(mask)] <- 1
    e <- terra::ext(mask)
    message("Searching for the optimal grid size...\n")
    mask2 <- mask
    mask2[] <- 0
    presences2 <- data
    presences2 <- terra::vect(presences2,
        geom = c("x", "y"),
        crs = terra::crs(mask)
    )
    grid <- list()
    DIM <- matrix(0, length(cell_size), 2)
    colnames(DIM) <- c("R", "C")
    for (i in 1:length(cell_size)) {
        mask3 <- mask2
        terra::res(mask3) <- cell_size[i]
        mask3 <- terra::extend(mask3, y = c(1, 1))
        DIM[i, ] <- dim(mask3)[1:2]
        terra::values(mask3) <- 1
        NAS <- c(terra::extract(mask3, presences2)[-1])
        if (any(is.na(NAS))) {
            while (any(is.na(NAS))) {
                terra::ext(mask3) <- terra::ext(mask3) + cell_size[i]
                terra::res(mask3) <- cell_size[i]
                DIM[i, ] <- dim(mask3)[1:2]
                terra::values(mask3) <- 1
                NAS <- terra::extract(mask3, presences2)[-1]
            }
        }
        grid[[i]] <- mask3
    }
    rm(list = c("mask3", "mask2", "mask"))
    for (i in 1:length(grid)) {
        if (n_part %% 2 == 0) {
            group <- c(
                rep(1:n_part, DIM[i, 2])[1:DIM[i, 2]],
                rep(c((n_part / 2 + 1):n_part, 1:(n_part / 2)), DIM[
                    i,
                    2
                ])[1:DIM[i, 2]]
            )
        }
        if (n_part %% 2 == 1) {
            group <- c(
                rep(1:n_part, DIM[i, 2])[1:DIM[i, 2]],
                rep(
                    c(as.integer(n_part / 2 + 1):n_part, 1:(n_part / 2)),
                    DIM[i, 2]
                )[1:DIM[i, 2]]
            )
        }
        terra::values(grid[[i]]) <- rep(group, length.out = terra::ncell(grid[[i]]))
    }
    part <- data.frame(matrix(0, nrow(presences2), length(grid)))
    for (i in 1:length(grid)) {
        part[, i] <- terra::extract(grid[[i]], presences2)[
            ,
            2
        ]
    }
    part <- dplyr::tibble(part)
    pp <- sapply(part[pa == 1, ], function(x) {
        length(unique(x))
    })
    pp <- pp == n_part
    pf <- sapply(part[pa == 1, ], table)
    if (is.list(pf) == TRUE) {
        pf <- which(sapply(pf, min) < min_occ)
    } else {
        pf <- which(apply(pf, 2, min) < min_occ)
    }
    pp[pf] <- FALSE
    cell_size <- cell_size[pp]
    grid <- grid[pp]
    part <- part[, pp]
    names(part) <- names(which(pp == TRUE))
    if (any(unique(pa) == 0)) {
        pa <- presences2$pr_ab
        pp <- sapply(part[pa == 0, ], function(x) {
            length(unique(x))
        })
        pp <- pp == n_part
        pf <- sapply(part[pa == 0, ], table)
        if (is.list(pf) == TRUE) {
            pf <- which(sapply(pf, min) < min_occ)
        } else {
            pf <- which(apply(pf, 2, min) < min_occ)
        }
        pp[pf] <- FALSE
        cell_size <- cell_size[pp]
        grid <- grid[pp]
        part <- part[, pp]
        names(part) <- names(which(pp == TRUE))
    }
    if (ncol(part) == 0) {
        message("It was not possible to find a good partition. Try to change values in 'n_part', or in 'min_res_mult', 'max_res_mult', or 'num_grids'")
        return(NA)
    }
    ncell <- data.frame(matrix(0, nrow(presences2), length(grid)))
    for (i in 1:length(grid)) {
        ncell[, i] <- terra::cellFromXY(grid[[i]], terra::geom(presences2)[
            ,
            c("x", "y")
        ])
    }
    sd_p <- rep(NA, length(grid))
    if (any(unique(pa) == 0)) {
        sd_a <- rep(NA, length(grid))
    }
    for (i in 1:ncol(part)) {
        if (any(unique(pa) == 0)) {
            sd_a[i] <- stats::sd(table(part[pa == 0, i]))
        }
        sd_p[i] <- stats::sd(table(part[pa == 1, i]))
    }
    Env.P <- terra::extract(env_layer, presences2)[-1]
    env_sim <- rep(NA, length(grid))
    for (i in 1:ncol(part)) {
        cmb <- unique(part[, i][[1]]) %>% combn(2)
        Env.P1 <- cbind(part[i], Env.P)
        Env.P1 <- Env.P1[complete.cases(Env.P1), ]
        Env.P1 <- split(Env.P1[, -1], Env.P1[, 1])
        euq_c <- list()
        for (r in 1:ncol(cmb)) {
            euq_c[[r]] <- euc_dist(Env.P1[[cmb[1, r]]], Env.P1[[cmb[2, r]]]) %>% mean()
        }
        env_sim[i] <- euq_c %>%
            unlist() %>%
            mean()
        rm(list = c("Env.P1"))
    }
    spa_auto <- rep(NA, length(grid))
    presences2 <- terra::geom(presences2)[, c("x", "y")] %>%
        as.data.frame()
    dist <- euc_dist(presences2, presences2)
    dist <- 1 / dist
    diag(dist) <- 0
    dist[which(dist == Inf)] <- 0
    for (p in 1:ncol(part)) {
        cmb <- unique(part[, p][[1]]) %>% utils::combn(2)
        imoran_grid_c <- rep(NA, ncol(cmb))
        dff <- dplyr::tibble(
            nrow = 1:nrow(part), data["pr_ab"],
            group = part[p][[1]]
        )
        for (c in 1:ncol(cmb)) {
            filt <- dff %>%
                dplyr::group_by(group, pr_ab) %>%
                dplyr::slice_sample(prop = prop) %>%
                dplyr::pull(nrow) %>%
                sort()
            odd <- which((part[p][[1]] == cmb[1, c])[filt])
            even <- which((part[p][[1]] == cmb[2, c])[filt])
            dist2 <- dist[filt, filt]
            dist2[odd, odd] <- 0
            dist2[even, even] <- 0
            mins <- apply(dist2, 2, function(x) {
                max(x, na.rm = TRUE)
            })
            for (i in 1:length(mins)) {
                dist2[, i] <- ifelse(dist2[, i] == mins[i], mins[i],
                    0
                )
            }
            if (nrow(data) < 3) {
                imoran_grid_c[c] <- NA
            } else {
                im <- sapply(data[filt, names(env_layer)], function(x) {
                    suppressMessages(morani(x, dist2,
                        na.rm = TRUE,
                        scaled = TRUE
                    ))
                })
                imoran_grid_c[c] <- mean(abs(im))
            }
        }
        spa_auto[p] <- mean(imoran_grid_c)
    }
    Opt <- if (any(unique(pa) == 0)) {
        data.frame(
            n_grid = 1:length(cell_size), cell_size = cell_size,
            round(
                data.frame(spa_auto, env_sim, sd_p, sd_a),
                3
            )
        )
    } else {
        data.frame(
            n_grid = 1:length(cell_size), cell_size = cell_size,
            round(data.frame(spa_auto, env_sim, sd_p), 3)
        )
    }
    Opt2 <- Opt
    rownames(Opt2) <- colnames(part)
    Dup <- if (any(unique(pa) == 0)) {
        !duplicated(Opt2[c("spa_auto", "env_sim", "sd_p", "sd_a")])
    } else {
        !duplicated(Opt2[c("spa_auto", "env_sim", "sd_p")])
    }
    Opt2 <- Opt2[Dup, ]
    while (nrow(Opt2) > 1) {
        if (nrow(Opt2) == 1) {
            break
        }
        Opt2 <- Opt2[which(Opt2$spa_auto <= summary(Opt2$spa_auto)[2]), ]
        if (nrow(Opt2) == 1) {
            break
        }
        Opt2 <- Opt2[which(Opt2$env_sim >= summary(Opt2$env_sim)[5]), ]
        if (nrow(Opt2) == 1) {
            break
        }
        Opt2 <- Opt2[which(Opt2$sd_p <= summary(Opt2$sd_p)[2]), ]
        if (nrow(Opt2) == 2) {
            break
        }
        if (any(unique(pa) == 0)) {
            Opt2 <- Opt2[which(Opt2$sd_a <= summary(Opt2$sd_a)[2]), ]
            if (nrow(Opt2) == 2) {
                break
            }
        }
        if (unique(Opt2$spa_auto) && unique(Opt2$env_sim) &&
            unique(Opt2$sd_p)) {
            Opt2 <- Opt2[nrow(Opt2), ]
        }
    }
    if (nrow(Opt2) > 1) {
        Opt2 <- Opt2[nrow(Opt2), ]
    }
    result <- data.frame(data, .part = c(part[, rownames(Opt2)])[[1]])
    result <- result %>% dplyr::select(-names(env_layer))
    colnames(result) <- c("pr_ab", "x", "y", ".part")
    result <- result[c("x", "y", "pr_ab", ".part")]
    grid <- grid[[Opt2$n_grid]]
    names(grid) <- ".part"
    out <- list(
        part = dplyr::tibble(result), best_part_info = dplyr::tibble(Opt2),
        grid = grid
    )
    return(out)
}
