sample_bg_upd <- function (data, x, y, n, method = "random", rlayer, maskval = NULL, 
    calibarea = NULL, rbias = NULL, sp_name = NULL) 
{
    . <- NULL
    if (!method[1] %in% c("random", "thickening", "biased")) {
        stop("argument 'method' was misused, available methods 'random', 'thickening'")
    }
    if (method[1] %in% c("biased") & is.null(rbias)) {
        stop("for using 'random' method a raster layer with biases data must be provided in 'rbias' argument")
    }
    if (class(rlayer)[1] != "SpatRaster") {
        rlayer <- terra::rast(rlayer)
    }
    if (!is.null(rbias)) {
        if (class(rbias)[1] != "SpatRaster") 
            rbias <- terra::rast(rbias)
    }
    rlayer <- rlayer[[1]]
    data <- data[, c(x, y)]
    rlayer[na.omit(terra::cellFromXY(rlayer, as.matrix(data)))] <- NA
    if (!is.null(calibarea)) {
        rlayer <- rlayer %>% terra::crop(., calibarea) %>% terra::mask(., 
            calibarea)
    }
    if (!is.null(maskval)) {
        if (is.factor(maskval)) {
            maskval <- which(levels(maskval) %in% as.character(maskval))
            rlayer <- rlayer * 1
        }
        filt <- terra::match(rlayer, maskval)
        rlayer <- terra::mask(rlayer, filt)
    }
    if (method[1] %in% c("biased")) {
        if (any(!(ext(rlayer)[1:4] %in% ext(rbias)[1:4])) | all(!res(rlayer) %in% 
            res(rbias))) {
            if (!all(res(rlayer) %in% res(rbias))) {
                rbias <- terra::resample(rbias, rlayer, method = "bilinear")
            }
            rbias2 <- rbias %>% terra::crop(., rlayer) %>% terra::mask(., 
                rlayer)
            rbias <- rbias2
            rm(rbias2)
        }
    }
    if (method[1] %in% c("biased")) {
        rlayer <- mask(rlayer, rbias)
    }
    ncellr <- terra::global(!is.na(rlayer), sum)
    if (any(method == "thickening")) {
        data2 <- terra::vect(data, geom = c(x, y), crs = crs(rlayer))
        if (is.na(method["width"])) {
            buf_with <- mean(terra::distance(data2))
        } else {
            buf_with <- as.numeric(method["width"])
        }
        buf <- terra::buffer(data2, buf_with, quadsegs = 10)
        buf_r <- terra::rasterize(buf, rlayer, background = 0, fun = "count")
        buf_r <- terra::mask(buf_r, rlayer)
    }
    if (ncellr < n) {
        message("Number of background-points exceeds number of cell will be returned ", 
            ncellr, " background-points")
        cell_samp <- terra::as.data.frame(rlayer, na.rm = TRUE, 
            cells = TRUE)[, "cell"]
        cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>% 
            data.frame() %>% dplyr::tibble()
    } else {
        cell_samp <- terra::as.data.frame(rlayer, na.rm = TRUE, 
            cells = TRUE)[, "cell"]
        if (any(method == "random")) {
            cell_samp <- sample(cell_samp, size = n, replace = FALSE, 
                prob = NULL)
        } else if (any(method == "thickening")) {
            cell_samp <- sample(cell_samp, size = n, replace = FALSE, 
                prob = buf_r[cell_samp][, 1])
        } else if (any(method == "biased")) {
            cell_samp <- sample(cell_samp, size = n, replace = FALSE, 
                prob = rbias[cell_samp][, 1])
        }
        cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>% 
            data.frame() %>% dplyr::tibble()
    }
    colnames(cell_samp) <- c(x, y)
    cell_samp$pr_ab <- 0
    if (!is.null(sp_name)) {
        cell_samp <- tibble(sp = sp_name, cell_samp)
    }
    return(cell_samp)
}
