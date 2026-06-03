# Dependencies should be declared in DESCRIPTION under Imports:
# data.table, lubridate, rgbif, rinat, BIEN, neonUtilities

#' @description Generic function to print messages
#' @param message Character string of message to be printed.
#' @return Prints message to console.
#' @keywords internal
progress <- function(message) {
  print(paste0(message))
}

#' @description Convert all columns in a data.table to character
#' @param data A data.table.
#' @return The input data.table with all columns converted to character.
#' @keywords internal
as_character_dt <- function(data) {
  for (col in names(data)) {
    data.table::set(data, j = col, value = as.character(data[[col]]))
  }

  data
}

#' @description Convert all data.tables in an occurrence list to character
#' @param occs A list of data.tables, raw data, metadata, & cleaned data.
#' @return A list of data.tables with all columns converted to character.
#' @keywords internal
as_character_list <- function(occs) {
  lapply(occs, as_character_dt)
}

#' @description Create an empty occurrence output list
#' @return A list of empty data.tables for raw data, metadata, & cleaned data.
#' @keywords internal
empty_occ_list <- function() {
  list(
    data.table::data.table(),
    data.table::data.table(),
    data.table::data.table()
  )
}

#' @description Create temporary occurrence cache file path
#' @param tmp_dir Character string giving the temporary cache directory.
#' @param sp Character string giving the species name.
#' @param date_range Character string giving the date range.
#' @param db Character string identifying the database.
#' @return Character string giving the temporary cache file path.
#' @keywords internal
make_tmp_file <- function(tmp_dir, sp, date_range, db) {
  date_range_fn <- gsub(",", "_", gsub("-", "", date_range))

  file.path(
    tmp_dir,
    paste0(gsub(" ", "_", sp), "_", date_range_fn, "_", db, ".rds")
  )
}

#' @description Flatten nested occurrence download results
#' @param x A nested list of data.tables, data.frames, or NULL values.
#' @return A flat list of data.tables or data.frames.
#' @keywords internal
flatten_occ_results <- function(x) {
  if (is.null(x)) {
    return(list())
  }

  if (inherits(x, "data.frame")) {
    return(list(x))
  }

  if (is.list(x)) {
    return(unlist(lapply(x, flatten_occ_results), recursive = FALSE))
  }

  list()
}

#' @description Row-bind nested occurrence download results
#' @param x A nested list of data.tables, data.frames, or NULL values.
#' @return A data.table.
#' @keywords internal
rbind_occ_results <- function(x) {
  x <- flatten_occ_results(x)

  if (!length(x)) {
    return(data.table::data.table())
  }

  data.table::rbindlist(x, fill = TRUE)
}

#' @description Finalize occurrence output
#' @param occs List containing raw data, metadata, and cleaned data.
#' @param taxon_key A data.table returned from `batch_upd_dates()`.
#' @param path Character string to the Data folder.
#' @param tmp_file Character string giving the cache file path.
#' @param save_cache Logical. If TRUE, save occurrence output to `tmp_file`.
#' @return Finalized occurrence list.
#' @keywords internal
finish_occurrence_output <- function(occs, taxon_key, path, tmp_file,
                                     save_cache = TRUE) {
  occs <- append_dirs(occs, taxon_key, path)

  if (nrow(occs[[2]]) == 0) {
    clean_file <- paste0(path, names(occs)[3])

    if (file.exists(clean_file)) {
      occs[[3]] <- data.table::fread(
        clean_file,
        colClasses = "character"
      )
    }
  }

  if (save_cache) {
    saveRDS(occs, tmp_file)
  }

  occs
}

#' @description Function to get the GBIF taxon key for a given taxon and rank
#' @param scientific_name Character string of taxon name.
#' @param taxon_rank Character string of rank. Must be `"GENUS"` or `"SPECIES"`.
#' @return A data.table of the taxon key for the given taxon and rank.
#' @export
get_gbif_keys <- function(scientific_name, taxon_rank) {
  taxon_rank <- toupper(taxon_rank)

  if (taxon_rank == "GENUS") {
    keys <- rgbif::name_backbone(
      name = scientific_name,
      rank = taxon_rank,
      strict = TRUE
    )$usageKey

    species <- rgbif::name_lookup(
      higherTaxonKey = keys,
      status = "ACCEPTED",
      rank = "SPECIES",
      limit = 99999
    )$data

    species <- data.table::as.data.table(species)
    species <- species[, c("species", "speciesKey", "genusKey"), with = FALSE]
  } else if (taxon_rank == "SPECIES") {
    keys <- rgbif::name_backbone(
      name = scientific_name,
      rank = taxon_rank,
      strict = TRUE
    )

    species <- data.table::as.data.table(keys)
    species <- species[, c("species", "speciesKey", "genusKey"), with = FALSE]
  } else {
    stop("confirm that taxon rank is 'GENUS' or 'SPECIES'")
  }

  species
}

#' @description Add date ranges to taxon keys based on existing metadata
#' @param taxon_key data.table of taxon keys & epithets from `get_gbif_keys()`.
#' @param path Character string to the Data folder.
#' @return A data.table of taxon keys with date ranges added for each database.
#' @export
batch_upd_dates <- function(taxon_key, path) {
  sp <- taxon_key$species
  key <- taxon_key$speciesKey

  metadata_path <- paste0(
    path,
    "Original/host-occurrences/hostpts_metadata.csv"
  )

  if (file.exists(metadata_path)) {
    metadata <- data.table::fread(metadata_path, colClasses = "character")
  } else {
    metadata <- data.table::data.table(
      speciesKey = character(),
      species = character(),
      date_range = character(),
      database = character(),
      dateretrieved = character(),
      start_date = character(),
      end_date = character()
    )
  }

  metadata$start_date <- lubridate::ymd(metadata$start_date)
  metadata$end_date <- lubridate::ymd(metadata$end_date)

  db_values <- c("gbif", "inat", "bien", "neon")

  result <- lapply(db_values, function(db) {
    metadata_db <- metadata[
      metadata$speciesKey == key & metadata$database == db,
    ]

    if (nrow(metadata_db) > 0) {
      date_range <- paste0(metadata_db$end_date + 1, ",", Sys.Date())
      dateretrieved <- metadata_db$dateretrieved
    } else {
      dateretrieved <- NA
      date_range <- switch(db,
        gbif = paste0("1600-01-01", ",", Sys.Date()),
        inat = paste0("2008-01-01", ",", Sys.Date()),
        bien = paste0("1600-01-01", ",", Sys.Date()),
        neon = paste0("2013-05-31", ",", Sys.Date()),
        stop("confirm that database is 'gbif', 'inat', 'neon', 'bien'")
      )
    }

    data.table::data.table(
      speciesKey = key,
      species = sp,
      date_range = date_range,
      database = db,
      dateretrieved = dateretrieved
    )
  })

  data.table::rbindlist(result, fill = TRUE)
}

#' @description Create date intervals to avoid API limits
#' @param date_range Character string of start and end dates separated by comma.
#' @param db Character string of database name.
#' @return A vector or list of dates.
#' @keywords internal
get_dates <- function(date_range, db) {
  start_date <- as.Date(strsplit(date_range, ",")[[1]][1])
  end_date <- as.Date(strsplit(date_range, ",")[[1]][2])

  if (start_date == "1600-01-01") {
    dates <- unique(c(
      start_date,
      as.Date("1994-12-31"),
      seq(as.Date("1994-12-31"), as.Date("2014-12-31"), by = "8 month"),
      seq(as.Date("2014-12-31"), Sys.Date(), by = "4 month"),
      Sys.Date()
    ))
  } else if (start_date == "2008-01-01") {
    dates <- list(
      years = seq(as.Date("2008-01-01"), as.Date("2017-01-01"), by = "year"),
      months = NULL,
      days = NULL
    )

    dates[[2]] <- seq(as.Date("2018-01-01"), Sys.Date(), by = "month")
    dates[[3]] <- seq(dates[[2]][length(dates[[2]])], Sys.Date(), by = "day")
    dates[[2]] <- dates[[2]][-length(dates[[2]])]
  } else if (db == "gbif") {
    dates <- seq(start_date + 1, end_date, by = "4 month")
    dates <- sort(c(dates, end_date))
  } else {
    dates <- seq(start_date, end_date, by = "day")
  }

  dates
}

#' @description Split an iNat query period into smaller date periods
#' @param date Date value representing the current query period.
#' @param date_unit Character string giving the current date resolution.
#' @param date_range Character string giving the full allowed date range.
#' @return A list with `dates` and `date_unit` for the next smaller resolution.
#' @keywords internal
split_inat_period <- function(date, date_unit, date_range) {
  range <- as.Date(strsplit(date_range, ",")[[1]])
  start_date <- range[1]
  end_date <- range[2]

  if (date_unit == "year") {
    dates <- seq(
      as.Date(paste0(lubridate::year(date), "-01-01")),
      as.Date(paste0(lubridate::year(date), "-12-01")),
      by = "month"
    )

    dates <- dates[
      dates >= as.Date(format(start_date, "%Y-%m-01")) &
        dates <= as.Date(format(end_date, "%Y-%m-01"))
    ]

    return(list(dates = dates, date_unit = "month"))
  }

  if (date_unit == "month") {
    month_start <- as.Date(format(date, "%Y-%m-01"))
    month_end <- lubridate::ceiling_date(month_start, "month") - 1

    dates <- seq(
      max(month_start, start_date),
      min(month_end, end_date),
      by = "day"
    )

    return(list(dates = dates, date_unit = "day"))
  }

  stop("Too many records for a single day. Cannot reduce date range further.")
}

#' @description Split a GBIF query period into smaller date periods
#' @param start_date Start date of the current query period.
#' @param end_date End date of the current query period.
#' @return A data.table with `start_date` and `end_date` columns.
#' @keywords internal
split_gbif_period <- function(start_date, end_date) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  if (start_date > end_date) {
    return(data.table::data.table())
  }

  if (as.integer(end_date - start_date) > 31) {
    month_starts <- seq(
      as.Date(format(start_date, "%Y-%m-01")),
      as.Date(format(end_date, "%Y-%m-01")),
      by = "month"
    )

    periods <- data.table::data.table(
      start_date = pmax(month_starts, start_date),
      end_date = pmin(
        lubridate::ceiling_date(month_starts, "month") - 1,
        end_date
      )
    )

    return(periods[periods$start_date <= periods$end_date])
  }

  if (start_date < end_date) {
    days <- seq(start_date, end_date, by = "day")

    return(data.table::data.table(
      start_date = days,
      end_date = days
    ))
  }

  stop("Too many GBIF records. Cannot reduce date range further.")
}

#' @description Add output file paths to occurrence list
#' @param occs List of occurrence data.
#' @param taxon_key A data.table of taxon key, species, database, & date range.
#' @param path Character string to the Data folder.
#' @return List of occurrence data with output paths as names.
#' @keywords internal
append_dirs <- function(occs, taxon_key, path) {
  sp <- taxon_key$species
  db <- taxon_key$database
  species_fn <- gsub(" ", "_", sp)

  if (db == "neon") {
    filenames <- c(
      "Original/host-occurrences/plot_db/neon/neon_raw.csv",
      "Original/host-occurrences/hostpts_metadata.csv",
      paste0(
        "Table/USA/host-occurrences/", species_fn, "/",
        species_fn, "_", db, "_clean.csv"
      )
    )
  } else {
    filenames <- c(
      paste0(
        "Original/host-occurrences/", species_fn, "/",
        species_fn, "_", db, "_raw.csv"
      ),
      "Original/host-occurrences/hostpts_metadata.csv",
      paste0(
        "Table/Global/host-occurrences/", species_fn, "/",
        species_fn, "_", db, "_clean.csv"
      )
    )
  }

  lapply(filenames, function(filename) {
    dir.create(
      paste0(path, dirname(filename)),
      recursive = TRUE,
      showWarnings = FALSE
    )
  })

  names(occs) <- filenames
  occs
}

#' @description Copy data from cloud to local folder and set permissions
#' @param cloud_path Character string to the cloud Data folder.
#' @param local_path Character string to the local Data folder.
#' @return Invisibly, copies directories and removes `.DS_Store` files.
#' @export
copy_existing_data <- function(cloud_path, local_path) {
  dirs_to_copy <- c(
    "Original/host-occurrences/",
    "Table/Global/host-occurrences/",
    "Table/USA/host-occurrences/"
  )

  lapply(dirs_to_copy, function(dir) {
    cloud_dir <- file.path(cloud_path, dir)
    local_dir <- file.path(local_path, dir)

    if (dir.exists(cloud_dir)) {
      dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(cloud_dir, dirname(local_dir), recursive = TRUE)

      Sys.chmod(
        list.files(local_dir, full.names = TRUE, recursive = TRUE),
        mode = "0777",
        use_umask = FALSE
      )
    }
  })

  hidden_files <- list.files(
    local_path,
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE
  )

  ds_store_files <- hidden_files[grepl(".DS_Store", hidden_files)]

  if (length(ds_store_files) > 0) {
    Sys.chmod(ds_store_files, mode = "0777", use_umask = FALSE)
    invisible(file.remove(ds_store_files))
  }
}

#' @description Read cached occurrence data if available
#' @param tmp_file Character string giving the expected cache file path.
#' @param tmp_dir Character string giving the temporary cache directory.
#' @param sp Character string giving the species name.
#' @param start_date Start date for the queried period.
#' @param end_date End date for the queried period.
#' @param db Character string identifying the database/source.
#' @param progress Function used to print progress messages.
#' @return Cached occurrence list read from an `.rds` file, or `NULL`.
#' @keywords internal
read_cached_occurrences <- function(
    tmp_file, tmp_dir, sp,
    start_date, end_date, db, progress = message) {
  cache_file <- if (file.exists(tmp_file)) {
    tmp_file
  } else {
    files <- list.files(
      tmp_dir,
      full.names = TRUE,
      pattern = paste0("^", gsub(" ", "_", sp), ".*_", db, "\\.rds$")
    )

    if (!length(files)) {
      return(NULL)
    }

    files[which.max(file.info(files)$ctime)]
  }

  progress(sprintf(
    "Reading in existing %s data for %s...",
    db,
    sp
  ))

  readRDS(cache_file)
}

#' @description Count occurrence records and report progress
#' @param db Character string identifying the database/source.
#' @param sp Character string giving the species name.
#' @param key GBIF species key. Required when `db = "gbif"`.
#' @param date Optional date filter.
#' @param date_range Optional character string used in progress messages.
#' @param date_unit Character string giving the iNat date resolution.
#' @param data Optional downloaded iNat data object used to report count.
#' @param progress Function used to print progress messages.
#' @return Integer count of occurrence records.
#' @keywords internal
occ_count_progress <- function(db, sp, key = NULL, date = NULL,
                               date_range = NULL, date_unit = NULL,
                               data = NULL, progress = message) {
  db <- match.arg(db, c("gbif", "inat"))

  if (db == "inat" && !is.null(data)) {
    period <- switch(date_unit,
      year = lubridate::year(date),
      month = paste0(lubridate::month(date), "-", lubridate::year(date)),
      day = as.character(date)
    )

    progress(sprintf(
      "Downloaded iNat records (n: %s) for %s",
      nrow(data),
      period
    ))

    return(nrow(data))
  }

  if (db == "gbif") {
    args <- list(
      speciesKey = key,
      hasCoordinate = TRUE,
      hasGeospatialIssue = FALSE
    )

    if (length(date) == 2) {
      args$eventDate <- paste(date, collapse = ",")
    }

    count <- do.call(rgbif::occ_count, args)

    absent_count <- do.call(
      rgbif::occ_count,
      c(args, list(occurrenceStatus = "ABSENT"))
    )
  }

  if (db == "inat") {
    args <- list(
      taxon_name = sp,
      quality = "research",
      geo = TRUE,
      maxresults = 0,
      meta = TRUE
    )

    if (!is.null(date)) {
      args$year <- lubridate::year(date)

      if (date_unit %in% c("month", "day")) {
        args$month <- lubridate::month(date)
      }

      if (date_unit == "day") {
        args$day <- lubridate::day(date)
      }
    }

    count <- do.call(rinat::get_inat_obs, args)$meta$found
  }

  if (is.null(date)) {
    if (db == "gbif") {
      progress(sprintf(
        "Downloading %s %s (%s): %s total presences, %s total absences...",
        db, sp, date_range, count, absent_count
      ))
    } else {
      progress(sprintf(
        "Downloading %s %s (%s): %s total presences...",
        db, sp, date_range, count
      ))
    }
  } else if (db == "gbif") {
    progress(sprintf(
      "%s: %s presences, %s absences",
      date_range, count, absent_count
    ))
  }

  count
}

#' @description Download iNat data for a given species over a date period
#' @param sp Character string giving the species name.
#' @param dates Date vector to query.
#' @param date_unit Character string giving the date resolution.
#' @param date_range Character string used in progress messages.
#' @param progress Function used to print progress messages.
#' @return A list of iNat occurrence data objects.
#' @keywords internal
get_inat_period <- function(sp, dates, date_unit, date_range,
                            progress = message) {
  lapply(dates, function(x) {
    tryCatch(
      {
        count <- occ_count_progress(
          db = "inat",
          sp = sp,
          date = x,
          date_range = date_range,
          date_unit = date_unit,
          progress = progress
        )

        Sys.sleep(1)

        if (count == 0) {
          return(NULL)
        }

        if (count > 9999) {
          split_dates <- split_inat_period(
            date = x,
            date_unit = date_unit,
            date_range = date_range
          )

          progress(sprintf(
            "Too many iNat records for %s; splitting to %s...",
            switch(date_unit,
              year = lubridate::year(x),
              month = paste0(lubridate::month(x), "-", lubridate::year(x)),
              day = as.character(x)
            ),
            split_dates$date_unit
          ))

          return(get_inat_period(
            sp = sp,
            dates = split_dates$dates,
            date_unit = split_dates$date_unit,
            date_range = date_range,
            progress = progress
          ))
        }

        args <- list(
          taxon_name = sp,
          quality = "research",
          geo = TRUE,
          maxresults = 9999,
          meta = FALSE,
          year = lubridate::year(x)
        )

        if (date_unit %in% c("month", "day")) {
          args$month <- lubridate::month(x)
        }

        if (date_unit == "day") {
          args$day <- lubridate::day(x)
        }

        data <- do.call(rinat::get_inat_obs, args)

        occ_count_progress(
          db = "inat",
          sp = sp,
          date = x,
          date_range = date_range,
          date_unit = date_unit,
          data = data,
          progress = progress
        )

        data
      },
      error = function(e) {
        if (grepl("Cannot reduce date range further", conditionMessage(e))) {
          stop(e)
        }

        NULL
      }
    )
  })
}

#' @description Download GBIF occurrence data for a period, splitting if needed
#' @param key GBIF species key.
#' @param sp Character string giving the species name.
#' @param start_date Start date of the query period.
#' @param end_date End date of the query period.
#' @param occurrence_status Occurrence status. Use `"ABSENT"` for absences.
#' @param publishing_org Optional GBIF publishing organization key.
#' @param progress Function used to print progress messages.
#' @return A data.table, nested list of data.tables, or `NULL`.
#' @keywords internal
get_gbif_period <- function(key, sp, start_date, end_date,
                            occurrence_status = NULL,
                            publishing_org = NULL,
                            progress = message) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  period <- paste0(start_date, ",", end_date)

  count_args <- list(
    speciesKey = key,
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE,
    eventDate = period
  )

  if (!is.null(occurrence_status)) {
    count_args$occurrenceStatus <- occurrence_status
  }

  if (!is.null(publishing_org)) {
    count_args$publishingOrg <- publishing_org
  }

  count <- do.call(rgbif::occ_count, count_args)

  if (is.null(occurrence_status)) {
    progress(sprintf("%s: %s presences", period, count))
  } else {
    progress(sprintf(
      "%s: %s %s records",
      period,
      count,
      tolower(occurrence_status)
    ))
  }

  if (count == 0) {
    return(NULL)
  }

  if (count > 99999) {
    split_dates <- split_gbif_period(start_date, end_date)

    progress(sprintf(
      "Too many GBIF records for %s; splitting into %s smaller periods...",
      period,
      nrow(split_dates)
    ))

    return(lapply(seq_len(nrow(split_dates)), function(i) {
      Sys.sleep(1)

      get_gbif_period(
        key = key,
        sp = sp,
        start_date = split_dates$start_date[i],
        end_date = split_dates$end_date[i],
        occurrence_status = occurrence_status,
        publishing_org = publishing_org,
        progress = progress
      )
    }))
  }

  data_args <- list(
    speciesKey = key,
    limit = 99999,
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE,
    eventDate = period
  )

  if (!is.null(occurrence_status)) {
    data_args$occurrenceStatus <- occurrence_status
  }

  if (!is.null(publishing_org)) {
    data_args$publishingOrg <- publishing_org
  }

  do.call(rgbif::occ_data, data_args)$data
}

#' @description Function to get GBIF occurrence data for a given species
#' @param taxon_key A data.table returned from `batch_upd_dates()`.
#' @param path Character string to the Data folder.
#' @return A list containing raw data, metadata, and cleaned occurrence data.
#' @export
get_gbif_pts <- function(taxon_key, path) {
  tmp_dir <- file.path(path, "Original", "host-occurrences", "tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  taxon_key <- taxon_key[taxon_key$database == "gbif", ]

  key <- taxon_key$speciesKey
  sp <- taxon_key$species
  dateretrieved <- taxon_key$dateretrieved
  db <- taxon_key$database
  date_range <- taxon_key$date_range

  start_date <- as.Date(strsplit(date_range, ",")[[1]][1])
  end_date <- as.Date(strsplit(date_range, ",")[[1]][2])

  tmp_file <- make_tmp_file(tmp_dir, sp, date_range, db)

  cached <- read_cached_occurrences(
    tmp_file = tmp_file,
    tmp_dir = tmp_dir,
    sp = sp,
    start_date = start_date,
    end_date = end_date,
    db = db,
    progress = progress
  )

  if (!is.null(cached)) {
    return(cached)
  }

  day_diff <- if (is.na(dateretrieved)) {
    7
  } else {
    as.integer(Sys.Date() - as.Date(dateretrieved))
  }

  if (day_diff > 6) {
    dates <- get_dates(date_range, db = db)
    dates <- dates[dates <= end_date]

    inat <- "28eb1a3f-1c15-4a95-931a-4af90ecb574d"
    index <- seq_along(dates)[-length(dates)]

    occ_count_progress(
      db = db,
      sp = sp,
      key = key,
      date_range = date_range,
      progress = progress
    )

    presences0 <- lapply(index, function(x) {
      Sys.sleep(1)

      get_gbif_period(
        key = key,
        sp = sp,
        start_date = dates[x],
        end_date = dates[x + 1],
        occurrence_status = NULL,
        progress = progress
      )
    })

    presences0 <- rbind_occ_results(presences0)

    presences1 <- if (dates[1] != "1600-01-01") {
      Sys.sleep(1)

      rbind_occ_results(get_gbif_period(
        key = key,
        sp = sp,
        start_date = "1600-01-01",
        end_date = "2007-12-31",
        occurrence_status = NULL,
        publishing_org = inat,
        progress = progress
      ))
    } else {
      NULL
    }

    absences <- rbind_occ_results(get_gbif_period(
      key = key,
      sp = sp,
      start_date = dates[1],
      end_date = dates[length(dates)],
      occurrence_status = "ABSENT",
      progress = progress
    ))

    occs <- data.table::rbindlist(
      list(presences0, absences, presences1),
      fill = TRUE
    )

    if (nrow(occs) > 0 && "key" %in% names(occs)) {
      occs <- occs[!duplicated(occs$key), ]
    }

    if (nrow(occs) > 0) {
      occs$eventDate2 <- lubridate::as_date(occs$eventDate)

      failed <- occs[is.na(occs$eventDate2), ]
      failed$eventDate2 <- lubridate::as_date(paste(
        failed$year,
        failed$month,
        failed$day,
        sep = "-"
      ))

      occs <- rbind(occs, failed)
      occs <- occs[!is.na(occs$eventDate2), ]
      occs$eventDate2 <- lubridate::ymd(occs$eventDate2)

      inat_time <- lubridate::ymd("2008-01-01")

      occs <- occs[
        occs$eventDate2 < inat_time |
          (occs$publishingOrgKey != inat & occs$eventDate2 > inat_time),
      ]
    }

    occs <- list(occs)

    if (nrow(occs[[1]]) == 0) {
      occs[[2]] <- data.table::data.table(
        total = 0,
        present = 0,
        absent = 0,
        speciesKey = key,
        species = sp,
        database = db,
        dateretrieved = Sys.Date(),
        start_date = start_date,
        end_date = end_date,
        date_range = date_range
      )
    } else {
      occs[[2]] <- unique(data.table::data.table(
        total = nrow(occs[[1]]),
        present = sum(occs[[1]]$occurrenceStatus == "PRESENT"),
        absent = sum(occs[[1]]$occurrenceStatus == "ABSENT"),
        speciesKey = occs[[1]]$speciesKey,
        species = occs[[1]]$species,
        dateretrieved = Sys.Date(),
        start_date = start_date,
        end_date = end_date,
        date_range = date_range,
        database = db
      ))

      occs[[3]] <- unique(data.table::data.table(
        fkey = occs[[1]]$gbifID,
        date = occs[[1]]$eventDate2,
        p_a = occs[[1]]$occurrenceStatus,
        lat = occs[[1]]$decimalLatitude,
        lon = occs[[1]]$decimalLongitude,
        db = db,
        sciname = occs[[1]]$species
      ))

      occs[[3]]$p_a[occs[[3]]$p_a == "PRESENT"] <- 1
      occs[[3]]$p_a[occs[[3]]$p_a == "ABSENT"] <- 0
      occs[[3]]$date <- lubridate::year(occs[[3]]$date)
    }

    occs <- as_character_list(occs)
  } else {
    occs <- empty_occ_list()
  }

  finish_occurrence_output(
    occs = occs,
    taxon_key = taxon_key,
    path = path,
    tmp_file = tmp_file
  )
}

#' @description Function to get iNat occurrence data for a given species
#' @param taxon_key A data.table returned from `batch_upd_dates()`.
#' @param path Character string to the Data folder.
#' @return A list containing raw data, metadata, and cleaned occurrence data.
#' @export
get_inat_pts <- function(taxon_key, path) {
  tmp_dir <- file.path(path, "Original", "host-occurrences", "tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  taxon_key <- taxon_key[taxon_key$database == "inat", ]

  key <- taxon_key$speciesKey
  sp <- taxon_key$species
  db <- "inat"
  dateretrieved <- unique(taxon_key$dateretrieved)
  date_range <- unique(taxon_key$date_range)

  start_date <- as.Date(strsplit(date_range, ",")[[1]][1])
  end_date <- as.Date(strsplit(date_range, ",")[[1]][2])

  tmp_file <- make_tmp_file(tmp_dir, sp, date_range, db)

  cached <- read_cached_occurrences(
    tmp_file = tmp_file,
    tmp_dir = tmp_dir,
    sp = sp,
    start_date = start_date,
    end_date = end_date,
    db = db,
    progress = progress
  )

  if (!is.null(cached)) {
    return(cached)
  }

  day_diff <- if (is.na(dateretrieved)) {
    7
  } else {
    as.numeric(Sys.Date() - as.Date(dateretrieved))
  }

  if (day_diff > 6) {
    dates <- get_dates(date_range, db = db)

    occ_count_progress(
      db = db,
      sp = sp,
      key = key,
      date_range = date_range,
      progress = progress
    )

    if (is.list(dates)) {
      occs <- list(
        get_inat_period(sp, dates$years, "year", date_range, progress),
        get_inat_period(sp, dates$months, "month", date_range, progress),
        get_inat_period(sp, dates$days, "day", date_range, progress)
      )
    } else {
      occs <- get_inat_period(sp, dates, "day", date_range, progress)
    }

    occs <- rbind_occ_results(occs)

    if ("coordinates_obscured" %in% names(occs)) {
      occs <- occs[occs$coordinates_obscured == "false"]
    }

    occs <- list(occs)

    if (nrow(occs[[1]]) == 0) {
      occs[[2]] <- data.table::data.table(
        total = 0,
        present = 0,
        absent = 0,
        speciesKey = key,
        species = sp,
        database = db,
        date_range = date_range,
        start_date = as.character(start_date),
        end_date = as.character(end_date),
        dateretrieved = Sys.Date()
      )
    } else {
      occs[[1]]$species <- sp

      occs[[2]] <- unique(data.table::data.table(
        total = nrow(occs[[1]]),
        present = nrow(occs[[1]]),
        absent = 0,
        speciesKey = key,
        species = occs[[1]]$species,
        database = db,
        date_range = date_range,
        start_date = as.character(start_date),
        end_date = as.character(end_date),
        dateretrieved = Sys.Date()
      ))

      if (
        "coordinates_obscured" %in% names(occs[[1]]) &&
          "taxon_geoprivacy" %in% names(occs[[1]])
      ) {
        occs[[1]] <- unique(occs[[1]][
          occs[[1]]$coordinates_obscured != "true" &
            occs[[1]]$taxon_geoprivacy != "obscured"
        ])
      }

      occs[[3]] <- unique(data.table::data.table(
        fkey = occs[[1]]$id,
        date = occs[[1]]$observed_on,
        p_a = 1,
        lat = occs[[1]]$latitude,
        lon = occs[[1]]$longitude,
        db = db,
        sciname = occs[[1]]$species
      ))

      occs[[3]]$date <- lubridate::ymd(occs[[3]]$date)
      occs[[3]]$date <- lubridate::year(occs[[3]]$date)
    }

    occs <- as_character_list(occs)
  } else {
    occs <- empty_occ_list()
  }

  finish_occurrence_output(
    occs = occs,
    taxon_key = taxon_key,
    path = path,
    tmp_file = tmp_file
  )
}

#' @description Function to get records from the BIEN database
#' @param taxon_key A data.table returned from `batch_upd_dates()`.
#' @param path Character string to the Data folder.
#' @return A list containing raw data, metadata, and cleaned occurrence data.
#' @export
get_bien_pts <- function(taxon_key, path) {
  tmp_dir <- file.path(path, "Original", "host-occurrences", "tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  taxon_key <- taxon_key[taxon_key$database == "bien", ]

  key <- taxon_key$speciesKey
  sp <- taxon_key$species
  db <- "bien"
  date_range <- taxon_key$date_range
  dateretrieved <- taxon_key$dateretrieved

  start_date <- as.Date(strsplit(date_range, ",")[[1]][1])
  end_date <- as.Date(strsplit(date_range, ",")[[1]][2])

  tmp_file <- make_tmp_file(tmp_dir, sp, date_range, db)

  cached <- read_cached_occurrences(
    tmp_file = tmp_file,
    tmp_dir = tmp_dir,
    sp = sp,
    start_date = start_date,
    end_date = end_date,
    db = db,
    progress = progress
  )

  if (!is.null(cached)) {
    return(cached)
  }

  last_updated <- BIEN::BIEN_metadata_database_version()$db_release_date
  last_updated <- as.Date(last_updated)

  if (is.na(dateretrieved)) {
    dateretrieved <- as.Date(last_updated) - 1
  }

  if (last_updated > dateretrieved) {
    occs <- BIEN::BIEN_occurrence_species(
      species = sp,
      natives.only = FALSE,
      observation.type = TRUE,
      collection.info = TRUE,
      only.geovalid = TRUE,
      cultivated = TRUE
    )

    occs <- data.table::as.data.table(occs)
    names(occs)[match("date_collected", names(occs))] <- "date_collected2"
    occs <- occs[occs$datasource != "GBIF" & occs$datasource != "FIA"]
    occs <- list(occs)

    if (nrow(occs[[1]]) == 0) {
      occs[[2]] <- data.table::data.table(
        total = 0,
        present = 0,
        absent = 0,
        speciesKey = key,
        species = sp,
        database = db,
        date_range = paste0(start_date[1], ",", end_date[1]),
        start_date = start_date,
        end_date = end_date,
        dateretrieved = Sys.Date()
      )
    } else {
      occs[[1]]$species <- sp

      occs[[2]] <- unique(data.table::data.table(
        total = nrow(occs[[1]]),
        present = nrow(occs[[1]]),
        absent = 0,
        speciesKey = key,
        species = occs[[1]]$scrubbed_species_binomial,
        database = db,
        date_range = paste0(start_date[1], ",", end_date[1]),
        start_date = start_date,
        end_date = end_date,
        dateretrieved = Sys.Date()
      ))

      id <- lubridate::ymd_hms(Sys.time())
      id <- paste0(
        lubridate::year(id), "_",
        lubridate::month(id), "_",
        lubridate::month(id), "_",
        lubridate::hour(id), "_",
        lubridate::minute(id)
      )
      id <- paste0(seq_len(nrow(occs[[1]])), "_", id)

      date_col <- if ("date_collected" %in% names(occs[[1]])) {
        "date_collected"
      } else {
        "date_collected2"
      }

      occs[[3]] <- unique(data.table::data.table(
        fkey = id,
        date = occs[[1]][[date_col]],
        p_a = 1,
        lat = occs[[1]]$latitude,
        lon = occs[[1]]$longitude,
        db = db,
        sciname = occs[[1]]$species
      ))

      occs[[3]]$date <- lubridate::ymd(occs[[3]]$date)
      occs[[3]]$date <- lubridate::year(occs[[3]]$date)
    }

    occs <- as_character_list(occs)
  } else {
    occs <- empty_occ_list()
  }

  finish_occurrence_output(
    occs = occs,
    taxon_key = taxon_key,
    path = path,
    tmp_file = tmp_file
  )
}

#' @description Function to get records from the NEON database
#' @param taxon_key A data.table returned from `batch_upd_dates()`.
#' @param path Character string to the Data folder.
#' @param token NEON API token. Defaults to `NULL`.
#' @param max_cores Maximum number of cores to use. Defaults to `1`.
#' @return A list containing raw data, metadata, and cleaned occurrence data.
#' @export
get_neon_pts <- function(taxon_key, path, token = NULL, max_cores = 1) {
  tmp_dir <- file.path(path, "Original", "host-occurrences", "tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  taxon_key <- taxon_key[taxon_key$database == "neon", ]

  key <- taxon_key$speciesKey
  sp <- taxon_key$species
  db <- "neon"
  date_range <- taxon_key$date_range
  dateretrieved <- taxon_key$dateretrieved

  start_date <- as.Date(strsplit(date_range, ",")[[1]][1])
  end_date <- as.Date(strsplit(date_range, ",")[[1]][2])

  tmp_file <- make_tmp_file(tmp_dir, sp, date_range, db)

  cached <- read_cached_occurrences(
    tmp_file = tmp_file,
    tmp_dir = tmp_dir,
    sp = sp,
    start_date = start_date,
    end_date = end_date,
    db = db,
    progress = progress
  )

  if (!is.null(cached)) {
    return(cached)
  }

  path2 <- paste0(path, "Original/host-occurrences/plot_db/neon/")
  dir.create(path2, recursive = TRUE, showWarnings = FALSE)

  neon_raw <- paste0(path2, "neon_raw.csv")
  neon_exists <- file.exists(neon_raw)

  if (neon_exists && !is.na(dateretrieved)) {
    start_date <- as.Date(file.info(neon_raw)$mtime)
  }

  progress(paste0(
    "Downloading & stacking NEON data for ", sp, " from ",
    start_date, " to ", end_date
  ))

  files <- neonUtilities::queryFiles(
    dpID = "DP1.10058.001",
    startdate = format(as.Date(start_date), "%Y-%m"),
    enddate = format(as.Date(end_date), "%Y-%m"),
    token = if (!is.null(token)) token else NA_character_
  )$files

  stackfiles <- list.files(
    path2,
    recursive = TRUE,
    full.names = TRUE
  )

  stackfiles <- stackfiles[!grepl("zip", stackfiles)]

  if (!all(basename(files) %in% basename(stackfiles))) {
    progress("Downloading missing NEON files...")

    neonUtilities::zipsByProduct(
      dpID = "DP1.10058.001",
      check.size = FALSE,
      savepath = path2,
      startdate = format(as.Date(start_date), "%Y-%m"),
      enddate = format(as.Date(end_date), "%Y-%m"),
      token = if (!is.null(token)) token else NA_character_
    )

    neonUtilities::stackByTable(
      filepath = paste0(path2, "filesToStack10058"),
      savepath = paste0(path2, "filesToStack10058"),
      folder = TRUE,
      saveUnzippedFiles = TRUE,
      nCores = pmin(parallel::detectCores() - 1, max_cores),
      dpID = "DP1.10058.001"
    )
  }

  div1_files <- list.files(
    path2,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "div_1"
  )

  occs <- lapply(div1_files, function(x) {
    data.table::fread(x, colClasses = "character")
  })

  occs <- data.table::rbindlist(occs, fill = TRUE)
  occs <- list(occs)

  matched <- occs[[1]][
    grepl(sp, occs[[1]]$scientificName, ignore.case = TRUE)
  ]

  occs[[2]] <- unique(data.table::data.table(
    uid = matched$uid,
    plotID = matched$plotID,
    decimalLatitude = matched$decimalLatitude,
    decimalLongitude = matched$decimalLongitude,
    endDate = matched$endDate,
    eventID = matched$eventID,
    sciname = sp
  ))

  all_plots <- unique(data.table::data.table(
    plotID = occs[[1]]$plotID,
    decimalLatitude = occs[[1]]$decimalLatitude,
    decimalLongitude = occs[[1]]$decimalLongitude,
    endDate = occs[[1]]$endDate,
    eventID = occs[[1]]$eventID
  ))

  occs[[2]] <- merge(
    occs[[2]],
    all_plots,
    by = c(
      "plotID", "decimalLatitude", "decimalLongitude",
      "endDate", "eventID"
    ),
    all = TRUE
  )

  occs[[2]]$p_a <- ifelse(is.na(occs[[2]]$sciname), 0, 1)

  occs[[2]]$uid <- ifelse(
    is.na(occs[[2]]$sciname),
    paste0(occs[[2]]$plotID, "__", occs[[2]]$eventID),
    occs[[2]]$uid
  )

  occs[[2]]$sciname <- sp

  occs[[2]] <- data.table::data.table(
    fkey = occs[[2]]$uid,
    date = occs[[2]]$endDate,
    p_a = occs[[2]]$p_a,
    lat = occs[[2]]$decimalLatitude,
    lon = occs[[2]]$decimalLongitude,
    db = db,
    sciname = occs[[2]]$sciname
  )

  occs[[3]] <- occs[[2]]

  occs[[2]] <- data.table::data.table(
    total = nrow(occs[[2]]),
    present = sum(occs[[2]]$p_a == 1),
    absent = sum(occs[[2]]$p_a == 0),
    speciesKey = key,
    species = sp,
    database = db,
    date_range = paste0(start_date[1], ",", end_date[1]),
    start_date = start_date,
    end_date = end_date,
    dateretrieved = Sys.Date()
  )

  occs[[3]]$date <- lubridate::ymd(occs[[3]]$date)
  occs[[3]]$date <- lubridate::year(occs[[3]]$date)

  occs <- as_character_list(occs)

  finish_occurrence_output(
    occs = occs,
    taxon_key = taxon_key,
    path = path,
    tmp_file = tmp_file
  )
}

#' @description Export occurrence outputs to disk and update metadata
#' @param occs An occurrence list returned by a database download function.
#' @param path Character string to the Data folder.
#' @return Invisibly returns `TRUE`.
#' @export
export_occurrences <- function(occs, path) {
  for (i in seq_along(occs)) {
    if (nrow(occs[[i]]) == 0) {
      next
    }

    output_file <- paste0(path, names(occs)[i])

    dir.create(
      dirname(output_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    new_data <- data.table::as.data.table(occs[[i]])

    if (file.exists(output_file)) {
      old_data <- data.table::fread(
        output_file,
        colClasses = "character"
      )

      if (basename(output_file) == "hostpts_metadata.csv") {
        old_key <- paste(old_data$speciesKey, old_data$database)
        new_key <- paste(new_data$speciesKey, new_data$database)
        old_data <- old_data[!old_key %in% new_key]
      }

      out <- data.table::rbindlist(
        list(old_data, new_data),
        fill = TRUE
      )
    } else {
      out <- new_data
    }

    data.table::fwrite(unique(out), output_file)
  }

  invisible(TRUE)
}

#' @description Export occurrence outputs returned across multiple databases
#' @param occs_by_db A named list returned by `get_species_pts()`.
#' @param path Character string to the Data folder.
#' @return Invisibly returns `TRUE`.
#' @export
export_species_occurrences <- function(occs_by_db, path) {
  lapply(occs_by_db, export_occurrences, path = path)

  invisible(TRUE)
}

#' @description Download occurrence data for one species from one database
#' @param scientific_name Character string giving the species name.
#' @param taxon_rank Character string giving the taxon rank.
#' @param path Character string to the Data folder.
#' @param db Character string identifying the database.
#' @param token Optional NEON API token. Used only when `db = "neon"`.
#' @param max_cores Maximum number of cores to use for NEON downloads.
#' @param export Logical. If TRUE, export outputs to disk.
#' @return A list containing raw data, metadata, & cleaned occurrence data.
#' @export
get_species_db_pts <- function(scientific_name, taxon_rank = "SPECIES", path,
                               db = c("gbif", "inat", "bien", "neon"),
                               token = NULL, max_cores = 1,
                               export = FALSE) {
  db <- match.arg(db)

  taxon_key <- batch_upd_dates(
    get_gbif_keys(scientific_name, taxon_rank),
    path
  )

  out <- switch(db,
    gbif = get_gbif_pts(taxon_key, path),
    inat = get_inat_pts(taxon_key, path),
    bien = get_bien_pts(taxon_key, path),
    neon = get_neon_pts(
      taxon_key = taxon_key,
      path = path,
      token = token,
      max_cores = max_cores
    )
  )

  if (isTRUE(export)) {
    export_occurrences(out, path)
  }

  out
}

#' @description Download occurrence data for one species across databases
#' @param scientific_name Character string giving the species name.
#' @param taxon_rank Character string giving the taxon rank.
#' @param path Character string to the Data folder.
#' @param db Character vector of databases to query. Defaults to all databases.
#' @param token Optional NEON API token. Used only when `"neon"` is included.
#' @param max_cores Maximum number of cores to use for NEON downloads.
#' @param export Logical. If TRUE, export outputs to disk.
#' @return A named list of occurrence outputs, one element per database.
#' @export
get_species_pts <- function(scientific_name, taxon_rank = "SPECIES", path,
                            db = c("gbif", "inat", "bien", "neon"),
                            token = NULL, max_cores = 1,
                            export = FALSE) {
  db <- match.arg(db, several.ok = TRUE)

  taxon_key <- batch_upd_dates(
    get_gbif_keys(scientific_name, taxon_rank),
    path
  )

  out <- lapply(db, function(x) {
    switch(x,
      gbif = get_gbif_pts(taxon_key, path),
      inat = get_inat_pts(taxon_key, path),
      bien = get_bien_pts(taxon_key, path),
      neon = get_neon_pts(
        taxon_key = taxon_key,
        path = path,
        token = token,
        max_cores = max_cores
      )
    )
  })

  names(out) <- db

  if (isTRUE(export)) {
    export_species_occurrences(out, path)
  }

  out
}

#' Example usage of `get_species_pts()`
#' @examples
species <- c(
  "Ailanthus altissima",
  "Acer rubrum",
  "Juglans nigra",
  "Vitis vinifera",
  "Vitis labrusca"
)

for (sp in species) {
  start_time <- Sys.time()
  get_species_pts(
    scientific_name = sp,
    taxon_rank = "SPECIES",
    path = "C:/Users/blaginh/Desktop/Data/",
    db = c("gbif", "inat", "bien", "neon"),
    export = TRUE,
    token = gsub(
      "NEON: ", "",
      readLines("C:/Users/blaginh/Desktop/neon_token.txt")
    ),
    max_cores = 4
  )
  total_time <- Sys.time() - start_time
  time_units <- attr(total_time, "units")
  progress(sprintf("Finished %s in %.2f %s", sp, total_time, time_units))
}
