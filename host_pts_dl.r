# Load libraries
require(data.table)
require(lubridate)
require(rgbif)
require(rinat)
require(BIEN)
require(neonUtilities)


#' @description  Generic function to print messages
#' @param message character string of message to be printed
#' @return prints message to console
#' @export progress
#' @examples progress("Getting GBIF presences for Ailanthus altissima from
#' 2010-01-01 to 2019-12-31")

progress <- function(message) {
  print(paste0(message))
}

#' @import rgbif
#' @import data.table
#' @param scientific_name character string of taxon name
#' @param taxon_rank character string of taxon rank, can only be "GENUS" or
#' "SPECIES"
#' @return a data.table of the taxon key for the given taxon name and taxon rank
#' @export get_gbif_keys
#' @examples get_gbif_keys("Ailanthus altissima", "SPECIES")

get_gbif_keys <- function(scientific_name, taxon_rank) {
  # convert taxon_rank to uppercase
  taxon_rank <- toupper(taxon_rank)
  # check that taxon_rank is valid
  if (taxon_rank == "GENUS") {
    # get the usage key for the genus
    keys <- rgbif::name_backbone(name = scientific_name, rank = taxon_rank,
                                 strict = TRUE)$usageKey
    # use genus usage key to find all species within the genus
    species <- rgbif::name_lookup(higherTaxonKey = keys, status = "ACCEPTED",
                                  rank = "SPECIES", limit = 99999)$data
    species <- data.table::as.data.table(species)
    species <- (species[, .(species, speciesKey, genusKey)]) # nolint
  } else if (taxon_rank == "SPECIES") {
    # get the usage key for the species
    keys <- rgbif::name_backbone(name = scientific_name, rank = taxon_rank,
                                 strict = TRUE)
    keys <- data.table::as.data.table(keys)
    species <- keys[, .(species, speciesKey, genusKey)] # nolint
  } else {
    stop("confirm that taxon rank is 'GENUS' or 'SPECIES'")
  }
  return(species)
}

#' @description Function to add date range to taxon key for each species based
#' on the existing metadata file. If the taxon key is not in the metadata file,
#' the date range will be set to 1600-01-01 to the current date.
#' @param path character string to the Data folder
#' @param taxon_key data.table of taxon key and species name from get_gbif_keys
#' @return a data.table of the taxon key with the date range added
#' @export upd_dates
#' @examples
#' folder <- paste0(system.file("extdata",package="gbif4SDM"),"/")
#' taxon_key <- data.table::fread(paste0(folder,"hostpts_metadata.csv"))
#' upd_dates(folder, key = taxon_key$speciesKey, sp = taxon_key$species,
#' db = taxon_key$database)

batch_upd_dates <- function(path, taxon_key) {
  sp <- taxon_key$species
  key <- taxon_key$speciesKey
  # Read in the existing metadata file as data.table
  metadata <- data.table::fread(paste0(path,
                                       "Original/host-occurrences/hostpts_metadata.csv"), # nolint
                                colClasses = "character")

  # Convert start_date and end_date to date format
  metadata$start_date <- lubridate::ymd(metadata$start_date)
  metadata$end_date <- lubridate::ymd(metadata$end_date)

  # Pull out the date_range for each species based on the existing metadata
  db_values <- c("gbif", "inat", "bien", "neon")
  result <- lapply(db_values, function(db) {
    # Filter by key
    metadata <- metadata[metadata$speciesKey == key & metadata$database == db, ]
    dateretrieved <- max(metadata$dateretrieved)

    # If metadata is not empty, take max end_date
    if(nrow(metadata) > 0) {
      date_range <- paste0(max(metadata$end_date) + 1, ",", Sys.Date())
    } else {
      # If database is gbif date_range is set to 1600-01-01 to current date
      if (db == "gbif") {
        date_range <- paste0("1600-01-01", ",", Sys.Date())
        # Else if database is inat date_range is set to 2008-01-01 to Sys.Date()
      } else if (db == "inat") {
        date_range <- paste0("2008-01-01", ",", Sys.Date())
      } else if (db == "bien") {
        date_range <- paste0("NA", ",", Sys.Date())
      } else if (db == "neon") {
        date_range <- paste0("2013-05-31", ",", Sys.Date())
      } else {
        stop("confirm that database is 'gbif', 'inat', 'neon', 'bien'")
      }
    }

    data.table::as.data.table(cbind(speciesKey = key, species = sp,
                                    date_range, database = db,
                                    dateretrieved = dateretrieved))
  })
  taxon_key_upd <- data.table::rbindlist(result, fill = TRUE)
  return(taxon_key_upd)
}



#' @description Function to create a vector of dates to avoid API limits
#' @param start_date character string of start date
#' @param end_date character string of end date
#' @param db character string of the database
#' @return a vector of dates
#' @export get_dates

get_dates <- function(date_range, db) {
  start_date <- as.Date(strsplit(date_range, ",")[[1]][1])
  end_date <- as.Date(strsplit(date_range, ",")[[1]][2])

  if (start_date == "1600-01-01") {
    dates <- c(start_date, as.Date("1960-12-31"))
    dates <- c(dates, seq((dates)[length(dates)], as.Date(Sys.Date()),
                          by = "1 year")[-1])
    dates <- c(dates, Sys.Date())
  } else if (start_date == "2008-01-01") {
    dates <- list(years = seq(as.Date("2008-01-01"), as.Date("2018-01-01"),
                              by = "year"), months = NULL, days = NULL)
    dates[[2]] <- seq(dates[[1]][length(dates[[1]])], Sys.Date(), by = "month")
    dates[[3]] <- seq(dates[[2]][length(dates[[2]])], Sys.Date(), by = "day")
    dates[[2]] <- dates[[2]][-length(dates[[2]])]
    names(dates)
  } else if (db == "gbif") {
    dates <- seq(start_date+1, end_date, by = "month")
    dates <- sort(c(dates, end_date))

  } else {
    dates <- seq(start_date, end_date, by = "day")
  }
  return(dates)
}

#' @description Function to create directories
#' @param path character string to the Data folder
#' @param taxon_key data.table of taxon key and species name from get_gbif_keys
#' @param occs list of occurrence data from get_gbif_pts, get_inat_pts,
#' get_bien_pts, get_neon_pts
#' @return list of directories

append_dirs <- function(occs, taxon_key, path) {
  sp <- taxon_key$species
  db <- taxon_key$database
  species_fn <- gsub(" ", "_", sp)
  if (db == "neon") {
    filenames <- c(
      paste0(
        "Original/host-occurrences/plot_db/neon/neon_raw.csv"
      ),
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
  lapply(seq_along(filenames), function(x) {
    dir.create(paste0(path, dirname(filenames[x])),
      recursive = TRUE,
      showWarnings = TRUE
    )
  })
  names(occs) <- filenames
  return(occs)
}

#' @title Download GBIF occurrence data for a given species
#' @description Function to get the number of records for each species
#' @param taxon_key data.table returned from batch_upd_dates
#' @param path character string to the Data folder
#' @return a data.table of the taxon key with the number of records added
#' @export get_inat_counts

get_gbif_pts <- function(taxon_key, path) {
  # Retain gbif metadata
  taxon_key <- taxon_key[taxon_key$database == "gbif", ]
  
  # Setup sp name, key, and date_range
  key <- taxon_key$speciesKey
  sp <- taxon_key$species
  date_range <- taxon_key$date_range
  dates <- get_dates(date_range, db = "gbif")

  # Setup gbif function params
  inat <- "28eb1a3f-1c15-4a95-931a-4af90ecb574d"
  index <- seq_along(dates)[-length(dates)]

  # Retrieve gbif presences
  presences0 <- lapply(index, function(x) {
    Sys.sleep(1) # pause for 1 seconds to avoid overloading the server
    data <- rgbif::occ_data(
      speciesKey = key, limit = 99999,
      hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
      eventDate = paste0(dates[x], ",", dates[x + 1])
    )$data
    progress(paste0("Getting GBIF presences for ", sp, " from ", dates[x],
                    " to ", dates[x + 1]))
    return(data)
  })

  # Retrieve historical inat presences (periodically added to gbif)
  if (dates[1] != "1600-01-01") {
    presences1 <- rgbif::occ_data(
      speciesKey = key, limit = 99999,
      hasCoordinate = TRUE, hasGeospatialIssue = FALSE, eventDate =
        paste0("1600-01-01,2007-12-31"), publishingOrg = inat
    )$data
  } else {
    presences1 <- NULL
  }

  # Retrieve gbif absences
  absences <- rgbif::occ_data(
    speciesKey = key,
    limit = 99999, hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
    eventDate = paste0(dates[1], ",", dates[length(dates)]),
    occurrenceStatus = "ABSENT"
  )$data

  # Remove NULL values from presences and absences
  presences0 <- presences0[!sapply(presences0, is.null)]
  # Rbind using data.table::rbindlist
  presences0 <- data.table::rbindlist(presences0, fill = TRUE)
  occs <- data.table::rbindlist(list(presences0, absences, presences1),
    fill = TRUE
  )

  # Convert date format to 'yyyy-mm-dd'
  occs$eventDate2 <- lubridate::as_date(occs$eventDate)
  failed <- occs[is.na(occs$eventDate2), ] # correct failed dates
  failed$eventDate2 <- lubridate::as_date(paste(failed$year, failed$month,
    failed$day,
    sep = "-"
  ))
  occs <- rbind(occs, failed)
  occs <- occs[!is.na(occs$eventDate2), ]
  occs$eventDate2 <- lubridate::ymd(occs$eventDate2)

  # Remove any inat records added to gbif after the inat release date
  inat_time <- lubridate::ymd("2008-01-01")
  occs <- list(occs[occs$eventDate2 < inat_time |
                      (occs$publishingOrgKey != inat &
                        occs$eventDate2 > inat_time), ])

  if (nrow(occs[[1]]) == 0) {
    occs[[2]] <- data.table::as.data.table(cbind(
      total = 0, present = 0,
      absent = 0, speciesKey = key, species = sp, database =
        "gbif", date_range = paste0(dates[1], ",", dates[length(dates)]),
      start_date = as.character(dates[1]), end_date =
        as.character(dates[length(dates)], dateretrieved = dates[length(dates)])
    ))

  } else {
    # Create metadata file
    occs[[1]]$date_range <- date_range[length(dates)]
    occs[[2]] <- unique(occs[[1]][, .(
      total = .N, present = # nolint
        sum(occurrenceStatus == "PRESENT"), absent = # nolint
        sum(occurrenceStatus == "ABSENT"), speciesKey, # nolint
      species, database = "gbif", date_range = paste0(
        dates[1], ",", # nolint
        dates[length(dates)]
      ), start_date = dates[1],
      end_date = dates[length(dates)], dateretrieved = dates[length(dates)]
    )])
    # Clean occs data for sdm
    occs[[3]] <- unique(occs[[1]][, .(
      fkey = gbifID, date = eventDate2, # nolint
      p_a = occurrenceStatus, lat = decimalLatitude, lon = # nolint
        decimalLongitude, db = "gbif", sciname = species
    )]) # nolint
    # Convert p_a 'PRESENT' to 1 and 'ABSENT' to 0
    occs[[3]][p_a == "PRESENT", p_a := 1] # nolint
    occs[[3]][p_a == "ABSENT", p_a := 0] # nolint
    # Convert date format to 'yyyy-mm-dd-hh-mm-ss' and 'yyyy'
    occs[[3]]$date <- lubridate::year(occs[[3]]$date)
  }
  occs <- lapply(seq_along(occs), function(x) {
    occs[[x]] <- occs[[x]][, names(occs[[x]]) := lapply(
      .SD, # nolint
      as.character
    )]
  })
  append_dirs(occs, taxon_key, path)
  return(occs)
}


species <- c("Ailanthus altissima","Acer rubrum", "Juglans nigra", "Vitis vinifera", "Vitis labrusca")
path <- "Z:/pops_pesthostuse/pops.sdm/Data/"

lapply(species, function(x){
  taxon_key <- get_gbif_keys(x, "SPECIES")
  taxon_key <- batch_upd_dates(path = path, taxon_key = taxon_key)
  start_time <- Sys.time()
  occs <- get_gbif_pts(taxon_key, path)
  export_data(occs,path)
  end_time <- Sys.time()
  end_time - start_time
})

#' @description Generic function that will export/append raw and cleaned
#' occurrence data and metadata for a given species and database
#' @param occs is a list of data.tables containing the raw and cleaned
#' occurrence records return from get_gbif_pts, get_inat_pts, get_bien_pts,
#' get_blm_pts, get_fia_pts, get_neon_pts
#' @param path character string to the Data folder
#' @export export_data

export_data <- function(occs, path) {
  if (nrow(occs[[1]]) == 0) {
    progress(paste0(
      "No new ", occs[[2]]$database[1],
      " records found for this species."
    ))
    data <- data.table::fread(paste0(path, names(occs)[2]),
      colClasses = "character"
    )
    data <- unique(data.table::rbindlist(list(data, occs[[2]]),
      fill = TRUE
    ))
    data.table::fwrite(data, paste0(path, names(occs)[2]),
      row.names = FALSE
    )
  } else {
    x<-1
    # Use lapply to read in the existing data and append the new data
    lapply(seq_along(occs), function(x) {
      # Check if file exists
      if (file.exists(paste0(path, names(occs)[x]))) {
        # Read in existing data
        data <- data.table::fread(paste0(path, names(occs)[x]),
          colClasses = "character"
        )
        # Append new data
        data <- data.table::rbindlist(list(data, occs[[x]]), fill = TRUE)
        # Remove duplicate records
        data <- unique(data)
        # Export data
        data.table::fwrite(data, paste0(path, names(occs)[x]),
          row.names = FALSE
        )
      } else {
        # Export data
        data.table::fwrite(occs[[x]], paste0(
          path, names(occs)[x]
        ), row.names = FALSE)
      }
    })
  }
}

#' @description Function to get inat occurrence data for a given species
#' @param taxon_key of taxon key from batch_upd_dates
#' @param path character string to the Data folder
#' @return a data.table of the taxon key with the number of records added
#' @export get_inat_counts

get_inat_pts <- function(taxon_key, path) { # nolint
  # Retain inat metadata
  taxon_key <- taxon_key[taxon_key$database == "inat", ]
  key <- taxon_key$speciesKey
  sp <- taxon_key$species
  date_range <- taxon_key$date_range
  dates <- get_dates(date_range, db = "inat")

  # if dates is a list of dates
  if (is.list(dates)) {
    presences0 <- lapply(dates$years, function(x) {
      tryCatch(
        {
          Sys.sleep(3) # pause for 3 seconds to avoid overloading the server
          data <- rinat::get_inat_obs(
            taxon_name = sp, quality = "research",
            geo = TRUE, maxresults = 9999, meta = FALSE, year =
              lubridate::year(x)
          )
          progress(paste0("Getting iNat presences for ", sp, " from ",
                          lubridate::year(x)))
          return(data)
        },
        error = function(e) {
          return(NULL)
        }
      )
    })
    presences1 <- lapply(dates$months, function(x) {
      tryCatch(
        {
          Sys.sleep(3) # pause for 3 seconds to avoid overloading the server
          data <- rinat::get_inat_obs(
            taxon_name = sp, quality = "research",
            geo = TRUE, maxresults = 9999, meta = FALSE, year =
              lubridate::year(x), month = lubridate::month(x)
          )
          progress(paste0("Getting iNat presences for ", sp, " from ",
                          lubridate::year(x), "-", lubridate::month(x)))
          return(data)
        },
        error = function(e) {
          return(NULL)
        }
      )
    })
    presences2 <- lapply(dates$days, function(x) {
      tryCatch(
        {
          Sys.sleep(3) # pause for 3 seconds to avoid overloading the server
          data <- rinat::get_inat_obs(
            taxon_name = sp, quality = "research",
            geo = TRUE, maxresults = 9999, meta = FALSE, year =
              lubridate::year(x), month = lubridate::month(x), day =
              lubridate::day(x)
          )
          progress(paste0("Getting iNat presences for ", sp, " from ", x))
          return(data)

        },
        error = function(e) {
          return(NULL)
        }
      )
    })
    # Use lapply to rbindlist presence data
    occs <- list(presences0, presences1, presences2)
    occs <- lapply(seq_along(occs), function(x) {
      data.table::rbindlist(occs[[x]], fill = TRUE)
    })
  } else {
    occs <- lapply(dates, function(x) {
      tryCatch(
        {
          Sys.sleep(3) # pause for 3 seconds to avoid overloading the server
          data <- rinat::get_inat_obs(
            taxon_name = sp, quality = "research",
            geo = TRUE, maxresults = 9999, meta = FALSE, year =
              lubridate::year(x), month = lubridate::month(x), day =
              lubridate::day(x)
          )
          progress(paste0("Getting iNat presences for ", sp, " from ", x))
          return(data)
        },
        error = function(e) {
          return(NULL)
        }
      )
    })
  }

  occs <- unique(data.table::rbindlist(occs, fill = TRUE))
  occs <- list(occs)

  # split date_range into start_date and end_date
  dates <- strsplit(date_range, ",")[[1]]

  # Setup metadata file, raw data, and cleaned data
  if (nrow(occs[[1]]) == 0) {
    occs[[2]] <- data.table::as.data.table(cbind(
      total = 0, present = 0,
      absent = 0, speciesKey = key, species = sp, database =
        "inat", date_range = paste0(dates[1], ",", dates[2]),
      start_date = as.character(dates[1]), end_date =
        as.character(dates[2]), dateretrieved = as.character(dates[2])
    ))
  } else {
    # Create metadata file
    occs[[1]]$species <- sp
    occs[[2]] <- unique(occs[[1]][, .(
      total = .N, present = .N, # nolint
      absent = 0, speciesKey = key, species, database = "inat", # nolint
      date_range = paste0(dates[1], ",", dates[2]), start_date = dates[1],
      end_date = dates[2], dateretrieved = as.character(dates[2])
    )])
    # Clean occs data for sdm
    occs[[1]] <- unique(occs[[1]][coordinates_obscured != "true" & # nolint
      taxon_geoprivacy != "obscured"]) # nolint
    occs[[3]] <- unique(occs[[1]][, .(
      fkey = id, date = observed_on, # nolint
      p_a = 1, lat = latitude, lon = longitude, db = "inat", sciname = # nolint
        species
    )]) # nolint
    # Convert date format to 'yyyy-mm-dd-hh-mm-ss' and 'yyyy'
    occs[[3]]$date <- lubridate::ymd(occs[[3]]$date)
    occs[[3]]$date <- lubridate::year(occs[[3]]$date)
  }
  occs <- lapply(seq_along(occs), function(x) {
    occs[[x]] <- occs[[x]][, names(occs[[x]]) := lapply(
      .SD, # nolint
      as.character
    )]
  })
  occs <- append_dirs(occs, taxon_key, path)
  return(occs)
}

#' @description Function to get the number of records for each species
#' @param path character string to the Data folder
#' @param taxon_key data.table from batch_upd_dates
#' @return a data.table with the number of records added
#' @export get_bien_pts

taxon_key <- get_gbif_keys("Ailanthus altissima", "SPECIES")
taxon_key <- batch_upd_dates(path, taxon_key)

get_bien_pts <- function(path, taxon_key) {
  # Retain bien metadata
  taxon_key <- taxon_key[taxon_key$database == "bien", ]
  key <- taxon_key$speciesKey
  sp <- taxon_key$species
  date_range <- taxon_key$date_range

  if (grepl("NA", date_range)) {
    start_date <- "NA"
  } else {
    start_date <- as.Date(strsplit(date_range, ",")[[1]][1])
  }
  end_date <- as.Date(strsplit(date_range, ",")[[1]][2])

  occs <- BIEN::BIEN_occurrence_species(
    species = sp, natives.only = FALSE,
    observation.type = TRUE, collection.info = TRUE, only.geovalid = TRUE
  )

  occs <- data.table::as.data.table(occs)
  names(occs)[match("date_collected", names(occs))] <- "date_collected2"
  occs <- occs[datasource != "GBIF" & datasource != "FIA", ]
  occs <- list(occs)

  if (nrow(occs[[1]]) == 0) {
    occs[[2]] <- data.table::as.data.table(cbind(
      total = 0, present = 0,
      absent = 0, speciesKey = key, species = sp, database =
        "bien", date_range = paste0(start_date[1], ",", end_date[1]),
      start_date = start_date, end_date = end_date,
      dateretrieved = end_date
    ))
  } else {
    # Create metadata file
    occs[[1]]$species <- sp
    occs[[2]] <- unique(occs[[1]][, .(
      total = .N, present = .N, # nolint
      absent = 0, speciesKey = key, species = scrubbed_species_binomial, database = "bien", # nolint
      date_range = paste0(start_date[1], ",", end_date[1]), start_date = start_date,
      end_date = end_date,
      dateretrieved = end_date
    )])
    # Clean occs data for sdm
    id <- lubridate::ymd_hms(Sys.time())
    id <- paste0(
      lubridate::year(id), "_", lubridate::month(id), "_",
      lubridate::month(id), "_",
      lubridate::hour(id), "_", lubridate::minute(id)
    )
    id <- paste0(seq(1:nrow(occs[[1]])), "_", id)
    occs[[3]] <- unique(occs[[1]][, .(
      fkey = id, date = date_collected, # nolint
      p_a = 1, lat = latitude, lon = longitude, db = "bien", sciname = species
    )]) # nolint
    # Convert date format to 'yyyy-mm-dd-hh-mm-ss' and 'yyyy'
    occs[[3]]$date <- lubridate::ymd(occs[[3]]$date)
    occs[[3]]$date <- lubridate::year(occs[[3]]$date)
  }
  occs <- lapply(seq_along(occs), function(x) {
    occs[[x]] <- occs[[x]][, names(occs[[x]]) := lapply(.SD, as.character)]
  })
  occs <- append_dirs(occs, taxon_key, path)
  return(occs)
}

taxon_key <- get_gbif_keys("Ailanthus altissima", "SPECIES")
taxon_key <- batch_upd_dates(path, taxon_key)

#' @description Function to get the number of records for each species
#' @param taxon_key character string of the taxon key from get_gbif_keys
#' @param path character string to the Data folder
#' @return a data.table with the number of records added
#' @export get_neon_pts

get_neon_pts <- function(taxon_key, path) { # nolint
  # Retain neon metadata, key, species, and date_range
  taxon_key <- taxon_key[taxon_key$database == "neon", ]
  key <- taxon_key$speciesKey
  sp <- taxon_key$species
  date_range <- taxon_key$date_range

  # Obtain start_date and end_date from date_range
  start_date <- format(as.Date(strsplit(date_range, ",")[[1]][1]), "%Y-%m")
  end_date <- format(as.Date(strsplit(date_range, ",")[[1]][2]), "%Y-%m")
  
  # Create directory for neon data
  path2 <- paste0(path, "Original/host-occurrences/plot_db/neon/")
  dir.create(path2, recursive = TRUE, showWarnings = FALSE)

  # Download neon data using start_date and end_date
  neonUtilities::zipsByProduct(
    dpID = "DP1.10058.001", check.size = FALSE,
    savepath = path2, startdate = start_date,
    enddate = end_date
  )

  # Stack neon data
  neonUtilities::stackByTable(
    filepath = paste0(path2, "filesToStack10058"),
    savepath = paste0(path2,"filesToStack10058"),
    folder = TRUE, saveUnzippedFiles = FALSE, 
    nCores = parallel::detectCores() - 1,
    dpID = "DP1.10058.001"
  )

  # Pull in the plant data and combine
  div1_files <- list.files(path2, recursive = TRUE, full.names = TRUE,
                           pattern = "div_1")
  occs <- lapply(div1_files, function(x) {
    z <- data.table::fread(x)
    z$dateretrieved <- lubridate::ymd_hms(file.info(x)$mtime)
    z$dateretrieved <- lubridate::as_date(z$dateretrieved)
    z
  })

  occs <- data.table::rbindlist(occs, fill = TRUE)
  occs <- list(occs)
  occs[[2]] <- unique(occs[[1]][
    grepl(sp, occs[[1]]$scientificName,
      ignore.case = TRUE
    ),
    .(uid, plotID, decimalLatitude,
      decimalLongitude, endDate, eventID,
      sciname = sp
    )
  ])
  occs[[2]] <- merge(occs[[2]], unique(occs[[1]][, .(
    plotID, decimalLatitude,
    decimalLongitude, endDate,
    eventID
  )]),
  by = c(
    "plotID", "decimalLatitude", "decimalLongitude",
    "endDate", "eventID"
  ), all = TRUE
  )
  occs[[2]]$p_a <- ifelse(is.na(occs[[2]]$sciname), 0, 1)
  occs[[2]]$uid <- ifelse(is.na(occs[[2]]$sciname),
    paste0(occs[[2]]$plotID, "__", occs[[2]]$eventID),
    occs[[2]]$uid
  )
  occs[[2]]$sciname <- sp
  occs[[2]] <- occs[[2]][, .(
    fkey = uid, date = endDate, p_a = p_a,
    lat = decimalLatitude, lon = decimalLongitude,
    db = "neon", sciname
  )]
  occs[[3]] <- occs[[2]]
  # create metadata file
  start_date <- paste0(start_date, "-01")
  end_date <- paste0(end_date, "-01")
  occs[[2]] <- occs[[2]][, .(
    total = .N, present = sum(p_a == 1),
    absent = sum(p_a == 0), speciesKey = key,
    species = sp, database = "neon",
    date_range = paste0(
      start_date[1], ",",
      end_date[1]
    ),
    start_date = start_date, end_date = end_date,
    dateretrieved = unique(occs[[1]]$dateretrieved)
  )]
  occs[[3]]$date <- lubridate::ymd(occs[[3]]$date)
  occs[[3]]$date <- lubridate::year(occs[[3]]$date)
  occs[[3]]$dateretrieved <- unique(occs[[1]]$dateretrieved)
  occs <- lapply(seq_along(occs), function(x) {
    occs[[x]] <- occs[[x]][, names(occs[[x]]) := lapply(.SD, as.character)]
  })
  occs <- append_dirs(occs, taxon_key, path)
  return(occs)
}

```