# MAINTAIN A LOCAL REPOSITORY OF OISST
# AND COMPUTE GLOBAL SST WARMING TRENDS

# - requires "sstcoremulti", a C++ program to eficiently transform time slices (daily sst data) into pixel time-series (sst cores); "sstcoremulti" can be found in the same GitHub repository
# - must be run before "oisst_warming_trends_analysis.R", which will look for the "info" and "stat" files produced by this script

# DATASET:
# NOAA'S 1/4  arc-degree Daily
# 53 Optimum Interpolation SST version 2 (dOISST.v2)
# --
# Banzon, V., Smith, T. M., Chin, T. M., Liu, C., and Hankins, W. (2016). A long-term record of
# blended satellite and in situ sea-surface temperature for climate monitoring, modeling and
# environmental studies. Earth System Science Data 8, 165–176. doi:10.5194/essd-8-165-2016.
# --

# Local repository consists of all available daily oisst data since 1981-09-01
# 2 folders are created:
#   - "oisst_repository" is where the local repository is saved; inside...
#     - "daily" contains 1 file per day (basically a simplified version of the original data)
#     - "cores" contains 1 file per ocean pixel, with the entire time-series available
#   - "oisst_download_outputs" is where "info" and "stat" files are saved to allow easy reconstruction of the original map
# Repository data is saved in binary, with vals only for ocean pixels

# The "stat" files saved to "oisst_download_outputs" are:
#   - "slope_per_decade.RDS", the slope of the linear regression of the seasonally detrended sst time-series (aka warming rate, in degrees C per decade)
#   - "se_corr.RDS", the standard error of the slope, with a correction to the number of degrees of freedom in order to account for temporal autocorrelation

# NOTE: The code can be re-run to update the local oisst repository. In that case only oisst data not yet locally available will be downloaded. Rebuilding the cores (takes a long time) only happens after at least 90 new days of data have been added to the local repository, but  the user can force it by overriding the if statement. Updating of the "stat" files also only happens when cores are rebuilt.

libraries <- c("ncdf4", "raster", "doParallel", "stringr", "lubridate", "tibble", "dplyr", "readr", "purrr", "xts")
for (l in libraries) library(l, character.only = TRUE)

CORES <- 140
registerDoParallel(CORES)

## sstcoremulti ----
## this is a C++ app needed for efficient extraction of sst cores from the local daily oisst repository
has.sstcoremulti <- tryCatch(system("sstcoremulti", intern = TRUE, ignore.stderr = TRUE))
has.sstcoremulti <- !inherits(has.sstcoremulti, "try-error")
if(!has.sstcoremulti) {
  stop("sstcoremulti not available\nsst cores can't be built or maintained but oisst download can proceed\ncompile sstcoremulti to make sst cores")
}

## folders ----
url <- "https://www.ncei.noaa.gov/thredds/fileServer/OisstBase/NetCDF/AVHRR/XXXXX/"
rep <- "oisst_repository/"
if (!dir.exists(rep)) dir.create(rep)
daily_folder  <- str_c(rep, "daily")
cores_folder  <- str_c(rep, "cores")
if (!dir.exists(daily_folder)) dir.create(daily_folder)
if (!dir.exists(cores_folder)) dir.create(cores_folder)
fn_tmp     <- "tmp"
out_folder <- "oisst_download_outputs/"
if (!dir.exists(out_folder)) dir.create(cores_folder)
fn_info1   <- str_c(out_folder, "oisst_info1.RData")
fn_info2   <- str_c(out_folder, "oisst_info2.RData")


## date range ----
## # download data from t0 to t1 (or the closest available)
t0 <- as_date("1981-09-01")
t1 <- Sys.Date()

## lat range ----
# discard data beyond these latitudes
ylims <- c(-90, 90)

## file details ----
# a tibble with info about all sst files to be downloaded
x <- tibble(
  date1  = seq.Date(t0, t1, "day"),
  date2  = str_replace_all(date1, "-", ""),
  date3  = str_sub(date2, 1, 6),
  fnmWEB = str_c("avhrr-only-v2.", date2, ".nc"),
  fnmTMP = str_c(daily_folder, "/", date1, ".tmp"),
  fnmSST = str_c(daily_folder, "/", date1, "_sst"),
  url    = str_c(str_replace_all(url, "XXXXX", date3), fnmWEB))

## sst helper data ----
## contains data and functions needed for easy use of the oisst repository

# check if the fn_info1 file exists
# if yes, load it
# if not, prepare it now
if (!file.exists(fn_info1)) {
  cat("building fn_info1 file\n")
  download.file(x$url[1], fn_tmp, "curl")
  
  # download 1 oisst file to serve as template
  tmp <- nc_open(fn_tmp)
  tmp <- ncvar_get(tmp, varid = "sst", start = c(1,1,1,1), count = c(-1,-1,-1,-1))
  
  # make it a 0 (sea) NA (land) raster
  tmp <- tmp %>% t %>% raster %>% flip(direction = "y")
  extent(tmp) <- c(0, 360, -90, 90)
  projection(tmp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  names(tmp) <- "daily.sst"
  tmp@data@unit <- "degrees C"
  tmp[!is.na(tmp[])] <- 0
  sea_land_mask <- tmp
  unlink(fn_tmp)
  
  # make a vector with as many entries as pixels in the sea_land_mask raster
  # TRUE  means not NA, i.e., a sea pixel
  # FALSE means NA, i.e., a land pixel
  notNA <- !is.na(sea_land_mask[])
  
  # make a raster similar to the sea_land_mask, but now the value stored in sea pixels
  # is not 0, but instead their order (excluding land pixels)
  order_mask <- sea_land_mask
  i <- 1:sum(notNA)
  order_mask[notNA] <- i
  order_mask <- crop(order_mask, extent(c(0, 360, ylims)))
  
  # functions for data conversion
  ## convert numeric to binary
  ## add 30C so that there are no negative numbers and multiply by 100 so that when converted to integer we still retain 2 decimal places (as in the original data)
  rTObin  <- function(x) ((x + 30) * 100)
  ## convert binary back to numeric
  binTOr  <- function(x) ((x / 100) - 30)
  ## take a vector of values and rebuild the map
  ## values can be sst, stats or any other data, as long as they are numeric and with the same length the one required to fill all pixels that should have data
  rebuild <- function(vals, conv = NULL) {
    stopifnot(sum(notNA) == length(vals))
    if (!is.null(conv)) vals <- conv(vals)
    sea_land_mask[notNA] <- vals
    sea_land_mask
  }
  
  info1 <- tibble(
    item = c("sea_land_mask", "notNA", "order_mask", "rTObin", "binTOr", "rebuild"),
    info = c(
      "raster, 0 (sea) NA (land)",
      "vector, logical, TRUE (sea) FALSE (land)",
      "raster, pixel order (sea) NA (land)",
      "function, convert numeric to binary",
      "function, convert binary to numeric",
      "function, take vals and rebuild map"))
  
  save(list = c(info1$item, "info1"), file = fn_info1)
}
load(fn_info1)

## download files ----
## this is done sequentially so that it is easy to spot the last available file in the server
## only the first clean run of the script will take longer
## subsequent runs will only download files not yet in the database
old <- dir(daily_folder)
cat("downloading sst files\n")
void <- foreach(i = 1:nrow(x), .inorder = FALSE) %dopar% {
  fnmSST <- x$fnmSST[i]
  
  if (!file.exists(fnmSST)) {
    fnmURL <- x$url[i]
    fnmWEB <- x$fnmWEB[i]
    fnmTMP <- x$fnmTMP[i]
    date1  <- x$date1[i]
    
    download.file(fnmURL, fnmTMP, "curl", quiet = TRUE)
    if (file.info(fnmTMP)$size > 1000000) {
      tmp <- nc_open(fnmTMP)
      sst <- ncvar_get(tmp, varid = "sst", start = c(1,1,1,1), count = c(-1,-1,-1,-1))
      sst  <- sst %>% t %>% raster %>% flip(direction = "y")
      sst  <- values(sst)
      
      notNA2 <- !is.na(sst)
      if (!identical(notNA, notNA2)) stop()
      
      sst <- sst[notNA2]
      sst <- rTObin(sst)
      sst <- as.integer(round(sst))
      
      f <- file(fnmSST, "wb")
      writeBin(sst, f, size = 2)
      close(f)
    }
    unlink(fnmTMP)
  }
}
new <- dir(daily_folder)
new <- new[!(new %in% old)]
MSG1 <- if (length(new)) {
  str_c("sst files downloaded:\n   ", str_c(new, collapse = "\n   "))
}else{
  "local repository was already up-to-date\n"
}

if (has.sstcoremulti) {
  ## (re)build sst cores
  ## either make new cores from scratch
  ##
  ## or, if they are already present, update with new data
  ##   this is only performed automatically if there are more than 3 months
  ##   of new data to be appended to the existing sst cores
  
  sst_files  <- dir(daily_folder, pattern = "_sst", full.names = TRUE)
  core_files <- dir(cores_folder, pattern = "RDS",  full.names = TRUE)
  
  # compute the size of each core to closely
  # match the number of processor threads available
  datbin <- readBin(con = sst_files[1], what = "integer", n = 1440 * 720, size = 2)
  maxInd <- length(datbin)
  
  core_size <- (maxInd / CORES) / 1000
  core_size <- ceiling(core_size) * 1000
  
  MSG2 <- ""
  
  ## rebuild cores ----
  ## zero core files available = rebuilding from scratch
  if (!length(core_files)) {
    MSG2 <- "sst cores rebuilt from scratch\n"
    
    unlink(cores_folder, recursive = TRUE)
    dir.create(cores_folder)
    
    # write a list of all sst files available
    sst_list <- "sst_list"
    unlink(sst_list)
    write_lines(sst_files, path = sst_list)
    
    # gather details
    dat <- tibble(
      ind0 = seq(1, maxInd, by = core_size),
      ind1 = ind0 + core_size - 1,
      ind = "")
    dat$ind1[nrow(dat)] <- maxInd
    for (i in 1:nrow(dat)) dat$ind[i] = str_c("'", str_c(dat$ind0[i]:dat$ind1[i], collapse = ","), "'")
    
    dat$call <- str_c("sstcoremulti ", sst_list, " ", dat$ind)
    dat$core <- str_c("core", formatC(1:nrow(dat), width = 3, flag = "0"))
    dat$fn   <- str_c(cores_folder, "/", dat$core, ".RDS")
    
    fn0 <- first( sst_files) %>% gsub("[^0-9-]", "", .)
    fn1 <- last(  sst_files) %>% gsub("[^0-9-]", "", .)
    l   <- length(sst_files)
    
    void <- foreach(i = 1:nrow(dat)) %dopar%{
      d <- dat[i,]
      
      if (!file.exists(d$fn)) {
        x <- system(d$call, intern = TRUE, ignore.stderr = TRUE)
        x <- str_split(x, ",") %>% do.call(rbind, .)
        x <- apply(x, 2, as.numeric)
        x <- list(
          ind0 = d$ind0,
          ind1 = d$ind1,
          fn0  = fn0,
          fn1  = fn1,
          l    = l,
          mat  = x)
        
        saveRDS(x, file = d$fn)
      }
    }
    unlink(sst_list)
  }else{
    ## or update cores ----
    
    # number of sst files available (i.e., days)
    slice_n <- sst_files %>% length
    
    # number of entries already in the core files
    x <- readRDS(core_files[1])
    core_n <- x$l
    
    # update core files if new slices are available
    if (slice_n < core_n) stop("critical error:\n  there are more entries in the core files than there are sst slices!!!\n  likely the entire dataset will have to be rebuilt")
    
    if (slice_n == core_n) MSG2 <- "sst cores are already up-to-date\n"
    
    if (slice_n > core_n & slice_n - core_n < 30) MSG2 <- "sst cores not up-to-date, but diff not large enough to justify updating\nif you still want to update the cores run the code manually\n"
    
    if (slice_n - core_n >= 90) {
      n_missing  <- slice_n - core_n
      f_missing1 <- tail(slice_files, n_missing) # missing files using method 1
      f_missing2 <- seq_along(slice_files) > which(slice_files == x$fn1)
      f_missing2 <- slice_files[which(f_missing2)] # missing files using method 2
      if (!identical(f_missing1, f_missing2)) stop("couldn't compute which slice files are missing from the cores")
      
      dat <- list()
      for (i in seq_along(f_missing1)) {
        dat[[i]] <- readBin(con = f_missing1[i], what = "integer", n = 1440 * 720, size = 2)
      }
      dat <- do.call(rbind, dat)
      
      void <- foreach(ii = seq_along(core_files)) %dopar% {
        x <- readRDS(core_files[ii])
        DAT   <- dat[, x$ind0:x$ind1, drop = FALSE]
        x$mat <- rbind(x$mat, DAT)
        x$l   <- nrow(x$mat)
        x$fn1 <- last(f_missing1) %>% gsub("[^0-9-]", "", .)
        saveRDS(x, file = core_files[ii])
      }
      MSG2 <- "sst cores updated\n"
    }
  }
  
  ## core helper data ----
  ## contains data and functions needed for easy use of the sst cores repository
  
  # check if the fn_info1 file exists
  # if yes, load it
  # if not, prepare it now
  if (!file.exists(fn_info2)) {
    cat("building fn_info2 file\n")
    
    id_core <- function(i) {
      i <- sort(i)
      tibble(id = i, core = ceiling(i / 5000), ind = i %% 5000)
    }
    
    apply2coreS <- function(oisst_cores_folder, FUNS, fun_detrend = NULL, smooth = NULL, fullYears = TRUE) {
      n_cores <- length(dir(oisst_cores_folder))
      x  <- readRDS(str_c(oisst_cores_folder, "core001.RDS"))
      d1 <- seq.Date(as_date(x$fn0), as_date(x$fn1), by = "days")
      d2 <- seq.Date(as_date(str_c(str_sub(x$fn0, 1, 5), "01-01")), as_date(str_c(str_sub(x$fn1, 1, 5), "12-31")), by = "days")
      d3 <- d2[!(d2 %in% d1)] %>% str_sub(1, 4) %>% unique
      fullYears <- if (fullYears) {
        !(str_sub(d1, 1, 4) %in% d3)
      }else{
        rep(TRUE, length(d))
      }
      registerDoParallel(min(n_cores, 140))
      all_results <- foreach(i = 1:n_cores, .inorder = FALSE) %dopar% {
        list(i = i, res = apply2core(i, oisst_cores_folder, FUNS, fun_detrend, smooth, fullYears))
      }
      all_results <- purrr::map(all_results, "res")[map_dbl(all_results, "i")]
      ALL_RESULTS <- list()
      for (i in seq_along(FUNS)) {
        vals <- purrr::map(all_results, i)
        ALL_RESULTS[[i]] <- unlist(vals)
      }
      list(n = sum(fullYears), r = ALL_RESULTS)
    }
    
    apply2core <- function(core_number, oisst_cores_folder, FUNS, fun_detrend, smooth, fullYears) {
      core <- str_c(oisst_cores_folder, "core", formatC(core_number, width = 3, flag = "0"), ".RDS")
      if (!file.exists(core)) stop("target file is missing")
      x  <- readRDS(core)
      t0 <- gsub("[^0-9-]", "", x$fn0) %>% as_date
      t1 <- gsub("[^0-9-]", "", x$fn1) %>% as_date
      tt <- tibble(
        date    = seq.Date(t0, t1, "day"),
        date_sh = str_sub(date, 6, 10),
        keep    = date_sh != "02-29")
      x <- x$mat %>% binTOr
      if (nrow(x) != nrow(tt)) stop("dates not properly recognised")
      
      x  <- x [fullYears, ]
      tt <- tt[fullYears, ]
      
      x <- as.data.frame(x) %>% as.list
      
      if (!is.null(fun_detrend)) x <- purrr::map(x, ~fun_detrend(.x, tt, smooth))
      
      results <- list()
      for (i in 1:length(FUNS)) {
        results[[i]] <- map_dbl(x, ~get(FUNS[i])(.x))
      }
      purrr::map(results, unname)
    }
    
    fun_detrend <- function(temp_core, tt, smooth) {
      tt$temp <- temp_core
      clim <- filter(tt, keep) %>%
        select(date_sh, temp) %>%
        group_by(date_sh) %>%
        summarise(avg = mean(temp))
      
      CLIM <- tibble(
        date_sh = c(clim$date_sh, "02-29"),
        avg = c(clim$avg, filter(clim, date_sh %in% c("02-28", "03-01"))$avg %>% mean)) %>%
        arrange(date_sh)
      
      if (!is.null(smooth)) {
        tmp <- rbind(
          add_column(CLIM, yr = rep(1, nrow(CLIM))),
          add_column(CLIM, yr = rep(2, nrow(CLIM))),
          add_column(CLIM, yr = rep(1, nrow(CLIM))))
        tmp$smooth <- rollmean(tmp$avg, k = smooth, fill = NA, align = "center")
        CLIM$avg <- filter(tmp, yr == 2)$smooth
      }
      
      tt$clim <- CLIM$avg[match(tt$date_sh, CLIM$date_sh)]
      as.numeric(tt$temp - tt$clim)
    }
    
    slope_per_decade <- function(vals) {
      vals <- lm(vals ~ seq_along(vals)) %>% coefficients
      vals[2] * 365.25 * 10
    }
    
    se_corr <- function(vals) {
      # standard error of the slope [[corrected]]
      # 
      # correction following: 
      # Foster, G., and Rahmstorf, S. (2011). Global temperature evolution 1979–2010. 
      # Environ Res Lett 6, 044022. doi:10.1088/1748-9326/6/4/044022
      # 
      # decay
      #   lag2_correlation / lag1_correlation
      #
      # effective degrees of freedom [eff_df]
      #   1 + ((2 * lag1_correlation) / (1 - decay))
      #
      # se_corrected
      #   se * sqrt(eff_df)
      
      se <- sd(vals) / sqrt(length(vals))
      lag1_corr <- cor(head(vals, -1), tail(vals, -1))
      lag2_corr <- cor(head(vals, -2), tail(vals, -2))
      decay  <- lag2_corr / lag1_corr
      eff_df <- 1 + ((2 * lag1_corr) / (1 - decay))
      se * sqrt(eff_df)
    }
    
    fun_list <- c("slope_per_decade", "se_corr", "id_core", "fun_detrend", "apply2core", "apply2coreS")
    
    info2 <- tibble(
      item = fun_list,
      info = c(
        "function, slope of the linear regression in degC per decade",
        "function, standard error, corrected for temporal autocorrelation",
        "function, find in which sst_core a given pixel index is stored",
        "function, detrend using doy climatology, with or without smoothing by rolling mean",
        "function, apply a function to all columns of a single sst_core",
        "function, apply a function to all sst_cores"))
    
    fun_list <- c("slope_per_decade", "se_corr")
    
    save(list = c(info2$item, "fun_list", "info2"), file = fn_info2)
  }
  load(fn_info2)
  
  if (MSG2 == "sst cores updated\n") {
    ## compute stats ----
    # slope stats from detrended data using smoothed climatology
    stats <- apply2coreS(oisst_cores_folder = cores_folder, FUNS = fun_list, fun_detrend = fun_detrend, smooth = 30)
    n     <- stats$n
    stats <- stats$r
    names(stats) <- fun_list
    for (i in 1:length(stats)) {
      stat <- list(n = n, vals = stats[[i]])
      saveRDS(stat, file = str_c(out_folder, fun_list[i], ".RDS"))
    }
  }
}

cat(MSG1)
cat(MSG2)
