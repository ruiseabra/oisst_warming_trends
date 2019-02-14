# CHECK IF YOUR SYSTEM HAS ALL THE REQUIREMENTS
# TO RUN THIS ANALYSIS

## parallel computing ----
# The code was designed to run on a 160 thread local server with 1.5 TB of RAM. 
# It makes frequent use of parallel computations and deals with several "dificult" computations simply by using brute-force. 
# Please bear in mind that if this code is to be run in a regular machine more efficient memory management and a lot of patience may be required...

### NECESSARY ACTION:
### adjust the number of cores in the "registerDoParallel" statement to match your machine's


## sstcoremulti ----
# "sstcoremulti" is a small C++ program to eficiently transform time slices (daily oisst data) into pixel time-series (oisst cores)
# "sstcoremulti" can be found in the same GitHub repository --> check the sstcoremulti folder from the code and instruction on how to compile it
# if you don't want to install (or can't) "sstcoremulti" in your machine then you need to edit the code so that you produce 1 time-series of all oisst data per ocean pixel
has.sstcoremulti <- tryCatch(system("sstcoremulti", intern = TRUE, ignore.stderr = TRUE))
has.sstcoremulti <- !inherits(has.sstcoremulti, "try-error")
if(!has.sstcoremulti) {
  stop("sstcoremulti not available\nsst cores can't be built or maintained but oisst download can proceed\ncompile sstcoremulti to make sst cores")
}else{
  message("your machine has sstcoremulti installed")
}

## gshhs coastline ----
# a binary GSHHG (Global Self-consistent, Hierarchical, High-resolution Geography) database file is required for plotting land
# the GSHHG database has been released in the public domain and may be downloaded from https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/
has.gshhs <- dir("oisst_warming_trends_inputs/", recursive = TRUE, full.names = TRUE, pattern = "gshhs_l.b")
if (!length(has.gshhs)) {
  stop("GSHHS coastline file is missing from the appropriate folder\nPlease refer to the instructions bellow for where to download it from and which folder to save it to")
}else{
  message("your machine has gshhs coastline data in the correct folder")
}
#   1 - download the binary zip file (e.g. gshhg-bin-2.3.7.zip)
#   2 - create a folder inside this project's folder named "oisst_warming_trends_inputs" (if you following these instructions AFTER trying to run the code and getting an error message, then the forlder already exists)
#   3 - unzip it into the "oisst_warming_trends_inputs" folder
#   4 - done!
# the "low" resolution file provides enough detail for the plots here produced

libraries <- c("doParallel", "stringr", "lubridate", "tibble", "dplyr", "readr", "purrr", "magrittr", "PBSmapping", "maptools", "spatstat", "xts", "ncdf4", "raster", "igraph", "scales", "Hmisc", "colorspace", "RColorBrewer", "cowplot")
libraries <- libraries[!(libraries %in% rownames(installed.packages()))]
if (length(libraries)) {
  stop(str_c("the following libraries are missing from this machine\n  ", str_c(libraries, collapse = ", "), "\nthey must be installed before proceeding"))
}else{
  message("all the necessary libraries are installed in this machine")
}
