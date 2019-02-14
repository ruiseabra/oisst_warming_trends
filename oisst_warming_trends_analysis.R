# PROCESS PREVIOUS GLOBAL SST WARMING TRENDS
# AND PRODUCE FINAL PLOTS AND STATS

# - must be run after "oisst_download.R", as it depends on the "info" and "stat" files produced by that script

# This routine requires a binary GSHHG (Global Self-consistent, Hierarchical, High-resolution Geography) database file. 
# The GSHHG database has been released in the public domain and may be downloaded from https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/. 
#   1 - download the binary zip file (e.g. gshhg-bin-2.3.7.zip)
#   2 - create a folder inside this project's folder named "oisst_warming_trends_inputs" (if you following these instructions AFTER trying to run the code and getting an error message, then the forlder already exists)
#   3 - unzip it into the "oisst_warming_trends_inputs" folder
#   4 - done!
# the "low" resolution file provides enough detail for the plots here produced

# NOTE: the computation of distances to land takes a long time (10-15 mins), but after cdone once the result is saved in a raster file and in subsequent runs of the script it is simply loaded back


rm(list = ls())

libraries <- c("stringr", "lubridate", "tibble", "dplyr", "readr", "purrr", "magrittr", "PBSmapping", "maptools", "spatstat", "raster", "igraph", "scales", "Hmisc", "colorspace", "RColorBrewer", "cowplot")
for (l in libraries) library(l, character.only = TRUE)

#### 0----- load ------0 ####

#### options ----
opt <- list()

# geographical range
opt$coords <- c(xmin = -180, xmax = 180, ymin = -60, ymax = 60)

# how far from the coast do we consider the upwelling area?
opt$Lim <- 500 # in km
opt$Bks <- 25  # with of distance breaks

# significance level for the significance level of slope of the warming rate
opt$sig <- 0.05

# function to convert raster data into a ggplot-ready tibble
opt$raster2gg <- function(r, CROP = FALSE) {
  if (CROP) r <- crop(r, opt$coords)
  r <- rasterToPoints(r) %>% as_tibble
  colnames(r) <- c("x", "y", "z")
  r
}

# folders and files
opt$ff <- list(I = "/oisst_warming_trends_inputs/",
               O = "/oisst_warming_trends_outputs/",
               downO  = "/oisst_download_outputs/",
               stats  = "/final_stats/",
               plots  = "/final_plots/") %>% purrr::map(~str_c(getwd(), .x))

dir.create(opt$ff$I, showWarnings = FALSE)
dir.create(opt$ff$O, showWarnings = FALSE)
gshhs <- dir(opt$ff$I, recursive = TRUE)
opt$ff$gshhsI <- dir(opt$ff$I, recursive = TRUE, full.names = TRUE, pattern = "gshhs_l.b")
if (!length(opt$ff$gshhsI)) stop("GSHHS coastline file is missing from the appropriate folder\nPlease refer to the instructions on the beggining of this script for where to download it from and which folder to save it to")
opt$ff$gshhsO = str_c(opt$ff$O, "gshhs_l.RData")

dir.create(opt$ff$plots, showWarnings = FALSE)
dir.create(opt$ff$stats, showWarnings = FALSE)

# load helper data and functions
load(str_c(opt$ff$downO, "oisst_info1.RData"))

#### coastline ----
if (file.exists(opt$ff$gshhsO)) {
  load(opt$ff$gshhsO)
}else{
  coast <- opt$ff$gshhsI %>%
    importGSHHS(xlim = opt$coords[1:2], ylim = opt$coords[3:4], maxLevel = 1) %>%
    PolySet2SpatialPolygons(close_polys = TRUE)
  projection(coast) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
  
  coast_all <- coast@polygons %>%
    purrr::map(~.x@Polygons[[1]]@coords) %>%
    map2(., 1:length(.), ~tibble(x = .x[,1], y = .x[,2], id = .y)) %>%
    do.call(rbind, .)
  
  large_land <- map_lgl(coast@polygons, ~.x@area >= (1 / 10000))
  
  coast <- coast@polygons[large_land] %>%
    purrr::map(~.x@Polygons[[1]]@coords) %>%
    map2(., 1:length(.), ~tibble(x = .x[,1], y = .x[,2], id = .y)) %>%
    do.call(rbind, .)
  
  coast <- list(all = coast_all, small = coast)
  save(coast, file = opt$ff$gshhsO)
}

dat <- list(ggcoastAll = coast$all, ggcoast = coast$small)

#### distance 2 land ----
FILE <- str_c(opt$ff$O, "distance_to_land_area.grd")
if (file.exists(FILE)) {
  distance_to_land <- raster(FILE)
}else{
  original <- rotate(sea_land_mask)
  
  # remove the Great Lakes and inner seas
  inland_waters <- original
  inland_waters[!is.na(values(inland_waters))] <- 1
  inland_waters <- suppressWarnings(clump(inland_waters))
  inland_waters[values(inland_waters) <= 1] <- NA
  inland_waters[!is.na(values(inland_waters))] <- 1
  
  just_sea <- original
  just_sea[values(inland_waters) == 1] <- NA
  
  land <- original
  land[] <- 0
  land[ is.na(just_sea)] <- 1
  land[!is.na(just_sea)] <- NA
  land <- crop(land, extent(opt$coords))
  
  distance_to_land <- distance(land)
  distance_to_land <- round(distance_to_land / 1000)
  distance_to_land[is.na(distance_to_land)] <- 0
  # we ignored small landmasses in order to compute distance to the coast, but distance values on those pixels must remain zero nontheless
  land <- crop(original, extent(opt$coords)) %>% values %>% is.na %>% which
  distance_to_land[land] <- 0
  writeRaster(distance_to_land, FILE, overwrite = TRUE)
}
dat$dist    <- distance_to_land

#### read rates ----
x <- readRDS(str_c(opt$ff$downO, "slope_per_decade.RDS"))
dat$n <- x$n
dat$slope <- rebuild(x$vals) %>% rotate %>% crop(opt$coords)

# se_corrected (is slope significant?)
x <- readRDS(str_c(opt$ff$downO, "se_corr.RDS"))
dat$se_corrected <- rebuild(x$vals) %>% rotate %>% crop(opt$coords)

# sig se_corrected
x <- (dat$slope / dat$se_corrected)
x[] <- pt(abs(x[]), df = dat$n - 2, lower.tail = FALSE) * 2
dat$p <- x
x[] <- ifelse(x[] < opt$sig, 1, 0)
x[dat$slope[] < 0 & x[] == 1] <- -1
dat$p_slope <- x

# area per pixel
x <- area(dat$slope)
x[is.na(dat$slope)] <- NA
dat$area <- x

#### areas ----
## defining EBUS areas and other areas relevant for the analysis
areas <- x <- list()
dat$vars <- vars  <- c("dist", "slope", "p_slope", "p", "area")

### GLOBAL
# Global
# i.e. everywhere
x$ext <- opt$coords
x$nm  <- "Global"
x$nm2 <- "Global\nocean"

for (var in vars) x[[var]] <- crop(dat[[var]], x$ext)
areas[[x$nm]] <- x

# Global non coastal
# i.e. everywhere farther than opt$Lim (500) kms from the coast
x$ext <- opt$coords
x$nm  <- "Global_non_coastal"
x$nm2 <- "Global\non-coastal"
for (var in vars) x[[var]] <- crop(dat[[var]], x$ext)
# coastal data will be set to NA later
areas[[x$nm]] <- x

# Global coastal
# i.e. everywhere at or nearer than opt$Lim (500) kms from the coast
x$ext <- opt$coords
x$nm  <- "Global_coast"
x$nm2 <- "Coastal"
for (var in vars) x[[var]] <- crop(dat[[var]], x$ext)
# set non-coastal pixels to NA
not_coastal <- values(x$dist) == 0 | values(x$dist) > opt$Lim
for (var in vars) x[[var]][not_coastal] <- NA
areas[[x$nm]] <- x

# Global non ebus
# i.e. locations outside EBUS at or nearer than opt$Lim (500) kms from the coast
x$ext <- opt$coords
x$nm  <- "Global_non_ebus"
x$nm2 <- "Outside\nEBUS"
# data inside EBUS will be set to NA after the EBUS are defined
areas[[x$nm]] <- x

# Global ebus
# i.e. locations inside EBUS at or nearer than opt$Lim (500) kms from the coast
x$ext <- opt$coords
x$nm  <- "Global_ebus"
x$nm2 <- "Inside\nEBUS"
# data inside EBUS will be set to NA after the EBUS are defined
areas[[x$nm]] <- x

### EBUS
# i.e. locations inside each EBUS at or nearer than opt$Lim (500) kms from the coast

# California
x$ext <- extent(-136,-110,24,46)
x$nm2 <- x$nm <- "California"
for (var in vars) x[[var]] <- crop(areas$Global_coast[[var]],  x$ext)
# excluding the Gulf of California
remove <- suppressWarnings(clump(x$dist))
remove <- which(values(remove) == 2)
for (var in vars) x[[var]][remove] <- NA
areas[[x$nm]] <- x

# Humboldt
x$ext <- extent(-90,-68,-41,-6)
x$nm2 <- x$nm <- "Humboldt"
for (var in vars) x[[var]] <- crop(areas$Global_coast[[var]],  x$ext)
areas[[x$nm]] <- x

# Canary
x$ext <- extent(-26,-6,14,44)
x$nm2 <- x$nm <- "Canary"
for (var in vars) x[[var]] <- crop(areas$Global_coast[[var]],  x$ext)
# this extent includes sea pixels that should be excluded
# along the northern coast of Spain
ext2 <- extent(-7.617352,-6,43,44)
for (var in vars) x[[var]][cellsFromExtent(x[[var]], ext2)] <- NA
areas[[x$nm]] <- x

# Benguela
x$ext <- extent(1,20,-36,-13)
x$nm2 <- x$nm <- "Benguela"
for (var in vars) x[[var]] <- crop(areas$Global_coast[[var]],  x$ext)
areas[[x$nm]] <- x

# save the indexes of specific areas for later reference
dat$ind$global    <- grep(regex("Global$"), names(areas))
dat$ind$coast     <- grep("Global_coast", names(areas))
dat$ind$non_coast <- grep("Global_non_coastal", names(areas))

# exclude coastal data from the "Global_non_coastal" layer
i <- dat$ind$coast
non_nas <- !is.na(values(areas[[i]]$dist))
for (var in vars) {
  cells_in_coast <- cellsFromExtent(areas$Global_non_coastal[[var]], areas[[i]]$ext)[non_nas]
  areas$Global_non_coastal[[var]][cells_in_coast] <- NA
}

# Excluding ebus from the "Global_non_ebus" layer
dat$ind$non_ebus <- grep("Global_non_ebus", names(areas))
dat$ind$ebus     <- grep("Global", names(areas), invert = TRUE)
for (i in dat$ind$ebus) {
  non_nas <- !is.na(values(areas[[i]]$dist))
  for (var in vars) {
    cells_in_ebus <- cellsFromExtent(areas$Global_non_ebus[[var]], areas[[i]]$ext)[non_nas]
    areas$Global_non_ebus[[var]][cells_in_ebus] <- NA
  }
}

# Excluding coastal data from the "Global_ebus" layer
i <- dat$ind$non_ebus
non_nas <- !is.na(values(areas[[i]]$dist))
for (var in vars) {
  cells_out_ebus <- cellsFromExtent(areas$Global_ebus[[var]], areas[[i]]$ext)[non_nas]
  areas$Global_ebus[[var]][cells_out_ebus] <- NA
}

# uniform plot areas
# ext2 ensures all EBUS are plotted with the same lon and lat width, to maintain approximate ratio and size
x <- areas[dat$ind$ebus] %>% purrr::map("ext") %>% lapply(as.vector) %>% as_tibble
loc <- colnames(x)
x <- x %>% t %>% as_tibble
colnames(x) <- areas[[1]]$ext %>% names
x$loc <- loc
x <- x %>% mutate(x = xmax - xmin, y = ymax - ymin, xc = xmin + x / 2, yc = ymin + y / 2)
X <- max(x$x)
Y <- max(x$y)
x <- x %>% mutate(xmin = xc - X / 2, xmax = xc + X / 2, ymin = yc - Y / 2, ymax = yc + Y / 2)

ind <- dat$ind$ebus
for (i in 1:length(ind)) areas[[ind[i]]]$ext2 <- dplyr::select(x[i,], xmin, xmax, ymin, ymax) %>% as_vector %>% extent

dat$areas <- areas

# rectangles for signaling EBUS in plots
x <- purrr::map(dat$areas[dat$ind$ebus], "ext")
dat$rect <- tibble(
  loc  = names(x),
  xmin = map_dbl(x, ~.x@xmin),
  xmax = map_dbl(x, ~.x@xmax),
  ymin = map_dbl(x, ~.x@ymin),
  ymax = map_dbl(x, ~.x@ymax))

#### stats ----
x <- tibble(
  loc  = map_chr(dat$areas, "nm"),
  loc2 = map_chr(dat$areas, "nm2"))

# area-weighted computation of stats ensures that pixels at higher latitudes are not over-represented
funs <- c(
  sd  = function(x) sqrt(wtd.var(na.omit(x$slope[]), na.omit(x$area[]))),
  min = function(x) min(na.omit(x$slope[])),
  q25 = function(x) wtd.quantile(na.omit(x$slope[]), na.omit(x$area[]), 0.25),
  avg = function(x) wtd.mean(    na.omit(x$slope[]), na.omit(x$area[])),
  q75 = function(x) wtd.quantile(na.omit(x$slope[]), na.omit(x$area[]), 0.75),
  max = function(x) max(na.omit(x$slope[])))

for (i in seq_along(funs)) x[[names(funs[i])]] <- map_dbl(dat$areas, ~funs[[i]](.x))

dat$stats <- x

#### gather data ----
# gather global, non_ebus and ebus data into one single tibble
ind <- c(dat$ind$non_coast, dat$ind$non_ebus, dat$ind$ebus)
df  <- list()
for (i in 1:length(ind)) {
  x <- list()
  nas    <- dat$areas[[ind[i]]]$slope %>% values %>% is.na
  coords <- dat$areas[[ind[i]]]$slope %>% coordinates %>% as_tibble %>% filter(!nas)
  for (v in dat$vars) x[[v]] <- values(dat$areas[[ind[i]]][[v]])[!nas]
  x$lat    <- coords$y
  x$lon    <- coords$x
  x$area   <- dat$areas[[ind[i]]]$slope %>% (raster::area) %>% values %>% .[!nas] %>% round
  x <- as_tibble(x)
  x$id <- ind[i]
  df[[i]] <- x
}
df <- do.call("rbind", df) %>% mutate(area_ratio = area / max(area))

# add location specifiers
loc    <- set_names(map_chr(dat$areas, "nm"), seq(dat$areas %>% length))
df$loc <- loc[df$id]

# ebus or not ebus (logical)
df$ebus <- ifelse(df$id == dat$ind$non_coast, NA, df$id %in% dat$ind$ebus)

# distance breaks
breaks <- seq(0, opt$Lim, by = opt$Bks)
nms    <- head(breaks, -1) + (opt$Bks / 2)
bks    <- cut(df$dist, breaks = breaks)
bks    <- nms[as.numeric(bks)]
df$bks <- bks

# capped version (slope)
# q <- wtd.quantile(x = df$slope, probs = c(0.025, 0.975), weight = df$area_ratio) %>% abs %>% max
# q == 0.3526074, so we set q to 0.4 to tidy the layout of color scales
q <- 0.4
y <- ifelse(df$slope > q, q, df$slope)
y <- ifelse(y < -q, -q, y)
df$fill_slope <- y

# the percentile of each slope against all COASTAL slopes
# quantiles are also area-weighted
qt <- ewcdf(df$slope, weights = df$area_ratio)
df$qt_slope <- qt(df$slope) %>% rescale

dat$df <- df

#### 0----- plots ------0 ####
colsABCD <- brewer.pal(4, "Set1")
COLS1    <- diverge_hsv(100)

Dist     <- "Distance to land (km)"
Breaks   <- c(0, opt$Lim / 2, opt$Lim)


#### map global ####
ABCD  <- mutate(dplyr::select(dat$rect, loc, xmin, ymax), lab = c("Ca", "Hu", "Cn", "Be"))
colnames(ABCD) <- c("loc", "x", "y", "lab")
ABCD$x <- ABCD$x - 3

p <- ggplot(dat$df) + 
  geom_raster(aes(lon, lat, fill = fill_slope)) + 
  geom_raster(data = filter(dat$df, !p_slope), aes(lon, lat), fill = "grey") +
  stat_contour(aes(lon, lat, z = p), breaks = opt$sig, size = 0.1, col = "black") +
  geom_polygon(data = dat$ggcoast, aes(x = x, y = y, group = id), color = NA, fill = "gray20") +
  scale_fill_gradientn(colours = COLS1, name = "Warming\nrate\n(ºC/dec)") +
  stat_contour(data = opt$raster2gg(dat$dist), aes(x, y, z = z), breaks = opt$Lim, color = "black", size = 0.1) +
  geom_rect(data = dat$rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = colsABCD, size = 0.75) +
  geom_text(data = ABCD, aes(x = x, y = y, label = lab), color = colsABCD, size = 3, vjust = 1, hjust = 1, fontface = "bold") +
  theme_minimal() + 
  coord_cartesian() +
  xlab("") + ylab("") +
  scale_x_discrete(limits = seq(opt$coords[1], opt$coords[2], by = 30)) +
  scale_y_discrete(limits = seq(opt$coords[3], opt$coords[4], by = 30)) +
  theme(legend.key.width  = unit(3, "mm"), 
        legend.key.height = unit(4, "mm"), 
        legend.title = element_text(size = 8), 
        legend.text  = element_text(size = 7),
        legend.background = element_rect(colour = NA, fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(str_c(opt$ff$plots, "fig1.pdf"), p, device = "pdf", width = 7, height = 3)

#### map detail ####
ebus.fun1 <- function(DAT, i, COL, sig) {
  AREA <- dat$areas[[dat$ind$ebus[i]]]
  x <- DAT[DAT$id == dat$ind$ebus[i], ]
  thisCOL <- rescale(COL, from = x$fill %>% range)

  lim <- opt$raster2gg(crop(dat$dist, AREA$ext))
  
  p <- ggplot(x) + 
    coord_cartesian(xlim = AREA$ext2[1:2], ylim = AREA$ext2[3:4], expand = FALSE) +
    scale_x_discrete(limits = pretty(seq(AREA$ext2[1], AREA$ext2[2], by = 10), 2)) +
    scale_y_discrete(limits = pretty(seq(AREA$ext2[3], AREA$ext2[4], by = 10), 3)) +
    theme_minimal() + 
    xlab("") + ylab("") +
    geom_raster(aes(lon, lat, fill = fill))
  
  if (sig) p <- p +
    geom_raster(data = filter(x, !p_slope), aes(lon, lat), fill = "grey") +
    stat_contour(aes(x = lon, y = lat, z = p), breaks = opt$sig, size = 0.4, color = "black")
  
  p <- p +
    geom_polygon(data = dat$ggcoastAll, aes(x = x, y = y, group = id), color = NA, fill = "gray20") +
    stat_contour(data = lim, aes(x, y, z = z), breaks = opt$Lim, color = "darkgrey", size = 0.4) +
    theme(
      plot.title   = element_text(size = 10, colour = colsABCD[i], hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA))
  
  list(p = p, col = thisCOL, nm = AREA$nm, x = x)
}

ebus.fun2 <- function(DAT, COL, sig) {
  X <- list()
  
  x <- ebus.fun1(DAT, 1, COL, sig)
  COL1 <- x$col
  p <- x$p +
    scale_fill_gradientn(colours = names(COL1), values = COL1, name = "", guide = FALSE)
  X[[x$nm]] <- p
  
  x <- ebus.fun1(DAT, 2, COL, sig)
  COL2 <- x$col
  p <- x$p +
    scale_fill_gradientn(colours = names(COL2), values = COL2, name = "", guide = FALSE)
  X[[x$nm]] <- p
  
  x <- ebus.fun1(DAT, 3, COL, sig)
  COL3 <- x$col
  p <- x$p +
    scale_fill_gradientn(colours = names(COL3), values = COL3, name = "", guide = FALSE)
  X[[x$nm]] <- p
  
  x <- ebus.fun1(DAT, 4, COL, sig)
  COL4 <- x$col
  p <- x$p +
    scale_fill_gradientn(colours = names(COL4), values = COL4, name = "", guide = FALSE)
  X[[x$nm]] <- p
  
  X <- list(maps = X)
  
  # legend
  p <- ggplot(DAT) + 
    geom_raster(aes(lon, lat, fill = fill)) +
    scale_fill_gradientn(colours = names(COL), values = seq(0, 1, length.out = length(COL)), name = "") +
    theme_void() + 
    theme(
      legend.key.width = unit(3, "mm"),
      legend.key.height = unit(4, "mm"), 
      legend.title = element_text(size = 8), 
      legend.text = element_text(size = 7))
  
  X$legend <- get_legend(p)
  
  X
}

map.ebus <- function(type, sig, COLS) {
  DAT  <- dat$df
  DAT$fill <- DAT[[type]]
  if (type == "qt_slope") DAT$fill <- DAT$fill * 100
  colSEQ <- set_names(seq(min(DAT$fill), max(DAT$fill), length.out = 100), COLS)
  ebus.fun2(DAT, colSEQ, sig)
}

FIG2 <- list()
FIG2$slope <- map.ebus(type = "fill_slope", sig = TRUE,  COLS = COLS1)
FIG2$quant <- map.ebus(type = "qt_slope",   sig = FALSE, COLS = COLS1)

#### warming rates ####
warming.v.distance.fun <- function(x3, COL, YLIM) {
  ggplot(x3) + 
    geom_hline(yintercept = 0, col = "black") +
    geom_ribbon(aes(bks, ymin = q1, ymax = q3, fill = col, group = id), color = NA, alpha = 0.2) +
    geom_line(aes(bks, q2, color = col, group = id)) +
    scale_color_manual(values = c("black", COL)) + 
    scale_fill_manual( values = c("black", COL)) + 
    scale_x_continuous(breaks = Breaks) +
    coord_cartesian(xlim = c(0, opt$Lim), ylim = YLIM) +
    theme_classic() + 
    guides(color = FALSE, fill = FALSE, size = FALSE) + 
    labs(x = Dist, y = "", title = "Warming rates\n(°C/decade)") + 
    theme(plot.title = element_text(size = 9, hjust = 0.5), 
          axis.title = element_text(size = 6, hjust = 0.5),
          axis.ticks.length = unit(-0.1, "cm"), 
          axis.text.x  = element_text(size = 7, margin = margin(0.2, 0, 0, 0, "cm")), 
          axis.text.y  = element_text(size = 7, margin = margin(0, 0.2, 0, 0, "cm")),
          panel.border = element_rect(colour = "black", fill = NA))
}

DAT <- filter(dat$df, loc != "Global_non_coastal")
x <- DAT %>% 
  mutate(col = as.character(id - min(DAT$id) + 1)) %>% 
  group_by(id, bks) %>% 
  summarise(
    loc = first(loc),
    col = first(col),
    q1  = wtd.quantile(slope, area, 0.25),
    q2  = wtd.quantile(slope, area, 0.50),
    q3  = wtd.quantile(slope, area, 0.75))
YLIM <- range(x$q1, x$q3)

EBUS <- names(dat$areas[dat$ind$ebus])

p <- list()

i <- 1
x_subset <- filter(x, loc == EBUS[i] | loc == "Global_non_ebus")
COL1 <- colsABCD[i]
p[[i]] <- warming.v.distance.fun(x_subset, COL1, YLIM)

i <- 2
x_subset <- filter(x, loc == EBUS[i] | loc == "Global_non_ebus")
COL2 <- colsABCD[i]
p[[i]] <- warming.v.distance.fun(x_subset, COL2, YLIM)

i <- 3
x_subset <- filter(x, loc == EBUS[i] | loc == "Global_non_ebus")
COL3 <- colsABCD[i]
p[[i]] <- warming.v.distance.fun(x_subset, COL3, YLIM)

i <- 4
x_subset <- filter(x, loc == EBUS[i] | loc == "Global_non_ebus")
COL4 <- colsABCD[i]
p[[i]] <- warming.v.distance.fun(x_subset, COL4, YLIM)

FIG2$slope_dist <- p

wrm  <- FIG2$slope_dist
slp  <- FIG2$slope$maps
slpL <- FIG2$slope$legend
qnt  <- FIG2$quant$maps
qntL <- FIG2$quant$legend

p <- plot_grid(wrm[[1]], wrm[[2]], wrm[[3]], wrm[[4]], NULL, NULL,
               slp[[1]], slp[[2]], slp[[3]], slp[[4]], NULL, slpL,
               qnt[[1]], qnt[[2]], qnt[[3]], qnt[[4]], NULL, qntL, 
               ncol = 6, rel_widths = c(1,1,1,1,0.5,1), 
               labels = c("A", "B", "C", "D", "", "", "E", "F", "G", "H", "", "", "I", "J", "K", "L", "", ""))
ggsave(str_c(opt$ff$plots, "fig2.pdf"), p, device = "pdf", width = 7, height = 7)


#### Cn + Be rates ####
XLIM <- c(-22, -14)
YLIM <- c( 17,  24)
x <- crop(dat$areas$Canary$slope, c(XLIM, YLIM))
MIN <- which.min(x[]) %>% xyFromCell(x, .)
MAX <- which.max(x[]) %>% xyFromCell(x, .)
minmax <- rbind(MIN, MAX) %>% as_tibble
txt <- str_c(
  "min: ", 
  raster::extract(x, MIN) %>% round(3), 
  ", max: ", 
  raster::extract(x, MAX) %>% round(3),
  ", linear dist = ",
  (pointDistance(MIN, MAX, lonlat = TRUE) / 1000)  %>% round,
  " km")
p <- ggplot(opt$raster2gg(x)) + 
  coord_fixed(xlim = XLIM, ylim = YLIM) +
  geom_raster(aes(x, y, fill = z)) +
  scale_fill_gradientn(colours = COLS1) +
  geom_polygon(data = dat$ggcoast, aes(x = x, y = y, group = id), color = "black", fill = "gray80") +
  geom_point(data = minmax, aes(x, y), size = 2) +
  ggtitle(txt) +
  xlab("") + ylab("") +
  guides(fill = FALSE)
ggsave(str_c(opt$ff$plots, "minmax_Canary.pdf"), p, device = "pdf")

XLIM <- dat$areas$Benguela$ext[1:2]
YLIM <- dat$areas$Benguela$ext[3:4]
x <- crop(dat$areas$Benguel$slope, c(XLIM, YLIM))
MIN <- which.min(x[]) %>% xyFromCell(x, .)
MAX <- which.max(x[]) %>% xyFromCell(x, .)
minmax <- rbind(MIN, MAX) %>% as_tibble
txt <- str_c(
  "min: ", 
  raster::extract(x, MIN) %>% round(3), 
  ", max: ", 
  raster::extract(x, MAX) %>% round(3),
  ", linear dist = ",
  (pointDistance(MIN, MAX, lonlat = TRUE) / 1000)  %>% round,
  " km")
p <- ggplot(opt$raster2gg(x)) + 
  coord_fixed(xlim = XLIM, ylim = YLIM) +
  geom_raster(aes(x, y, fill = z)) +
  scale_fill_gradientn(colours = COLS1) +
  geom_polygon(data = dat$ggcoast, aes(x = x, y = y, group = id), color = "black", fill = "gray80") +
  geom_point(data = minmax, aes(x, y), size = 2) +
  ggtitle(txt) +
  xlab("") + ylab("") +
  guides(fill = FALSE)
ggsave(str_c(opt$ff$plots, "minmax_Benguela.pdf"), p, device = "pdf")

#### 0----- stats ------0 ####

#### ocean warming rates ####
# Average warming rates of 
#   a) whole ocean (including coasts), 
#   b) coastal, non-ebus,
#   c) coastal, ebus
tab <- dplyr::select(dat$stats, min, q25, avg, sd, q75, max) %>% round(4)
tab <- cbind(dat$stats$loc, tab)
colnames(tab)[1] <- "loc"
write_csv(tab, str_c(opt$ff$stats, "average_warming_rates.csv"))

#### warming rates vs dist ####
# Warming rates vs distance to coast, for...
#   a) the nearest pixels, 
#   b) the pixels with lowest and highest rates,
#   c) the farthest pixels
DAT <- filter(dat$df, loc != "Global_non_coastal")
tab <- DAT %>% 
  group_by(id, bks) %>% 
  summarise(loc = first(loc),
            avg = wtd.mean(slope, area_ratio) %>% round(2),
            sd  = sqrt(wtd.var(slope, area_ratio)) %>% round(2)) %>%
  group_by(id) %>% 
  summarise(
    loc     = first(loc),
    km25    = first(avg), 
    km25sd  = first(sd), 
    min     = min(avg), 
    max     = max(avg), 
    km500   = last(avg),
    km500sd = last(sd)) %>% 
  dplyr::select(-id)

tab2 <- DAT %>% 
  filter(loc != "Global_non_ebus") %>% 
  group_by(bks) %>% 
  summarise(loc = first(loc),
            avg = wtd.mean(slope, area_ratio) %>% round(2),
            sd  = sqrt(wtd.var(slope, area_ratio)) %>% round(2)) %>%
  summarise(
    loc     = first(loc),
    km25    = first(avg), 
    km25sd  = first(sd), 
    min     = min(avg), 
    max     = max(avg), 
    km500   = last(avg),
    km500sd = last(sd))

tab2$loc <- "EBUS"

tab <- rbind(tab, tab2)

write_csv(tab, str_c(opt$ff$stats, "warming_rates_vs_dist.csv"))

#### . potential warming ####
# Temperature increase during 37yrs if same nearshore warming rate as Global_non_ebus
DAT <- filter(dat$df, loc != "Global_non_coastal")

wrEBUS <- DAT %>% 
  filter(bks == opt$Bks / 2) %>% 
  group_by(id) %>% 
  summarise(loc = first(loc),
            min = min(slope),
            avg = wtd.mean(slope, area_ratio),
            sd  = sqrt(wtd.var(slope, area_ratio)),
            max = max(slope)) %>%
  dplyr::select(-id)
wrALLebus <- DAT %>% 
  filter(loc != "Global_non_ebus", bks == opt$Bks / 2) %>% 
  summarise(loc = first(loc),
            min = min(slope) %>% round(2),
            avg = wtd.mean(slope, area_ratio),
            sd  = sqrt(wtd.var(slope, area_ratio)),
            max = max(slope))
wrALLebus$loc <- "EBUS"
wr <- rbind(wrEBUS, wrALLebus)

# generate n estimates of slope from the avg ± sd slope of ref (global non ebus) and the avg ± sd slope of an EBUS
# for each n pairs, use the resulting delta-slope to compute the potential additional warming
# finally, summarise the n estimates of potential additional warming as avg ± sd potential warming
x <- readRDS("oisst_repository/cores/core001.RDS")
d1 <- seq.Date(as_date(x$fn0), as_date(x$fn1), by = "days")
d2 <- seq.Date(as_date(str_c(str_sub(x$fn0, 1, 5), "01-01")), as_date(str_c(str_sub(x$fn1, 1, 5), "12-31")), by = "days")
d3 <- str_sub(d1, 1, 4) %>% unique
d4 <- d2[!(d2 %in% d1)] %>% str_sub(1, 4) %>% unique
d5 <- d3[!(d3 %in% d4)]
yrs <- c(first(d5), last(d5)) %>% as.numeric
yrs <- yrs[2] - yrs[1] + 1

wr$t_sd <- wr$t_avg <- 0

n   <- 100000
ref <- filter(wr, loc == "Global_non_ebus")
set.seed(1)
ref <- rnorm(n, ref$avg, ref$sd)
for (i in 2:nrow(wr)) {
  sub <- wr[i,]
  set.seed(i)
  sub <- rnorm(n, sub$avg, sub$sd)
  dif <- ref - sub
  dif <- dif * yrs / 10 # because warming rates are in degrees / decade
  wr$t_avg[i] <- mean(dif)
  wr$t_sd[i]  <- sd(dif)
}

wr <- mutate(wr, t_lo = t_avg - t_sd, t_hi = t_avg + t_sd)
wr <- mutate_at(wr, 2:5, ~round(.x, 2))
wr <- mutate_at(wr, 6:ncol(wr), ~round(.x, 1))

write_csv(wr, str_c(opt$ff$stats, "potential_warming.csv"))

#### . cumulative density ####
# slope_sum50 = warming slope when cumulative area proportion >= 0.5
# sum... = cumulative area proportion when warming slope >= 0 or the mean Ocean warming slope
x <- DAT
x$slope <- round(x$slope, 3)
X <- list()
for (i in 1:n_distinct(x$ebus)) {
  XLIM <- tibble(slope = c(-Inf, +Inf), cumsum = 0:1, col = rep(as.character(i), 2))
  X[[i]] <- filter(x, ebus == unique(x$ebus)[i]) %>%
    arrange(slope) %>%
    mutate(cumsum = cumsum(area) / sum(area), col = as.character(i)) %>%
    dplyr::select(slope, cumsum, col) %>%
    rbind(XLIM) %>%
    arrange(slope)
}
x <- do.call(rbind, X)
x$cumsum <- x$cumsum * 100


# slope at which cumsum becomes equal or greater than 0.5
x1 <- filter(x, cumsum >= 50) %>% group_by(col) %>% summarise(slope = head(slope, 1))
# cumsum at which slope becomes equal or greater than 0
x2 <- filter(x, slope >= 0) %>% group_by(col) %>% summarise(cumsum = head(cumsum, 1))
Ocean <- values(dat$slope) %>% mean(na.rm = TRUE) %>% round(2)
x3 <- filter(x, slope >= Ocean) %>% group_by(col) %>% summarise(cumsum = head(cumsum, 1))

X <- cbind(x1, slope0 = x2$cumsum, slopeOcean = x3$cumsum) %>% as_tibble
nm <- colnames(X)
nm[2] <- "cumsum50"
colnames(X) <- nm
X$ebus <- X$col %>% as.numeric %>% unique(DAT$ebus)[.]

tab <- mutate(X,
              slope_sum50  = round(cumsum50, 2),
              sum_slope0   = round(slope0, 1),
              sum_slopeAll = round(slopeOcean, 1)) %>%
  dplyr::select(-(1:4))
txt <- "Cumulative area"

write_csv(wr, str_c(opt$ff$stats, "cumulative_area.csv"))

#### . ebus slope quantiles ####
# slope quantiles inside each EBUS, for the entire area or just the nearshore pixels
stat.quantiles <- function(near = FALSE) {
  ind <- dat$ind$ebus
  DAT <- dat$df
  if (near) DAT <- filter(DAT, bks == opt$Bks / 2)
  tab <- group_by(DAT, loc) %>% summarise(
    min = min(qt_slope), 
    q25 = wtd.quantile(qt_slope, area_ratio, 0.25), 
    avg = wtd.mean(qt_slope, area_ratio),
    sd  = sqrt(wtd.var(qt_slope, area_ratio)),
    q75 = wtd.quantile(qt_slope, area_ratio, 0.75),
    max = max(qt_slope),
    aCool50 = sum(area[qt_slope < 0.5]) / sum(area)) 
  tab <- dplyr::select(tab, -loc) %>% 
    multiply_by(100) %>%
    round(1) %>%
    add_column(loc = tab$loc) %>%
    dplyr::select(loc, min, q25, avg, sd, q75, max, aCool50)
  
  tab2 <- filter(DAT, id %in% ind) %>% summarise(
    min = min(qt_slope), 
    q25 = wtd.quantile(qt_slope, area_ratio, 0.25), 
    avg = wtd.mean(qt_slope, area_ratio),
    sd  = sqrt(wtd.var(qt_slope, area_ratio)),
    q75 = wtd.quantile(qt_slope, area_ratio, 0.75),
    max = max(qt_slope),
    aCool50 = sum(area[qt_slope < 0.5]) / sum(area)) %>% 
    do.call(rbind, .) %>%
    multiply_by(100) %>%
    round(1)
  
  rbind(tab, c("EBUS", as.numeric(tab2))) %>% as_tibble
}

write_csv(stat.quantiles(),     str_c(opt$ff$stats, "slope_quantiles_all.csv"))
write_csv(stat.quantiles(TRUE), str_c(opt$ff$stats, "slope_quantiles_near.csv"))
