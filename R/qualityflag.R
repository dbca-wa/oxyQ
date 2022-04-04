#' Quality flagger function
#'
#' \code{flagR} assesses the raw field sample data from a sonde and flags with a
#'       numerical code any suspect samples in an output csv. Suspect issues are:
#'       \itemize{
#'         \item Sample stability check (depth/Vpos +/-10% variation) - `QUAL`
#'         code 01000
#'         \item Sample site check for out of sequence or duplicate (within 0.1m)
#'         records - `QUAL` codes 00100 and 00010
#'         \item Transect check for expected and or incorrectly named samples
#'         - `QUAL` cade 00001 }
#'
#' The \code{flagR} function output csv has the original GPS co-ordinates
#' transformed from WGS84 to GDA94 MGA50. All eastings and northings are
#' therefore measured in metres. \code{flagR} also outputs an Esri shape file
#' with an attribute table created from the flagged raw data input. Lastly,
#' salinity, temperature and density profile plots, facetted by transect are
#' output in png format. All outputs are written to the working directory.
#'
#' @param file A character string of full file path to raw csv data.
#'
#' @return A csv of the original data with quality issues flagged, a shape file
#'       of the same and salinity, temperature and density profile plots.
#'
#' @import readr
#' @import sf
#'
#' @export
#'
#' @examples
#' \dontrun{
#' flagR("./my_raw_sonde_data.csv")
#' }
flagR <- function(file){
  x <- imp_raw_csv(file)

  x <- old_data_updatR(x)

  x <- depth_checkR(x)

  x <- site_checkR(x)

  x <- transect_checkR(x)

  profile_plotR(x, file)

  stub <- basename(file) %>%
    gsub(pattern = ".csv", replacement = "_R1.csv", file)

  csv_fp <- paste0("./", stub)
  shp_fp <- gsub(pattern = ".csv", replacement = ".shp", stub)

  # write csv R version
  x %>%
    readr::write_csv(file = csv_fp)

  # write shape file
  x %>%
    sf::st_as_sf(coords = c("easting_m", "northing_m"), crs = 28350) %>%
    sf::st_write(dsn = shp_fp)
}

#' Import raw csv data from sonde
#'
#' @param file A character string of full file path to raw csv data.
#'
#' @return A tibble of the raw data with a `QUAL` variable set to 10000.
#'
#'
#' @importFrom  readr read_csv locale
#' @importFrom janitor clean_names
#' @importFrom mgsub mgsub
#' @import sf
#' @import dplyr
imp_raw_csv <- function(file){
  # dlist <- vector(mode = "list", length = length(file))
  for(i in seq_along(file)){
    f <- file[i]
    d <- readr::read_csv(f, locale = readr::locale(encoding="latin1")) %>%
      janitor::clean_names()
    # function to split on last '_' in col names and junk sonde id
    fn <- function(x){strsplit(x, "_(?!.*_)", perl=TRUE)[[1]]}
    # rename cols 1st pass
    n <- names(d)
    names(d) <- sapply(lapply(n, FUN = fn), "[[", 1)
    d <- d %>%
      sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      sf::st_transform(crs = 28350) %>%
      dplyr::mutate(x = sf::st_coordinates(.)[,1],
                    y = sf::st_coordinates(.)[,2]) %>%
      sf::st_set_geometry(NULL)
    # rename cols 2nd pass
    names(d) <- mgsub::mgsub(names(d), pattern = var_old, replacement = var_new)
    # nm <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(f))
    d <- d %>%
      dplyr::mutate(QUAL = 10000)
    # dlist[[i]] <- d
    # names(dlist)[i] <- nm
  }
  return(d)
}

#' Update old sample site names to current versions
#'
#' @param x A tibble output from running \code{\link{imp_raw_csv}}
#'
#' @return A tibble with site names in current format
#'
#' @importFrom dplyr mutate
old_data_updatR <- function(x){
  x <- x %>%
    dplyr::mutate(site = ifelse(grepl("c0", site),
                                paste0("cs",
                                       substr(site, 2, 3)),
                                site),
                  site = ifelse(grepl("go", site),
                                paste0("gs0",
                                       substr(site, 4, 4)),
                                site))
}

#' Helper function to reformat a file name.
#'
#' @param file A character string of full file path to raw csv data.
#'
#' @return A more human friendly plot title.
#'
#' @importFrom lubridate as_date
#' @importFrom stringr str_split
#' @importFrom tools file_path_sans_ext
plot_nameR <- function(file){
  date <- lubridate::as_date(substr(basename(file), 1, 8))
  loc <- stringr::str_split(tools::file_path_sans_ext(basename(file)), "_")[[1]][2]
  sname <- paste(date, loc)
  return(sname)
}

#' Checks stability of sample by comparison of depth to vpos measurements on sonde.
#'
#' @param x A tibble as exported from either \code{\link{imp_raw_csv}}
#'       or \code{\link{old_data_updatR}}
#'
#' @return A tibble where `QUAL` variable has 1000 added if > +/-10% variation.
#'
#' @import dplyr
depth_checkR <- function(x){
  out <- x %>%
    dplyr::mutate(percentage = round(((DEPTH_m - VPOS_m) / VPOS_m * 100), 2),
                  ratio = round((DEPTH_m / VPOS_m), 2),
                  QUAL = case_when(
                    abs(percentage) >= 10  ~ QUAL + 1000,
                    TRUE ~ QUAL)) %>%
    dplyr::select(-ratio, -percentage)
  return(out)
}

#' Checks site samples for correct sequence and duplicate depth measurements.
#'
#' @param x A tibble as exported from either \code{\link{imp_raw_csv}}
#'       or \code{\link{old_data_updatR}}
#'
#' @return A tibble where `QUAL` variable has 100 added for an out of sequence
#'       sample and/or 10 added for a duplicate.
#'
#' @import dplyr
site_checkR <- function(x){
  out <- dplyr::tibble()
  sites_to_do <- unique(x[["site"]])

  for(i in seq_along(sites_to_do)){
    s = sites_to_do[i]
    d <- x %>% dplyr::filter(site == s) %>% dplyr::pull(DEPTH_m)
    d <- c(d, 0)
    odf <- x %>%
      dplyr::filter(site == s) %>%
      dplyr::mutate(dd = diff(d),
                    dup = abs(-DEPTH_m - lead(-DEPTH_m)),
                    QUAL = ifelse(dd > 0, QUAL + 100, QUAL),
                    QUAL = case_when(
                      dup < 0.1 ~ QUAL + 10,
                      is.na(dup) ~ QUAL,
                      TRUE ~ QUAL)) %>%
      dplyr::select(-dd, -dup)

    out <- dplyr::bind_rows(out, odf)
  }
  return(out)
}


#' Checks named samples are from expected transects and if named correctly.
#'
#' @param x A tibble as exported from either \code{\link{imp_raw_csv}}
#'       or \code{\link{old_data_updatR}}
#'
#' @return A tibble where `QUAL` variable has 1 added for any sample not spatially
#'       located in appropriate transect or sample has been named incorrectly.
#'
#' @import sf
#' @import dplyr
transect_checkR <- function(x){
  # turn csv to shp
  x_shp <- sf::st_as_sf(x, coords = c("easting_m", "northing_m"), crs = 28350)

  # find outside extents and label tsects
  tcheck <- x_shp %>%
    dplyr::mutate(in_ext = lengths(st_within(x_shp, tpolys))) %>%
    sf::st_join(tpolys, join = st_within) %>%
    dplyr::select(-id, -grp, -buf) %>%
    dplyr::mutate(QUAL = case_when(
      is.na(tsect) ~ QUAL + 1,
      TRUE ~ QUAL)) %>%
    dplyr::mutate(tsamp = substr(tolower(site), 1, 3),
                  QUAL = case_when(
                    tsamp != tolower(tsect) ~ QUAL + 1,
                    TRUE ~ QUAL
                  ),
                  easting_m = sf::st_coordinates(.)[,1],
                  northing_m = sf::st_coordinates(.)[,2]) %>%
    sf::st_set_geometry(NULL) %>%
    dplyr::select(-tsamp, -tsect, -in_ext) %>%
    dplyr::relocate(any_of(c("QUAL", "rid")), .after = last_col())

  return(tcheck)
}

#' Constructs profile plots for salinity, temperature and density.
#'
#' @param x A tibble as exported from either \code{\link{imp_raw_csv}}
#'       or \code{\link{old_data_updatR}}
#' @param file A character string of full file path to raw csv data.
#'
#' @return Three individual profile plots, facetted by transect, for salinity,
#'       temperature and density (UNESCO formula).
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stringr str_detect
profile_plotR <- function(x, file){

  tname <- plot_nameR(file = file)
  fname <- tools::file_path_sans_ext(basename(file))

  facet_df <- dplyr::tibble(site = as.factor(x[["site"]]),
                            tsect = substr(site, 1, 3)) %>%
    unique() %>%
    dplyr::arrange(site)

  grp_facet <- facet_df %>%
    dplyr::left_join(unique(facet_df) %>%
                       dplyr::group_by(tsect) %>%
                       dplyr::mutate(off = 0:(n()-1))) %>%
    dplyr::select(-tsect)

  if(stringr::str_detect(x[["site"]][1], pattern = "c")){
    # caversham
    dat <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(site)) %>%
      dplyr::arrange(site) %>%
      dplyr::mutate(p = 1 + (DEPTH_m * 0.1),
                    depth = -DEPTH_m,
                    dens = marelac::sw_dens(S = SAL_psu, t = T_degC, p = p,
                                            method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(offD = off * 0.3,
                    offC = off * 0.2,
                    offS = off * 0.3,
                    densP = round(dens + offD, 2),
                    cP = round(T_degC + offS, 2),
                    salP = round(SAL_psu + offS, 2),
                    fac_grp = factor(substr(toupper(site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3"))) %>%
      dplyr::arrange(site, depth)

    depth <- c(0, dat[["DEPTH_m"]])

    # find index of deepest sites
    ind <- c(0, which(diff(depth)>0))

    lab_df <- dat %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = depth - 0.1,
                    fac_grp = factor(substr(toupper(site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3")),
                    site = toupper(site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(site)
  } else {
    # guildford
    dat <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(site)) %>%
      dplyr::arrange(site) %>%
      dplyr:: mutate(p = 1 + (DEPTH_m * 0.1),
                     depth = -DEPTH_m,
                     dens = marelac::sw_dens(S = SAL_psu, t = T_degC, p = p,
                                             method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(offD = off * 0.3,
                    offC = off * 0.2,
                    offS = off * 0.3,
                    densP = round(dens + offD, 2),
                    cP = round(T_degC + offS, 2),
                    salP = round(SAL_psu + offS, 2),
                    fac_grp = factor(substr(toupper(site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3"))) %>%
      dplyr::arrange(site, depth)

    depth <- c(0, dat[["DEPTH_m"]])

    # find index of deepest sites
    ind <- c(0, which(diff(depth)>0))

    lab_df <- dat %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = depth - 0.1,
                    fac_grp = factor(substr(toupper(site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3")),
                    site = toupper(site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(site)
  }

  d_plot <- ggplot(dat, aes(x = densP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    facet_grid(rows = vars(fac_grp)) +
    geom_text(data = lab_df, aes(label = site, x = densP, y = y)) +
    labs(title = tname,
         subtitle = expression(Offset~"0.3"*kg*"/"*m^3),
         y = "Depth (m)",
         x = expression(Density~"("*kg*"/"*m^3*")")) +
    theme_bw() +
    theme(legend.position="none")

  c_plot <- ggplot(dat, aes(x = cP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    facet_grid(rows = vars(fac_grp)) +
    geom_text(data = lab_df, aes(label = site, x = cP, y = y)) +
    labs(title = tname,
         subtitle = expression(Offset~"0.2"*degree*C),
         y = "Depth (m)",
         x = expression(Temperature~"("*degree*C*")")) +
    theme_bw() +
    theme(legend.position="none")

  s_plot <- ggplot(dat, aes(x = salP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    facet_grid(rows = vars(fac_grp)) +
    geom_text(data = lab_df, aes(label = site, x = salP, y = y)) +
    labs(title = tname,
         subtitle = "Offset 0.3ppt",
         y = "Depth (m)",
         x = "Salinity (ppt)") +
    theme_bw() +
    theme(legend.position="none")

  ggsave(d_plot, file = paste0("./", fname, "_density_profile_profile.png"))
  ggsave(c_plot, file = paste0("./", fname, "_temperature_profile.png"))
  ggsave(s_plot, file = paste0("./", fname, "_salinity_profile.png"))
}




