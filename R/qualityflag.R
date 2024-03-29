#' Quality flagger function
#'
#' \code{flagR} assesses the raw field sample data from a sonde and flags with a
#'       numerical code any suspect samples in an output csv. OK data has a `QUAL`
#'       code value of 10000. Quality issues are represented by codes that are
#'       added to this value. Quality issues are:
#'       \itemize{
#'         \item Sample stability check (depth/Vpos +/- 0.2m variation) - `QUAL`
#'         code 01000
#'         \item Sample site check for out of sequence record - `QUAL` code 00100
#'         \item Sample site check for duplicate record (within 0.1m) - `QUAL` code 00010
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
#' @param file A character string representing the full file path to raw csv data.
#'
#' @return A csv of the original data with quality issues flagged, a shape file
#'       of the same and salinity, temperature and density profile plots.
#'
#' @import readr
#' @import sf
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
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
  
  cat("Performing quality checks...\n")
  x <- depth_checkR(x)
  
  x <- site_checkR(x)
  
  x <- transect_checkR(x)
  
  plotR(x, file)
  
  stub <- basename(file) %>%
    gsub(pattern = ".csv", replacement = "_R1.csv", file)
  
  csv_fp <- paste0("./", stub)
  shp_fp <- gsub(pattern = ".csv", replacement = ".shp", stub)
  
  # write csv R version
  x %>%
    readr::write_csv(file = csv_fp)
  
  # write shape file
  x %>%
    dplyr::mutate(time = as.character(.data$time)) %>%
    sf::st_as_sf(coords = c("easting_m", "northing_m"), crs = 28350) %>%
    sf::st_write(dsn = shp_fp, quiet = TRUE)
}

#' Import raw csv data from sonde
#'
#' @param file A character string of full file path to raw csv data.
#'
#' @return A tibble of the raw data with a `QUAL` variable set to 10000.
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' @import readr
#' @importFrom janitor clean_names
#' @importFrom mgsub mgsub
#' @import sf
#' @import dplyr
imp_raw_csv <- function(file){
  # dlist <- vector(mode = "list", length = length(file))
  for(i in seq_along(file)){
    f <- file[i]
    d <- readr::read_csv(f, locale = readr::locale(encoding="latin1"),
                         show_col_types = FALSE) %>%
      janitor::clean_names()
    # function to split on last '_' in col names and junk sonde id
    fn <- function(x){strsplit(x, "_(?!.*_)", perl=TRUE)[[1]]}
    # rename cols 1st pass
    n <- names(d)
    names(d) <- sapply(lapply(n, FUN = fn), "[[", 1)
    d <- d %>%
      sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      sf::st_transform(crs = 7850) %>%
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
#' @import dplyr
old_data_updatR <- function(x){
  x <- x %>%
    dplyr::mutate(site = ifelse(grepl("c0", .data$site),
                                paste0("cs",
                                       substr(.data$site, 2, 3)),
                                .data$site),
                  site = ifelse(grepl("go", .data$site),
                                paste0("gs0",
                                       substr(.data$site, 4, 4)),
                                .data$site))
}

#' Helper function to reformat a file name.
#'
#' @param file A character string of full file path to raw csv data.
#'
#' @return A more human friendly plot title.
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' @importFrom lubridate as_date
#' @importFrom stringr str_split
#' @importFrom tools file_path_sans_ext
plot_nameR <- function(file){
  date <- lubridate::as_date(substr(basename(file), 1, 8))
  loc <- stringr::str_split(tools::file_path_sans_ext(basename(file)), "_")[[1]][2]
  rev <- stringr::str_split(tools::file_path_sans_ext(basename(file)), "_")[[1]][3]
  if(is.na(rev)){
    sname <- paste(date, loc)
  } else {
    sname <- paste(date, loc, rev)
  }
  
  return(sname)
}

#' Checks stability of sample by comparison of depth to vpos measurements on sonde.
#'
#' @param x A tibble as exported from either \code{\link{imp_raw_csv}}
#'       or \code{\link{old_data_updatR}}
#'
#' @return A tibble where `QUAL` variable has 1000 added if > 0.1m variation.
#'
#' @import dplyr
depth_checkR <- function(x){
  out <- x %>%
    dplyr::mutate(dif = round((.data$DEPTH_m - .data$VPOS_m), 2),
                  QUAL = case_when(
                    abs(dif) >= 0.1  ~ QUAL + 1000,
                    TRUE ~ QUAL)) %>%
    dplyr::select(-dif)
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
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' @import dplyr
site_checkR <- function(x){
  out <- dplyr::tibble()
  sites_to_do <- unique(x[["site"]])
  
  for(i in seq_along(sites_to_do)){
    s = sites_to_do[i]
    d <- x %>% dplyr::filter(.data$site == s) %>% dplyr::pull(.data$DEPTH_m)
    d <- c(d, 0)
    odf <- x %>%
      dplyr::filter(.data$site == s) %>%
      dplyr::mutate(dd = diff(d),
                    dup = abs(-.data$DEPTH_m - lead(-.data$DEPTH_m)),
                    QUAL = ifelse(dd > 0, .data$QUAL + 100, .data$QUAL),
                    QUAL = case_when(
                      dup < 0.1 ~ .data$QUAL + 10,
                      is.na(dup) ~ .data$QUAL,
                      TRUE ~ .data$QUAL)) %>%
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
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' @import sf
#' @import dplyr
transect_checkR <- function(x){
  # turn csv to shp
  x_shp <- sf::st_as_sf(x, coords = c("easting_m", "northing_m"), crs = 7850)
  
  # find outside extents and label tsects
  tcheck <- x_shp %>%
    dplyr::mutate(in_ext = lengths(st_within(x_shp, tpolys))) %>%
    sf::st_join(tpolys, join = st_within) %>%
    dplyr::select(-id, -grp, -buf) %>%
    dplyr::mutate(QUAL = case_when(
      is.na(.data$tsect) ~ .data$QUAL + 1,
      TRUE ~ .data$QUAL)) %>%
    dplyr::mutate(tsamp = substr(tolower(.data$site), 1, 3),
                  QUAL = case_when(
                    .data$tsamp != tolower(.data$tsect) ~ .data$QUAL + 1,
                    TRUE ~ .data$QUAL
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
#' \code{profile_plotR} creates salinity, temperature and density profile plots.
#'       It is used only on user edited quality flagged csv inputs and each
#'       iteration will retain the input 'R' value in the plot file name.
#'
#' @param file A character string. Full file path to raw csv data.
#'
#' @return Three individual profile plots, facetted by transect, for salinity,
#'       temperature and density (UNESCO formula).
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stringr str_detect
#'
#' @export
#'
#' @examples
#' \dontrun{
#' profile_plotR("./my_raw_sonde_data.csv")
#' }
profile_plotR <- function(file){
  
  x <- readr::read_csv(file, show_col_types = FALSE)
  
  tname <- plot_nameR(file = file)
  fname <- tools::file_path_sans_ext(basename(file))
  
  
  facet_df <- dplyr::tibble(site = as.factor(x[["site"]]),
                            tsect = substr(.data$site, 1, 3)) %>%
    unique() %>%
    dplyr::arrange(.data$site)
  
  grp_facet <- facet_df %>%
    dplyr::left_join(unique(facet_df) %>%
                       dplyr::group_by(.data$tsect) %>%
                       dplyr::mutate(off = 0:(n()-1)), by = c("site", "tsect")) %>%
    dplyr::select(-.data$tsect)
  
  if(stringr::str_detect(tolower(x[["site"]][1]), pattern = "c")){
    # caversham
    dat <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(.data$site)) %>%
      dplyr::arrange(.data$site) %>%
      dplyr::mutate(p = 1 + (.data$DEPTH_m * 0.1),
                    depth = -.data$DEPTH_m,
                    dens = marelac::sw_dens(S = .data$SAL_psu, t = .data$T_degC, p = .data$p,
                                            method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(offD = .data$off * 0.3,
                    offC = .data$off * 0.2,
                    offS = .data$off * 0.3,
                    densP = round(.data$dens + .data$offD, 2),
                    cP = round(.data$T_degC + .data$offC, 2),
                    salP = round(.data$SAL_psu + .data$offS, 2),
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3"))) %>%
      dplyr::arrange(.data$site, .data$depth)
    
    # no offset data
    dat_noff <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(.data$site)) %>%
      dplyr::arrange(.data$site) %>%
      dplyr::mutate(p = 1 + (.data$DEPTH_m * 0.1),
                    depth = -.data$DEPTH_m,
                    dens = marelac::sw_dens(S = .data$SAL_psu, t = .data$T_degC, p = .data$p,
                                            method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(densP = round(.data$dens, 2),
                    cP = round(.data$T_degC, 2),
                    salP = round(.data$SAL_psu, 2),
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3"))) %>%
      dplyr::arrange(.data$site, .data$depth)
    
    depth <- c(0, dat[["DEPTH_m"]])
    
    # find index of deepest sites
    ind <- c(0, which(diff(depth)>0))
    
    lab_df <- dat %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = .data$depth - 0.1,
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3")),
                    site = toupper(.data$site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(.data$site)
    
    # no offset labels
    lab_df_noff <- dat_noff %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = .data$depth - 0.1,
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3")),
                    site = toupper(.data$site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(.data$site)
  } else {
    # guildford
    dat <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(.data$site)) %>%
      dplyr::arrange(.data$site) %>%
      dplyr:: mutate(p = 1 + (.data$DEPTH_m * 0.1),
                     depth = -.data$DEPTH_m,
                     dens = marelac::sw_dens(S = .data$SAL_psu, t = .data$T_degC, p = .data$p,
                                             method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(offD = .data$off * 0.3,
                    offC = .data$off * 0.2,
                    offS = .data$off * 0.3,
                    densP = round(.data$dens + .data$offD, 2),
                    cP = round(.data$T_degC + .data$offC, 2),
                    salP = round(.data$SAL_psu + .data$offS, 2),
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3"))) %>%
      dplyr::arrange(.data$site, .data$depth)
    
    # no offset data
    dat_noff <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(.data$site)) %>%
      dplyr::arrange(.data$site) %>%
      dplyr:: mutate(p = 1 + (.data$DEPTH_m * 0.1),
                     depth = -.data$DEPTH_m,
                     dens = marelac::sw_dens(S = .data$SAL_psu, t = .data$T_degC, p = .data$p,
                                             method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(densP = round(.data$dens, 2),
                    cP = round(.data$T_degC, 2),
                    salP = round(.data$SAL_psu, 2),
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3"))) %>%
      dplyr::arrange(.data$site, .data$depth)
    
    
    depth <- c(0, dat[["DEPTH_m"]])
    
    # find index of deepest sites
    ind <- c(0, which(diff(depth)>0))
    
    lab_df <- dat %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = .data$depth - 0.1,
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3")),
                    site = toupper(.data$site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(.data$site)
    
    # no offset labels
    lab_df_noff <- dat_noff %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = .data$depth - 0.1,
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3")),
                    site = toupper(.data$site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(.data$site)
  }
  
  d_plot <- ggplot(dat, aes(x = densP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df, aes(label = site, x = densP, y = y)) +
    labs(title = tname,
         subtitle = expression(Offset~"0.3"*kg*"/"*m^3),
         y = "Depth (m)",
         x = expression(Density~"("*kg*"/"*m^3*")")) +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  # no offset plot
  d_noff <- ggplot(dat_noff, aes(x = densP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df_noff, aes(label = site, x = densP, y = y)) +
    labs(title = tname,
         subtitle = "No offset",
         y = "Depth (m)",
         x = expression(Density~"("*kg*"/"*m^3*")")) +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  c_plot <- ggplot(dat, aes(x = cP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df, aes(label = site, x = cP, y = y)) +
    labs(title = tname,
         subtitle = expression(Offset~"0.2"*degree*C),
         y = "Depth (m)",
         x = expression(Temperature~"("*degree*C*")")) +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  # no offset plot
  c_noff <- ggplot(dat_noff, aes(x = cP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df_noff, aes(label = site, x = cP, y = y)) +
    labs(title = tname,
         subtitle = "No offset",
         y = "Depth (m)",
         x = expression(Temperature~"("*degree*C*")")) +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  s_plot <- ggplot(dat, aes(x = salP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df, aes(label = site, x = salP, y = y)) +
    labs(title = tname,
         subtitle = "Offset 0.3ppt",
         y = "Depth (m)",
         x = "Salinity (ppt)") +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  # no offset plot
  s_noff <- ggplot(dat_noff, aes(x = salP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df_noff, aes(label = site, x = salP, y = y)) +
    labs(title = tname,
         subtitle = "No offset",
         y = "Depth (m)",
         x = "Salinity (ppt)") +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  suppressMessages(ggsave(d_plot, filename = paste0("./", fname, "_density_profile.png")))
  suppressMessages(ggsave(d_noff, filename = paste0("./", fname, "_density_profile_no_offset.png")))
  suppressMessages(ggsave(c_plot, filename = paste0("./", fname, "_temperature_profile.png")))
  suppressMessages(ggsave(c_noff, filename = paste0("./", fname, "_temperature_profile_no_offset.png")))
  suppressMessages(ggsave(s_plot, filename = paste0("./", fname, "_salinity_profile.png")))
  suppressMessages(ggsave(s_noff, filename = paste0("./", fname, "_salinity_profile_no_offset.png")))
}




#' Internal profile plot creator
#'
#' @param x A tibble as exported from either \code{\link{imp_raw_csv}}
#'       or \code{\link{old_data_updatR}}
#'
#' @param file A character string. Full file path to raw csv data.
#'
#' @return Three individual profile plots, facetted by transect, for salinity,
#'       temperature and density (UNESCO formula).
#'
#' @author Bart Huntley, \email{bart.huntley@@dbca.wa.gov.au}
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stringr str_detect
#' @importFrom ggh4x facet_grid2
plotR <- function(x, file){
  
  tname <- plot_nameR(file = file)
  fname <- tools::file_path_sans_ext(basename(file))
  
  
  facet_df <- dplyr::tibble(site = as.factor(x[["site"]]),
                            tsect = substr(.data$site, 1, 3)) %>%
    unique() %>%
    dplyr::arrange(.data$site)
  
  grp_facet <- facet_df %>%
    dplyr::left_join(unique(facet_df) %>%
                       dplyr::group_by(.data$tsect) %>%
                       dplyr::mutate(off = 0:(n()-1)), by = c("site", "tsect")) %>%
    dplyr::select(-tsect)
  
  if(stringr::str_detect(tolower(x[["site"]][1]), pattern = "c")){
    # caversham
    dat <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(.data$site)) %>%
      dplyr::arrange(.data$site) %>%
      dplyr::mutate(p = 1 + (.data$DEPTH_m * 0.1),
                    depth = -.data$DEPTH_m,
                    dens = marelac::sw_dens(S = .data$SAL_psu, t = .data$T_degC, p = .data$p,
                                            method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(offD = .data$off * 0.3,
                    offC = .data$off * 0.2,
                    offS = .data$off * 0.3,
                    densP = round(.data$dens + .data$offD, 2),
                    cP = round(.data$T_degC + .data$offC, 2),
                    salP = round(.data$SAL_psu + .data$offS, 2),
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3"))) %>%
      dplyr::arrange(.data$site, .data$depth)
    
    # no offset data
    dat_noff <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(.data$site)) %>%
      dplyr::arrange(.data$site) %>%
      dplyr::mutate(p = 1 + (.data$DEPTH_m * 0.1),
                    depth = -.data$DEPTH_m,
                    dens = marelac::sw_dens(S = .data$SAL_psu, t = .data$T_degC, p = .data$p,
                                            method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(densP = round(.data$dens, 2),
                    cP = round(.data$T_degC, 2),
                    salP = round(.data$SAL_psu, 2),
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3"))) %>%
      dplyr::arrange(.data$site, .data$depth)
    
    depth <- c(0, dat[["DEPTH_m"]])
    
    # find index of deepest sites
    ind <- c(0, which(diff(depth)>0))
    
    lab_df <- dat %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = .data$depth - 0.1,
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3")),
                    site = toupper(.data$site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(.data$site)
    
    # no offset labels
    lab_df_noff <- dat_noff %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = .data$depth - 0.1,
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("CD3", "CD2", "CD1", "CS0", "CU1",
                                                "CU2", "CU3")),
                    site = toupper(.data$site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(.data$site)
  } else {
    # guildford
    dat <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(.data$site)) %>%
      dplyr::arrange(.data$site) %>%
      dplyr:: mutate(p = 1 + (.data$DEPTH_m * 0.1),
                     depth = -.data$DEPTH_m,
                     dens = marelac::sw_dens(S = .data$SAL_psu, t = .data$T_degC, p = .data$p,
                                             method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(offD = .data$off * 0.3,
                    offC = .data$off * 0.2,
                    offS = .data$off * 0.3,
                    densP = round(.data$dens + .data$offD, 2),
                    cP = round(.data$T_degC + .data$offC, 2),
                    salP = round(.data$SAL_psu + .data$offS, 2),
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3"))) %>%
      dplyr::arrange(.data$site, .data$depth)
    
    # no offset data
    dat_noff <- x %>%
      dplyr::select(date, site, T_degC, DEPTH_m, SAL_psu) %>%
      dplyr::mutate(site = as.factor(.data$site)) %>%
      dplyr::arrange(.data$site) %>%
      dplyr:: mutate(p = 1 + (.data$DEPTH_m * 0.1),
                     depth = -.data$DEPTH_m,
                     dens = marelac::sw_dens(S = .data$SAL_psu, t = .data$T_degC, p = .data$p,
                                             method = "UNESCO")) %>%
      dplyr::left_join(grp_facet, by = "site") %>%
      dplyr::mutate(densP = round(.data$dens, 2),
                    cP = round(.data$T_degC, 2),
                    salP = round(.data$SAL_psu, 2),
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3"))) %>%
      dplyr::arrange(.data$site, .data$depth)
    
    
    
    depth <- c(0, dat[["DEPTH_m"]])
    
    # find index of deepest sites
    ind <- c(0, which(diff(depth)>0))
    
    lab_df <- dat %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = .data$depth - 0.1,
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3")),
                    site = toupper(.data$site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(.data$site)
    
    lab_df_noff <- dat_noff %>%
      dplyr::slice(ind) %>%
      dplyr::mutate(y = .data$depth - 0.1,
                    fac_grp = factor(substr(toupper(.data$site), 1, 3),
                                     levels = c("GD3", "GD2", "GD1", "GS0", "GU1",
                                                "GU2", "GU3")),
                    site = toupper(.data$site)) %>%
      dplyr::select(site, y, densP, cP, salP, fac_grp) %>%
      dplyr::arrange(.data$site)
  }
  
  d_plot <- ggplot(dat, aes(x = densP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df, aes(label = site, x = densP, y = y)) +
    labs(title = tname,
         subtitle = expression(Offset~"0.3"*kg*"/"*m^3),
         y = "Depth (m)",
         x = expression(Density~"("*kg*"/"*m^3*")")) +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  
  # no offset plot
  d_noff <- ggplot(dat_noff, aes(x = densP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df_noff, aes(label = site, x = densP, y = y)) +
    labs(title = tname,
         subtitle = "No offset",
         y = "Depth (m)",
         x = expression(Density~"("*kg*"/"*m^3*")")) +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  
  c_plot <- ggplot(dat, aes(x = cP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df, aes(label = site, x = cP, y = y)) +
    labs(title = tname,
         subtitle = expression(Offset~"0.2"*degree*C),
         y = "Depth (m)",
         x = expression(Temperature~"("*degree*C*")")) +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  # no offset plot
  c_noff <- ggplot(dat_noff, aes(x = cP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df_noff, aes(label = site, x = cP, y = y)) +
    labs(title = tname,
         subtitle = "No offset",
         y = "Depth (m)",
         x = expression(Temperature~"("*degree*C*")")) +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  s_plot <- ggplot(dat, aes(x = salP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df, aes(label = site, x = salP, y = y)) +
    labs(title = tname,
         subtitle = "Offset 0.3ppt",
         y = "Depth (m)",
         x = "Salinity (ppt)") +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  # no offset plot
  s_noff <- ggplot(dat_noff, aes(x = salP, y = depth, colour = site)) +
    geom_path() +
    geom_point() +
    ggh4x::facet_grid2(rows = vars(fac_grp), scales = "free_x", independent = "x") +
    geom_text(data = lab_df_noff, aes(label = site, x = salP, y = y)) +
    labs(title = tname,
         subtitle = "No offset",
         y = "Depth (m)",
         x = "Salinity (ppt)") +
    lims(y = c(-5, 0)) +
    theme_bw() +
    theme(legend.position="none")
  
  suppressMessages(ggsave(d_plot, filename = paste0("./", fname, "_density_profile.png")))
  suppressMessages(ggsave(d_noff, filename = paste0("./", fname, "_density_profile_no_offset.png")))
  suppressMessages(ggsave(c_plot, filename = paste0("./", fname, "_temperature_profile.png")))
  suppressMessages(ggsave(c_noff, filename = paste0("./", fname, "_temperature_profile_no_offset.png")))
  suppressMessages(ggsave(s_plot, filename = paste0("./", fname, "_salinity_profile.png")))
  suppressMessages(ggsave(s_noff, filename = paste0("./", fname, "_salinity_profile_no_offset.png")))
}

utils::globalVariables(c("densP", "site", "fac_grp", "y", "cP", "salP", "."))
