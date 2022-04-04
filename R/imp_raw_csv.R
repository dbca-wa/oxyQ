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
