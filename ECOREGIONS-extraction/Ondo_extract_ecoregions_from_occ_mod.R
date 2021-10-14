require(sf)
require(doParallel)
require(foreach)

#' Parallel version of Geometric binary predicates from the package `sf`
#'
#' Compute geometric binary predicates of intersection between 2 simple features geometry sets in parallel
#'
#' @param x An object of class `sf`, `sfc` or `sfg`
#' @param y An object of class `sf`, `sfc` or `sfg`
#' @param nchunks A numeric value specifying the number of subgeometries the x object should be split into.
#' @param ncores A numeric integer specifying the number of cpu cores to use for parallel processing.
#' @return A sparse list, with list element i an integer vector with all indices j for which intersection(x[i],y[i]) is `TRUE` if they intersect, hence integer(0) if none of them is `TRUE`.
#' @export
st_par_intersects <- function(x, y, nchunks=10, ncores=getOption("mc.cores")){

  if(missing(x))
    stop("x is missing.")

  if(missing(y))
    stop("y is missing.")

  if(!inherits(x, c("sf","sfc","sfg")))
    stop("x must be an sf or sfc object.")

  if(!inherits(y, c("sf","sfc","sfg")))
    stop("y must be an sf or sfc object.")

  #--------------------------------------------------
  #= 1. Split the geometric features in nchunks
  #--------------------------------------------------
  ngeom <- length(sf::st_geometry(x))
  iterations <- chunk2(1:ngeom, nchunks)

  #--------------------------------------------------
  #= 2. Apply st_intersects to each chunk in parallel
  #--------------------------------------------------
  on.exit({
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()}
  )
  ncores = if(is.null(ncores)) parallel::detectCores()-1 else min( min(2, ncores), parallel::detectCores()) # number of cores to use
  # update maximum allowed number of additional R processes
  options(mc.cores=ncores)

  if(ncores < 2)
    foreach::registerDoSEQ()
  else
    doParallel::registerDoParallel(ncores)

  out <- vector(length=length(iterations))

  out <- tryCatch({
    foreach::foreach(j=1:length(iterations), .packages="sf", .combine='c') %dopar% {
      sf::st_intersects(x[iterations[[j]], ], y)
    }
  },error=function(err){
    foreach::foreach(j=1:length(iterations), .combine='c') %do% {
      sf::st_intersects(x[iterations[[j]], ], y)
    }
  }, finally = {
    ## Stop the cluster
    doParallel::stopImplicitCluster()
    unregister()
  })

  out
}

#' helper function that split a vector x in n chunks of approximatively equal size
#' @export
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

# path to csv file containing species occurrence records (can be another format)
loc_data = "palm_occurences.txt"

# read ecoregions data (replace by the path to your shapefile)
ecoregions <- sf::st_read("Ecoregions2017.shp", quiet=TRUE)

# read input occurrence data
loc_dat <- try(read.table(loc_data,h=TRUE,sep = "\t"),silent=TRUE) # replace read.csv by read.table for non csv file.
if(inherits(loc_dat,"try-error"))
  stop("Unable to read the csv file of the input dataset.")
sf_loc_data <- loc_dat %>%
  as.data.frame(row.names=NULL)

# retrieve coordinates headers  
id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(sf_loc_data), value=TRUE)[1]
id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(sf_loc_data), value=TRUE)[1]
coordHeaders <- c(id_x_lon, id_y_lat)

# make a loop on all species by splitting the data (if a lot of data, should better use original script and loop through 
# input files, with a file per species)

# split data
sf_loc_data_sub <- split(sf_loc_data, sf_loc_data$accepted_name) 
#create table to store results
results <- data.frame(names(sf_loc_data_sub))
results$ECO <- rep(NA, length(results[,1]))
names(results)[1] <- "accepted_name"

# loop

for (s in 1:length(names(sf_loc_data_sub))) {
  
  sf_loc_data_spe <- sf_loc_data_sub[[s]]

  # convert to sf object
  sf_loc_data_spe <- sf_loc_data_spe %>%
  dplyr::rename(X=coordHeaders[1], Y=coordHeaders[2]) %>%
  dplyr::select(c(X,Y)) %>%
  sf::st_as_sf(., coords=c("X","Y"), crs=sf::st_crs(4326))  

num.point.by.ecoregions <- suppressMessages(lengths(st_par_intersects(ecoregions, sf_loc_data_spe, nchunks=100)))

num.ecoregions.occupied <- sum(num.point.by.ecoregions > 0)

results[s,2] <- num.ecoregions.occupied

# check if/which ecoregions are occupied (i.e. with at least one occurrence records inside)
#ecoregions.with.points <- ecoregions[num.point.by.ecoregions > 0,]
#ecoregions.with.points

}

write.table(results, "NumEcoregionsPerSpecies.txt", sep="\t", row.names = F)