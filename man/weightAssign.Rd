% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/instrinsicIsotopes.R
\name{weightAssign}
\alias{weightAssign}
\title{Calculate Weights for Isotope Assignments
weightAssign}
\usage{
weightAssign(
  knownLocs,
  isovalues,
  isoSTD,
  intercept,
  slope,
  odds = 0.67,
  relAbund,
  weightRange = c(-1, 1),
  sppShapefile = NULL,
  assignExtent = c(-179, -60, 15, 89),
  element = "Hydrogen",
  surface = FALSE,
  period = "Annual",
  verbose = 1,
  mapDirectory = NULL
)
}
\arguments{
\item{knownLocs}{matrix of capture locations of the same length as
\code{isovalues}}

\item{isovalues}{vector of tissue isotope values from known locations}

\item{isoSTD}{standard deviation from calibration}

\item{intercept}{intercept value from calibration}

\item{slope}{value from calibration}

\item{odds}{odds ratio to use to set likely and unlikely locations defaults
to 0.67}

\item{relAbund}{raster layer of relative abundance that sums to 1.}

\item{weightRange}{vector of length 2 within minimum and maximum values to
weight isotope and relative abundance. Default = c(-1,1)}

\item{sppShapefile}{A polygon spatial layer (sf - MULTIPOLYGON) defining
species range. Assignments are restricted to these areas.}

\item{assignExtent}{definition for the extent of the assignment. Can be used
in place of \code{sppShapefile} to limit assignment. Input should
follow \code{c(xmin,xmax,ymin,ymax)} in degrees longitude and latitude.}

\item{element}{The elemental isotope of interest. Currently the only
elements that are implemented are 'Hydrogen' (default) and 'Oxygen'}

\item{surface}{DEPRECATED function no longer returns surface water values.
Default is 'FALSE' which returns the precipitation isotopes ratio.}

\item{period}{The time period of interest. If 'Annual' returns a raster
of mean annual values in precipitation for the \code{element}. If
'GrowingSeason' returns growing season values in precipitation for
\code{element} of interest.}

\item{verbose}{takes values 0 or 1 (default). 0 prints no output during
run. 1 prints a message detailing where in the process the function
is.}

\item{mapDirectory}{Directory to save/read isotope map from. Can use relative
or absolute addressing. The default value (NULL) downloads to a temporary
directory, so we strongly recommend changing this from the default unless
you're sure you're not going to need these data more than once.}
}
\value{
returns an \code{weightAssign} object containing the following:
  \describe{
   \item{\code{top}}{data.frame with the optimal weightings}
   \item{\code{frontier}}{data.frame with values that fall along the Pareto
     frontier}
   \item{\code{performance}}{data.frame with error rate and assignment area
     for each weight combination}
}
}
\description{
The primary purpose of this function is to determine whether weighting likelihood based isotope assignments
and prior information, such as relative abundance can improve the model performance compared to the
isotope-only model. To do this, we raise the likelihood and prior values to powers from 0.1
to 10 and measure model performance using the assignment error rate and assignment area. Weights < 1 flatten
the likelihood/prior distributions (giving relatively more weight to smaller values) and weights > 1
sharpen the distributions (giving relatively less weight to smaller values. The \code{weightAssign} function
generates origin assignments using stable-hydrogen isotopes in tissue. If first generates
a probability surface of origin assignment from a vector of stable-isotope values for each animal/sample
captured at a known location. Probabilistic assignments are constructed by first converting observed
stable-isotope ratios (isoscape) in either precipitation or surface waters into a 'tissuescape' using
a user-provided intercept, slope and standard deviation. See
\href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0035137}{Hobson et. al. (2012)}.
}
\examples{
\donttest{
extensions <- c("shp", "shx", "dbf", "sbn", "sbx")
tmp <- tempdir()
for (ext in extensions) {
download.file(paste0(
              "https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity",
                     "/master/data-raw/Spatial_Layers/OVENdist.",
                     ext),
              destfile = paste0(tmp, "/OVENdist.", ext), mode = "wb")
}
OVENdist <- sf::st_read(paste0(tmp, "/OVENdist.shp"))
OVENdist <- OVENdist[OVENdist$ORIGIN==2,] # only breeding
sf::st_crs(OVENdist) <- sf::st_crs(4326)

download.file(paste0("https://raw.githubusercontent.com/SMBC-NZP/MigConnectivity",
                     "/master/data-raw/deltaDvalues.csv"),
              destfile = paste0(tmp, "/deltaDvalues.csv"))
OVENvals <- read.csv(paste0(tmp, "/deltaDvalues.csv"))

HBEFbirds <- OVENvals[grep("NH",OVENvals[,1]),]

# Create a spatial object of known capture sites
knownLocs <- sf::st_as_sf(data.frame(Long = rep(-73,nrow(HBEFbirds)),
                                    Lat = rep(43,nrow(HBEFbirds))),
                         coords = c("Long","Lat"),
                         crs = 4326)

#Get OVEN abundance from BBS estimates and read into R #
utils::download.file("https://www.mbr-pwrc.usgs.gov/bbs/ra15/ra06740.zip",
                     destfile = paste0(tmp, "/oven.zip"))
utils::unzip(paste0(tmp, "/oven.zip"), exdir = tmp)
oven_dist <- sf::st_read(paste0(tmp, "/ra06740.shp"))

# Empty raster with the same dimensions as isoscape and Ovenbird distribution

# We do this manually here but the weightedAssign function has been updated
# to ensure the isoscape and abundance rasts have the same extent using
# resampling to match  relAbund to the isoscape.
r <- terra::rast(nrow = 331, ncol = 870,
                 res = c(0.0833333, 0.0833333),
                 xmin = -125.1667, xmax = -52.66672,
                 ymin = 33.49995, ymax = 61.08327,
                 crs = sf::st_crs(4326)$wkt)

# rasterize the polygons from BBS - this is not needed if working with a
# rasterized surface
relativeAbun<-terra::rasterize(terra::vect(sf::st_transform(oven_dist,4326)),
                               r,
                               field = "RASTAT")

relativeAbund <- relativeAbun/terra::global(relativeAbun, sum,
                                            na.rm = TRUE)$sum


BE <- weightAssign(knownLocs = knownLocs,
                   isovalues = HBEFbirds[,2],
                   isoSTD = 12,
                   intercept = -10,
                   slope = 0.8,
                   odds = 0.67,
                   relAbund = relativeAbund,
                   weightRange = c(-1, 1),
                   sppShapefile = OVENdist,
                   assignExtent = c(-179,-60,15,89),
                   element = "Hydrogen",
                   period = "Annual")
}

}
\references{
Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
M. Gutierrez Ramirez, and P. P. Marra. 2019. The strength of migratory
connectivity for birds en route to breeding through the Gulf of Mexico.
Ecography 42: 658-669.

Rushing, C. S., P. P. Marra and C. E. Studds. 2017. Incorporating breeding
abundance into spatial assignments on continuous surfaces. Ecology and
Evolution 3: 3847-3855. \doi{10.1002/ece3.2605}

Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
M. Gutierrez Ramirez, and P. P. Marra. 2019. The strength of
migratory connectivity for birds en route to breeding through the Gulf of
Mexico. Ecography 42: 658-669.

Hobson, K. A., S. L. Van Wilgenburg, L. I. Wassenaar, and K. Larson. 2012.
Linking hydrogen isotopes in feathers and precipitation: sources of
variance and consequences for assignment to isoscapes. PLoS ONE 7: e35137.

Rushing, C. S., P. P. Marra, and C. E. Studds. 2017. Incorporating breeding
abundance into spatial assignments on continuous surfaces. Ecology and
Evolution 7: 3847-3855.
}
