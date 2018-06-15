#' Juvenile Pacific ocean perch CPUE 

#' This is a data set containing catch per unit of effort for juvenile (< 250 mm) Pacific 
#' ocean perch in the Gulf of Alaska bottom trawl survey from 1996-2009. The data originated
#' from RACEBASE and includes the midpoint of the latitude and longitude of the tow (lat & long)
#' the survey year, the thermocline depth,the bottom water temperature, the bottom water 
#' depth, the seafloor slope at the survey location, the combined catch of structure forming 
#' invertebrates (sponges and corals), the water temperature at the thermocline depth, the 
#' combined catch of shrimp species in the bottom trawl haul and finally catch per unit of 
#' effort for juvenile POP.

#'
#' @format A data frame with 4475 rows and 12 variables:
#' \describe{
#'   \item{hauljoin}{unique haul identifier from RACEBASE}
#'   \item{lat}{midpoint of haul latitude}
#'   \item{long}{midpoint of haul longitude}
#'   \item{year}{year of survey}
#'   \item{tdepth}{depth of thermocline, in meters}
#'   \item{btemp}{bottom water temperature, in degrees C}
#'   \item{bdepth}{bottom depth, in m}
#'   \item{slope}{slope of the seafloor at the station, in %}
#'   \item{inverts}{combined catch of corals and sponges, in log(CPUE)}
#'   \item{shrimp}{combined catch of shrimp species, in kg/ha}
#'   \item{juvenile_POP_CPUE}{catch of juvenile Pacific ocean perch, in kg/ha}
#'   ...
#' }
#' @source {RACEBASE}
"Juvenile_POP_Data"