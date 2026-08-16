# Documentation for the illustrative datasets shipped with the package.
# All four datasets are SYNTHETIC (see data-raw/make_data.R); they exist only to
# make examples and tests reproducible and must not be treated as field data.

#' Illustrative site-level bushmeat survey data
#'
#' Synthetic paired-transect data contrasting a lightly hunted *reference* site
#' with a *hunted* site, for use with the index functions [ot_pdc()] and [ot_hyco()].
#' **These are simulated values, not measurements from any real study.**
#'
#' @format A data frame with 12 rows and 5 variables:
#' \describe{
#'   \item{site_type}{`"reference"` or `"hunted"`.}
#'   \item{transect}{Transect identifier.}
#'   \item{density}{Wildlife density (individuals per km^2).}
#'   \item{yield_kg}{Harvested biomass (kg).}
#'   \item{effort_days}{Hunting effort (hunter-days).}
#' }
#' @source Simulated.
#' @seealso [ot_pdc()], [ot_hyco()]
"bushmeat_sites"

#' Illustrative individual-level age structure data
#'
#' Synthetic individual-level records of age class at a reference and a hunted
#' site, for use with [ot_asc()]. **Simulated, not field data.**
#'
#' @format A data frame with 240 rows and 2 variables:
#' \describe{
#'   \item{site_type}{`"reference"` or `"hunted"`.}
#'   \item{age_class}{`"juvenile"` or `"adult"`.}
#' }
#' @source Simulated.
#' @seealso [ot_asc()]
"bushmeat_ages"

#' Illustrative duiker demographic parameters
#'
#' Synthetic species-level life-history and offtake parameters for the
#' non-spatial harvest models [ot_pro()], [ot_pbr()], [ot_msy()] and [ot_samse()].
#' Ranges are loosely inspired by the duiker literature but the values are
#' **simulated and illustrative only.**
#'
#' @format A data frame with 3 rows and 8 variables:
#' \describe{
#'   \item{species}{Common name.}
#'   \item{density_k}{Carrying-capacity density (individuals per km^2).}
#'   \item{annual_take}{Observed annual offtake (individuals per km^2 per year).}
#'   \item{lifespan}{Age at last reproduction (years).}
#'   \item{b}{Mean annual female offspring per female.}
#'   \item{a}{Age at first reproduction (years).}
#'   \item{w}{Age at last reproduction (years).}
#'   \item{rmax}{Maximum instantaneous rate of increase (per year).}
#' }
#' @source Simulated.
#' @seealso [ot_pro()], [ot_pbr()], [ot_msy()], [ot_samse()]
"duiker_demography"

#' Illustrative hunting settlements
#'
#' Synthetic settlement locations and hunter populations for the spatial
#' source-sink model [ot_biode()]. Place names are Amazonian but the coordinates
#' and populations are **simulated and illustrative only.**
#'
#' @format A data frame with 4 rows and 4 variables:
#' \describe{
#'   \item{settlement}{Settlement name.}
#'   \item{x_km}{X coordinate (km).}
#'   \item{y_km}{Y coordinate (km).}
#'   \item{hunters}{Number of hunters.}
#' }
#' @source Simulated.
#' @seealso [ot_biode()]
"manu_settlements"
