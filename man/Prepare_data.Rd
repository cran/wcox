\name{Prepare_data}
\alias{Prepare_data}
\title{Prepare data for weight calculation using the population incidence.}
\usage{
Prepare_data(dat, population_incidence, breaks)
}
\arguments{
\item{dat}{Data.frame with one row per individual which at least includes a column \strong{d} with event indicator (1 for event, 0 for censored), a column \strong{y} with event/censoring time.}
\item{population_incidence}{A vector (in combination with breaks) or a data.frame (columns 1) 'start age group', 2) 'end age group', 3)'S_pop') with population incidence per 100,000 per interval k.}
\item{breaks}{Cut-points for the (age/time) groups. Only needed when population_incidence is a vector.}
}
\value{
Data.frame ready for weight calculation using function 'Calculate_weights()': one row per individual and a.o. columns \strong{id} unique ID; \strong{d} non-censoring indicator; \strong{k} interval of (age) group; \strong{S_k} population interval-based proportion of individuals experiencing the event in intervals later than k; \strong{S_k.} sample proportion of individuals experiencing the event in intervals later than k.
}
\description{
This function prepares the sample data for weight calculation using external information, i.e. the incidence in the population.
}
\details{
Weights are based on a comparison between the survival between sample and
population. Therefore, besides the sample data, the population incidence rate
(per 100 000) is needed as input, as well as the cut-offs of the
(age/time-to-event) groups for which this is available, unless this is provided in a data.frame.
}
