#' Prepare data to calculate weights.
#'
#' @description
#' This function prepares time-to-event data for weight calculation using
#' external information.
#'
#' @details
#' External information is the population incidence.
#'
#' @export
#' @importFrom stats aggregate end start
#' @importFrom dplyr mutate
#' @importFrom tidyr %>% spread


Prepare_data <- function(dat, population_incidence, breaks){

  ### Input:

#' @param dat Data.frame with one row per individual which at least includes
#' a column **d** with event indicator (1 for event, 0 for censored), a column
#' **y** with event/censoring time.
#' @param population_incidence A vector (in combination with breaks) or
#' a data.frame (columns 1) 'start age group', 2) 'end age group', 3)'S_pop')
#' with population incidence per 100,000 per interval k.
#' @param breaks Cut-points for the (age/time) groups. Only needed when
#' population_incidence is a vector.

  ### Output:

#' @return Data.frame one row per individual and a.o. columns *id* unique ID;
#' *d* non-censoring indicator; **k** interval of (age) group; **S_k**
#' population interval-based proportion of individuals experiencing the
#' event in intervals later than k; **S_k.** sample
#' proportion of individuals experiencing the event in intervals later
#' than k.

# --------------- Load packages.

#  require(survival)
#  require(dplyr)
#  require(tidyr)

  # --------------- Organize external information (population incidence rate).
  # N.B.: there are different options regarding the breaks (1-4).

  if(is.vector(population_incidence)){

    if(is.numeric(breaks)){ # Option 1: two vectors population_incidence and breaks.
      if(length(breaks)!=(length(population_incidence)+1) ) warning("Number of
                                                                  breaks should
                                                                  equal the
                                                                  number of
                                                                  groups plus
                                                                  one.")
      n_agegroups <- length(population_incidence)
      from <- breaks[-(n_agegroups + 1)]  # Remove last
      to <- breaks[-1] # Remove first
      dat_inc <- data.frame(start = from, end = to, incidence_rate =
                              population_incidence/100000)

    } else if(breaks == "5 years"){ # Option 2: 5 years age groups.
      n_agegroups <- length(population_incidence)
      breaks <- seq(0, 5*(n_agegroups), by = 5)
      from <- breaks[-(n_agegroups + 1)]  # Remove last
      to <- breaks[-1] # Remove first
      dat_inc <- data.frame(start = from, end = to, incidence_rate =
                              population_incidence/100000)

    } else if(breaks == "10 years"){ # Option 3: 10 years age groups.
      n_agegroups <- length(population_incidence)
      breaks <- seq(0, 10*(n_agegroups), by = 10)
      from <- breaks[-(n_agegroups + 1)]  # Remove last
      to <- breaks[-1] # Remove first
      dat_inc <- data.frame(start = from, end = to, incidence_rate =
                              population_incidence/100000)
    }

  } else { # Option 4: if population incidence is given as a data.frame.
    dat_inc <- data.frame(start = population_incidence[,1],
                          end = population_incidence[,2],
                          incidence_rate = population_incidence[,3]/100000)
    breaks <- population_incidence[,1]
  }

  dat_inc <- dat_inc %>% mutate( k = paste("(", start, ",", end, "]", sep=""),
                                 t = end - start)
  # t is the width of the age interval, in years.
  # Note that in the paper, this is referred to
  # as (alpha_k - alpha_({k-1}).

  # --------------- Calculate S_k population from population incidence.

  dat_inc$S_k <- c(exp(-(dat_inc$t * dat_inc$incidence_rate)))

  # --------------- Calculate S_k. based on sample data.

  dat$y_cat <- cut(dat$y, breaks = breaks)
  index <- which(is.na(dat$y_cat)==T)
  if(length(index)>0){
    dat<- dat[-index,]}  # Include only those that fall within an age category.

  # Store number of cases and unaffected individuals per age group.
  agegroups.info <- aggregate(id ~ y_cat + d, data = dat, length)

  agegroups.info <- agegroups.info %>% spread(d, -y_cat)
  colnames(agegroups.info) <- c("k","s_k","r_k")

  # Transform NA to 0 (its real meaning).
  agegroups.info$s_k <- ifelse(!is.na(agegroups.info$s_k), agegroups.info$s_k, 0)
  agegroups.info$r_k <- ifelse(!is.na(agegroups.info$r_k), agegroups.info$r_k, 0)

  # Calculate person years per outcome and age group.
  personyears <- aggregate(y ~ y_cat + d, data =
                             dat[,c("d","y_cat", "y"),], sum) %>%
    spread(d, -y_cat)
  colnames(personyears) <- c("k","q_k","p_k")
  agegroups.info <- merge(agegroups.info, personyears, by = "k")  # Merge.

  agegroups.info$q_k <- ifelse(!is.na(agegroups.info$q_k),
                               agegroups.info$q_k, 0)
  agegroups.info$p_k <- ifelse(!is.na(agegroups.info$p_k),
                               agegroups.info$p_k, 0)

  agegroups.info$years.before <- breaks[which(breaks != max(breaks))]

  # Person years among censored individuals.
  q <- agegroups.info$q_k - (agegroups.info$years.before * agegroups.info$s_k)

  # Person years among cases.
  p <- agegroups.info$p_k - (agegroups.info$years.before * agegroups.info$r_k)

  # Years per age group.
  # N.B. in the paper referred to as (alpha_k - alpha_{k-1}).
  t <- agegroups.info$t <- dat_inc$t

  # Number of censored.
  s <- agegroups.info$s_k

  # Number of cases.
  r <- agegroups.info$r_k

  # Age group label.
  k <- agegroups.info$k

  # Obtain incidence rate (mu) based on sample data.
  n_agegroups <- nrow(agegroups.info)
  for (k in 1:n_agegroups){
    if(k!=n_agegroups){calcsum <- sum(agegroups.info$r_k[(k+1):n_agegroups],
                                      agegroups.info$s_k[(k+1):n_agegroups])
    } else {calcsum <- 0}
    agegroups.info$mu_k[k] = agegroups.info$r_k[k]/
      (agegroups.info$p_k[k]+agegroups.info$q_k[k]+agegroups.info$t[k]*calcsum)

  }

  # Recalculate S_k. (in the sample!) based on the intervals.
  S_k. <- c()
  S_k. <- c(exp(-(agegroups.info$mu_k*t)))

  agegroups.info$S_k. <- S_k.
  agegroups.info$S_k <- dat_inc$S_k

  # --------------- Add S_k. to the sample data.

  dat_out <- merge(dat, agegroups.info, by.x = "y_cat", by.y = "k") %>%
    mutate(k = y_cat) %>%
    mutate(population_incidence = mu_k)

  dat_out$years.before <- NULL; dat_out$t <- NULL;
  dat_out$y_cat <- NULL; dat_out$mu_k <- NULL;
  dat_out$s_k <- dat_out$r_k <- dat_out$p_k <- dat_out$q_k <- NULL

  # --------------- Output the data.frame.

  dat_out

}
