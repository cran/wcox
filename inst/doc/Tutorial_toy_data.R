## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

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

  agegroups.info <- agegroups.info %>% tidyr::spread(d, -y_cat)
  colnames(agegroups.info) <- c("k","s_k","r_k")

  # Transform NA to 0 (its real meaning).
  agegroups.info$s_k <- ifelse(!is.na(agegroups.info$s_k), agegroups.info$s_k, 0)
  agegroups.info$r_k <- ifelse(!is.na(agegroups.info$r_k), agegroups.info$r_k, 0)

  # Calculate person years per outcome and age group.
  personyears <- aggregate(y ~ y_cat + d, data =
                             dat[,c("d","y_cat", "y"),], sum) %>%
    tidyr::spread(d, -y_cat)
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
#' Calculate inverse probability of selection weights.
#'
#' @description
#' This function calculates weights to correct for ascertainment bias in
#' time-to-event data where clusters are outcome-dependently sampled,
#' for example high-risk families in genetic epidemiological studies in
#' cancer research.
#'
#' @details
#' Weights are based on a comparison between the survival between sample and
#' population. Therefore, besides the sample data, the population incidence rate
#' (per 100 000) is needed as input, as well as the cut-offs of the
#' (age/time-to-event) groups for which this is available. The function provides
#' two options for the latter: cut-offs can be provided manually or using the
#' standard 5- or 10-years (age) categories (0-4, 5-9, ... or 0-9, 10-14, ...).
#' Note that resulting intervals are of the form [xx, xx).
#'
#' @export


Calculate_weights <- function(dat){

### Input:

#' @param dat Data.frame with one row per individual with columns *d*
#' non-censoring indicator; **k** interval of (age) group; **S_k**
#' population interval-based proportion of individuals experiencing the
#' event in intervals later than k; **S_k.** sample
#' proportion of individuals experiencing the event in intervals later
#' than k.

### Output:

#' @return Vector with weights.

  # --------------- Extract variables from input data.frame.

  # Group/interval.
  k <- dat$k

  # Population proportion of individuals experiencing the event in intervals
  # later than k.
  S_k <- dat$S_k

  # Sample proportion of individuals experiencing the event in intervals
  # later than k.
  S_k. <- dat$S_k.

  # --------------- Create empty containers.

  v <- w <- rep(NA, nrow(dat))

  # --------------- Calculate the weights.

  for (n in 1:nrow(dat)){

    # Weights for unaffected.
    v[n] = 1

    # Weights for cases.
    w[n] = ( (1 - S_k[n]) / S_k[n] ) * (S_k.[n] / ( 1 - S_k.[n]))
  }

  merged <- cbind(dat, w, v)
  merged$weight <- NA

  # --------------- Assign w for uncensored, v for censored (interval-wise).

  merged$weight[which(merged$d == 1)] <- merged$w[which(merged$d == 1)]
  merged$weight[which(merged$d == 0)] <- merged$v[which(merged$d == 0)]

  merged$w <- merged$v <- NULL  # Remove old variable.

  # --------------- Collect output.

  vec_weights <- merged$weight

  # If weight is 0, add very small value (coxph does not accept weights of 0).
  vec_weights[which(vec_weights==0)] <- 0.0000001

  # --------------- Print warning if weights are invalid.

  ifelse((sum(vec_weights<0)>0), print("Invalid (negative) weights!"),
         print("No negative weights"))

  # --------------- Return output.

  vec_weights

}







## ----Loading, message = F-----------------------------------------------------

library(wcox)

require(dplyr)
require(tidyr)
require(survival)

## ----Load data set, message = F-----------------------------------------------
# Load toy data.
data("fam_dat")

## ----Show first lines data----------------------------------------------------
# Show the first few lines of the toy data.
head(fam_dat)

## ----Distribution of family size----------------------------------------------
# Show the number of families and top rows.
cat( "There are", unique(fam_dat$family_id) %>% length() , "families in data set. ")

# Examine family sizes.
fam_sizes <- fam_dat %>% group_by(family_id) %>% count()

cat( "Family size is on average ", 
     fam_sizes$n %>% mean() %>% round(1), 
     " and ranges from ", 
     fam_sizes$n %>% min(), " to ", 
     fam_sizes$n %>% max(), ".", sep = "")

## ----External information-----------------------------------------------------

# Enter the incidence by age group per 100 000, starting with the youngest age group.
incidence_rate <- c(2882, 1766, 1367, 1155, 987, 845, 775, 798, 636, 650)

# Define the age group cutoffs.
breaks_manually <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
breaks_10yrs <- "10 years"

## ----Rename variables---------------------------------------------------------
# Rename variables.
my_dat <- rename(fam_dat, id = individual_id, d = event_indicator, y = age)

## ----Prepare data, message = F------------------------------------------------

# Using option 1 (manual cut-offs):
my_dat_out <- Prepare_data(dat = my_dat, population_incidence = incidence_rate, 
                           breaks = breaks_manually)

# Unhash to use option 2 (pre-set cut-offs):
# my_dat_out <- Prepare_data(dat = my_dat, population_incidence = incidence_rate, 
# breaks = "10 years")

## ----Look at prepared data----------------------------------------------------

# Select the newly add columns.
my_dat_out %>% select(id, k, d, S_k, S_k.) %>% arrange(id) %>% head(8)

## ----Calculate weights--------------------------------------------------------
# Calculate weights using the object that is output from the Prepare_data() function.

w <- Calculate_weights(dat = my_dat_out)

## ----Show weights-------------------------------------------------------------
# Show the weights.
my_dat_out %>% mutate(weight = w) %>% 
  select(id, d, weight) %>% 
  arrange(id) %>% 
  filter(id %in% c(1,3))

## ----Fit model using weights, message = F-------------------------------------
# Fit the model.
fit <- coxph(Surv(y, d==1) ~ x + cluster(family_id), weights = w, data =  my_dat_out)

fit

## ----Examine estimates of fitted model----------------------------------------
# Extract estimates.
fit.summary <- summary(fit)

# Summarize findings.
cat("Covariate effect (95% CI): ",
exp(fit.summary$coefficients[1]), " (",
exp(fit.summary$coefficients[1] - 1.96*fit.summary$coefficients[4]), ", ",
exp(fit.summary$coefficients[1] + 1.96*fit.summary$coefficients[4]), ").", sep = "")

