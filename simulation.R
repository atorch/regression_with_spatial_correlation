## See https://voxeu.org/article/standard-errors-persistence
## and https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3398303

library(data.table)
library(ggplot2)
library(MASS)

save_maps_single_simulation <- function(locations, distances, length_x1, length_epsilon, beta_x1) {

    covariance_x1 <- gaussian_kernel(distances, length=length_x1)
    covariance_epsilon <- gaussian_kernel(distances, length=length_epsilon)

    locations <- simulate_spatial_data(locations, covariance_x1, covariance_epsilon, beta_x1)

    p <- (ggplot(locations, aes(x=coord_x, y=coord_y, fill=epsilon)) +
          geom_raster() +
          scale_fill_gradient2("epsilon", low="#ef8a62", mid="white", high="#67a9cf", midpoint = 0) +
          xlab("pixel coordinate (easting)") +
          ylab("pixel coordinate (northing)") +
          theme_bw())
    filename <- "map_epsilon_single_simulation.png"  # TODO Put length epsilon in filename
    ggsave(filename=filename, plot=p, width=8, height=6, unit="in")

    p <- (ggplot(locations, aes(x=coord_x, y=coord_y, fill=y)) +
          geom_raster() +
          scale_fill_gradient2(low="#ef8a62", mid="white", high="#67a9cf", midpoint = 0) +
          xlab("pixel coordinate (easting)") +
          ylab("pixel coordinate (northing)") +
          theme_bw())
    filename <- "map_y_single_simulation.png"
    ggsave(filename=filename, plot=p, width=8, height=6, unit="in")
}

simulate_spatial_data <- function(locations, covariance_x1, covariance_epsilon, beta_x1) {
    locations$x1 <- mvrnorm(n=1, mu=rep(0, nrow(locations)), Sigma=covariance_x1)
    locations$epsilon <- mvrnorm(n=1, mu=rep(0, nrow(locations)), Sigma=covariance_epsilon)

    ## True regression function: y = X*beta + epsilon
    ## TODO Include an x2 in the true regression function, just to make it interesting?
    locations$y <- 10 + beta_x1*locations$x1 + locations$epsilon
    return(locations)
}

simulate_spatial_data_and_fit_model <- function(locations, covariance_x1, covariance_epsilon, beta_x1) {

    locations <- simulate_spatial_data(locations, covariance_x1, covariance_epsilon, beta_x1)

    model <- lm(y ~ x1, data=locations)
    summary(model)

    conf_interval <- as.vector(confint(model, "x1"))

    return(c(estimated_beta_x1=as.vector(coefficients(model)["x1"]),
             conf_interval_contains_beta_x1=(conf_interval[1] < beta_x1) & (conf_interval[2] > beta_x1)))

}

run_simulations <- function(locations, distances, length_x1, length_epsilon, beta_x1, n_sims) {

    covariance_x1 <- gaussian_kernel(distances, length=length_x1)
    covariance_epsilon <- gaussian_kernel(distances, length=length_epsilon)

    ## Variance-covariance matrix must be positive (semi-) definite
    all(eigen(covariance_epsilon)$epsilons > 0)

    simulations <- as.data.frame(t(replicate(n_sims, simulate_spatial_data_and_fit_model(locations, covariance_x1, covariance_epsilon, beta_x1))))
    simulations$length_epsilon <- length_epsilon

    return(simulations)
}

gaussian_kernel <- function(distance, sigma=2, length=5) {
    return(sigma^2 * exp(-distance^2 / (2 * length^2)))
}

## Quick and dirty spatial locations (points in a grid)
grid_width <- 20
locations <- expand.grid(coord_x=seq(1, grid_width), coord_y=seq(1, grid_width))

## Sample size for regressions is grid_width^2 = 400 points
dim(locations)
head(locations)

distances <- as.matrix(dist(locations[, c("coord_x", "coord_y")]))

## This is sqrt(2), the spatial distance from locations[1, ] to locations[22, ]
distances[1, 22]

beta_x1 <- 5

## Run a single simulation and visualize it
save_maps_single_simulation(locations, distances, length_x1=2, length_epsilon=2, beta_x1=beta_x1)

lengths_epsilon <- c(1/8, 0.25, 0.5, 1, 2, 4, 6, 8, 10)

simulations <- lapply(lengths_epsilon, function(length_epsilon) {
    run_simulations(locations, distances, length_x1=2, length_epsilon=length_epsilon, beta_x1=beta_x1, n_sims=500)
})

simulations <- do.call(rbind, simulations)

simulations$length_epsilon_subtitle <- sprintf("kernel length for epsilon: %s", simulations$length_epsilon)

p <- (ggplot(simulations, aes(x=estimated_beta_x1)) +
      geom_histogram(binwidth=0.1, color="black", fill="white") +
      geom_vline(xintercept=beta_x1, linetype=2, alpha=0.5) +
      facet_wrap(~ length_epsilon_subtitle) +
      ylab("") +
      theme_bw())
filename <- "sampling_distribution_beta_x1_with_varying_levels_of_spatial_correlation.png"
ggsave(filename=filename, plot=p, width=8, height=6, unit="in")

## Looks like estimated_beta_x1 is unbiased,
## but the 95% CI does not have anything close to 95% coverage
## when epsilon has a high degree of spatial correlation
summary(simulations)

simulations <- data.table(simulations)

coverage <- simulations[, list(pr_conf_interval_contains_beta_x1=mean(conf_interval_contains_beta_x1)), by=c("length_epsilon")]

p <- (ggplot(coverage, aes(x=length_epsilon, y=pr_conf_interval_contains_beta_x1)) +
      geom_point() +
      geom_hline(yintercept=0.95, linetype=2, alpha=0.5) +
      ylab("coverage of 95% confidence interval for beta") +
      xlab("kernel length for epsilon (controls spatial correlation)") +
      scale_y_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
      theme_bw())
filename <- "ci_coverage_for_beta_x1_with_varying_levels_of_spatial_correlation.png"
ggsave(filename=filename, plot=p, width=8, height=6, unit="in")

## TODO Test for spatial autocorrelation in errors

## TODO Bootstrap with spatial correlation?

## TODO Estimate kernel length by maximum likelihood?

## TODO Plot semi variogram (both of Y and of residuals!)
## What matters (for messing up CI coverage) is spatial correlation in residuals, not in X
