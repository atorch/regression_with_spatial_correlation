# regression_with_spatial_correlation

What happens when you fit a linear regression to spatially correlated data,
and calculate std errors and confidence intervals assuming the data are i.i.d.?

In this simulation, the researcher's goal is to produce unbiased estimates and valid
confidence intervals for beta_hat in a setting where Y = constant + x*beta + epsilon.

Related:
 - https://voxeu.org/article/standard-errors-persistence
 - https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3398303

![Epsilon map / raster](map_epsilon_single_simulation_length_epsilon_2.png)

![Sampling Distribution](sampling_distribution_beta_x1_with_varying_levels_of_spatial_correlation.png)

![CI Coverage](ci_coverage_for_beta_x1_with_varying_levels_of_spatial_correlation.png)
