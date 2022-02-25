# regression_with_spatial_correlation

What happens when you fit a linear regression to spatially correlated data,
and calculate std errors and confidence intervals assuming the data are i.i.d.?

In this simulation, the researcher's goal is to produce unbiased estimates and valid
confidence intervals for beta_hat in a setting where Y = constant + x*beta + epsilon.

Related:
 - https://voxeu.org/article/standard-errors-persistence
 - https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3398303
 - https://peterroelants.github.io/posts/gaussian-process-kernels/
 - https://www.stat.purdue.edu/~boli/stat512/lectures/topic3.pdf

![Epsilon raster high spatial correlation](map_epsilon_single_simulation_length_epsilon_2.png)

![Epsilon raster low spatial correlation](map_epsilon_single_simulation_length_epsilon_0.25.png)

![Covariance function](covariance_function_with_varying_levels_of_spatial_correlation.png)

![Sampling Distribution](sampling_distribution_beta_x1_with_varying_levels_of_spatial_correlation.png)

![CI Coverage](ci_coverage_for_beta_x1_with_varying_levels_of_spatial_correlation.png)
