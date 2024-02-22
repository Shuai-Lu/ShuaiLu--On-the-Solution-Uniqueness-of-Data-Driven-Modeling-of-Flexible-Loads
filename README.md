# On the Identification of Generalized Flexible Load


## Description
* main.m: start from this file.

* generate_data.m: This script generates (price, power) data by simulating a dynamic energy management system with adjustable loads and virtual batteries. The code calculates the optimal power consumption under different price signals.

* identification.m: This script identifies a physics-based generalized flexible model, including the forward problem modeling, the derivation of the concise form using KKT conditions, and the inverse optimization problem solution.

* plot_fig1.m: This script visualizes the convex hull of Γ (Conv(Γ)) and the set D under different numbers of samples.

* plot_fig2.m: This script plots the feasible region of the aggregated power of $Ω_{phy}$ under different numbers of virtual batteries.

## Authors
- Jiayi Ding, School of Electrical Engineering, Southeast University, Nanjing, China
- Shuai Lu, School of Electrical Engineering, Southeast University, Nanjing, China


## Usage
- Start from main.m


## Dependencies
- MATLAB
- YALMIP toolbox
- Gurobi solver
- Analyze N-dimensional Convex Polyhedra (Version 1.9.0.2 by Matt J)


## Note
- To reproduce the results presented in an associated paper, set `filename` to `'random_price_benchmark.mat'` in "generate_data.m".

