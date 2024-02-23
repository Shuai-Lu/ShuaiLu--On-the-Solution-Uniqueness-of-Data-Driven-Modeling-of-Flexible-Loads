# On the Identification of Generalized Flexible Load


## Description
* main.m: start from this file.

* generate_data.m: This script generates (price, power) data by simulating a dynamic energy management system with adjustable loads and virtual batteries. The code calculates the optimal power consumption under different price signals.

* identification.m: This script identifies a physics-based generalized flexible model, including the forward problem modeling, the derivation of the concise form using KKT conditions, and the inverse optimization problem solution.

* plot_fig1.m: This script visualizes the convex hull of $\Gamma$ (i.e., $\rm Conv(\Gamma)$) and the set $\Pi$ under different numbers of samples.

* plot_fig2.m: This script plots the feasible region of the aggregated power of $\rm{Ω}_{\rm phy}$ under different numbers of virtual batteries.

## Authors
- Jiayi Ding, School of Electrical Engineering, Southeast University, Nanjing, China
- Shuai Lu, School of Electrical Engineering, Southeast University, Nanjing, China


## Usage
- Start from main.m


## Dependencies
- MATLAB
- YALMIP toolbox
- Gurobi solver
- Multi-Parametric Toolbox 3 (MPT3)


## Note
- To reproduce the results presented in an associated paper, set `filename` to `'random_price_benchmark.mat'` in "generate_data.m".

