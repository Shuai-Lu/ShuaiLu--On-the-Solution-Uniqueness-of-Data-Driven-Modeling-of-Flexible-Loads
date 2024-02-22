clc; clear all; close all;
yalmip('clear');

%% generate (price, power) data
generate_data();

%% identify the Physics-Based Model
identification();

%% plot Conv(Gamma) and D under different number of samples
plot_fig1();

%% plot Ωphy under different numbers of virtual battery
plot_fig2();