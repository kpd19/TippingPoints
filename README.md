# Summary

This library contains code used in the manuscript titled Three ways to lose a friend: bifurcation-, noise-, and rate-induced tipping in a competitor-facilitator model. This paper analyzes the equilibrium behavior and uses stochastic models to quantify the risk of n-tipping an r-tipping under certain conditions. The modeling and analysis are implemented in R and Julia. 

## Requirements and Setup

The code was built using R version 4.3.2 and Julia version 1.12.5. R can be downloaded [here](https://www.r-project.org). Python can be installed [here](https://julialang.org/downloads/), and can be operated using the julialang extension in Visual Studio Code Editor, which can be downloaded [here](https://code.visualstudio.com/download). The code requires several packages that are not part of the base R installations. After installing R, navigate to the main repository directory and run the `installation.R` script to install necessary packages. The Julia packages used in this research are `Distributions`, `Random`, `KernelDensity`, `LinearAlgebra`, `Dates`, `CSV`, and `Tables`. 

## Stochastic simulations

The scripts in the `stochastic` directory simulate the model using two related methods, the Gillespie Method and the tau-leaping method. The script `CF_Gillespie_single.jl` simulates, plots, an saves a single stochastic realization using the Gillespie algorithm. The script `CF_tauleaping.jl` simulates, plots, and saves a single stochastic realization using the tau-leaping algorithm. The scripts `CF_Gillespie_kde.jl` and `CF_tauleaping_kde.jl` simulates the the model, and calculates the 2-D kernel density across simulations where neither species goes extinct. The script `CF_tauleaping_noise.jl` simulates the the model for varying correlated environmental noise, and records whether one or more species goes extinct during the simulation. 

## Analysis

The scripts in the `analyze` directory uses the output from the simulations to create the graphs in the manuscript and calculate the risk of r-tipping and n-tipping for the different cases. 
