# Multi-stage UFLP and plot

# using LinearAlgebra
using JuMP
import Gurobi
const ENV = Gurobi.Env()
using DataFrames
using Plots
using JLD2

# Import auxiliary functions
include("speed_setup_par.jl")
include("speed_functions_par.jl")
include("flatworld_plot.jl")

function output_muflp(s, L, N, par::ModelParams, directory, filename)
    Π, config, active_y, Xc, Yc, Xf, Yf, fc, A, _= solve_MILP_1sim(s, L, N, par);
    p = plot_solution(L, N, config, active_y, Xc, Yc, Xf, Yf, fc, A)
    savefig(p, directory * filename)
end

directory = "Plots/Others/"
isdir(directory) || mkpath(directory) # create the output directory if it doesn't exist already
s = 0; # seed

## K = 1 ad-hoc cases for plotting
# L = 100; N = 50; K = length(L);
# par = ModelParams(σᵥ=0.48, β=0.25, σᵩ=1.0, μᵩ=-2.0) 
# filename = "flatworld_speed_AES_MILP.pdf"
# output_muflp(34, L, N, par, directory, filename);

L = 40; N = 40; K = length(L);
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 4.5, μₙ = 0.0)
filename = "flatworld_speed_AES_MILP_sigmaV$(par.σᵥ).pdf"
output_muflp(s, L, N, par, directory, filename);

L = 40; N = 40; K = length(L);
par = ModelParams(σᵥ = 0.2, β = 0.25, σᵩ = 1.5, μₐ = 4.5, μₙ = 0.0)
filename = "flatworld_speed_AES_MILP_sigmaV$(par.σᵥ).pdf"
output_muflp(s, L, N, par, directory, filename);


##  K = 1
L = 30; N = 40; K = length(L);

# Case I: Low distance elasticity, no fixed cost of market entry
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 4.5, μₙ = 0.0) # Low demand than K=2,3 cases # -2.0
filename = "flatworld_K$(K)_Beta$(par.β)_fs0.pdf"
output_muflp(s, L, N, par, directory, filename);

# Case II: High distance elasticity, no fixed cost of market entry
par = ModelParams(β = 0.75, σᵩ = 1.5, μₐ = 4.5, μₙ = 0.0) # Low demand than K=2,3 cases
filename = "flatworld_K$(K)_Beta$(par.β)_fs0.pdf"
output_muflp(s, L, N, par, directory, filename);

# # Case III: Low distance elasticity, with unit fixed cost of market entry
# par = ModelParams(β = 0.25, σᵩ = 1.0, μₐ = 5.0) 
# filename = "flatworld_K$(K)_Beta$(par.β)_A$(par.μₐ).pdf"
# output_muflp(s, L, N, par, directory, filename);

# # Case IV: High distance elasticity, with unit fixed cost of market entry
# par = ModelParams(β = 0.75, σᵩ = 1.0, μₐ = 5.0) 
# filename = "flatworld_K$(K)_Beta$(par.β)_A$(par.μₐ).pdf"
# output_muflp(s, L, N, par, directory, filename);

# Case III: Low distance elasticity, with unit fixed cost of market entry
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 9.5) 
filename = "flatworld_K$(K)_Beta$(par.β).pdf"
output_muflp(s, L, N, par, directory, filename);

# Case IV: High distance elasticity, with unit fixed cost of market entry
par = ModelParams(β = 0.75, σᵩ = 1.5, μₐ = 9.5) 
filename = "flatworld_K$(K)_Beta$(par.β).pdf"
output_muflp(s, L, N, par, directory, filename);


## K = 2
L = [20 30]; N = 40; K = length(L); 

# Case I: Low distance elasticity
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 9.5) 
filename = "flatworld_K$(K)_Beta$(par.β).pdf"
output_muflp(s, L, N, par, directory, filename);

# Case II: High distance elasticity
par = ModelParams(β = 0.75, σᵩ = 1.5, μₐ = 9.5) 
filename = "flatworld_K$(K)_Beta$(par.β).pdf"
output_muflp(s, L, N, par, directory, filename);

## K = 3
L = [10 20 30]; N = 40; K = length(L); 

# Case I: Low distance elasticity
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 9.5) 
filename = "flatworld_K$(K)_Beta$(par.β).pdf"
output_muflp(s, L, N, par, directory, filename);

# Case II: High distance elasticity
par = ModelParams(β = 0.75, σᵩ = 1.5, μₐ = 9.5) 
filename = "flatworld_K$(K)_Beta$(par.β).pdf"
output_muflp(s, L, N, par, directory, filename);