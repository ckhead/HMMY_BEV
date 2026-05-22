# This code is for single-stage uncapacitated facility location problem using MILP (JuMP-Gurobi) to be compared with Arkolakis, Eckert and Shi (2023) using CDCP
# Vary heterogeneity parameters one at a time and save results

using JuMP
import Gurobi
const ENV = Gurobi.Env()
# Rowan Shi's CDCP.jl package
# import Pkg
# Pkg.add(url="https://github.com/rowanxshi/CDCP.jl#main")
using CDCP
using DataFrames
using IterTools # used for Brute force
using Plots
using JLD2
using Distributed, ThreadsX

# Import auxiliary functions
include("speed_setup_par.jl")
include("speed_functions_par.jl")
CPU = getCPU()

# par = ModelParams(σᵥ=0.25, β=0.25, σᵩ=1.0, μᵩ=-1.0, μₐ = 2.0, num_sims = 100)
# par = ModelParams(σᵥ=0.48, β=0.25, σᵩ=1.0, μᵩ=-2.0, num_sims = 100) 
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 4.5, μₙ = 0.0, num_sims = 300)

L_list = [20, 40, 80] # to be consistent with the grid in speed_AES_MILP_BF.jl
N = 40;

# Warmup once per L (the size of location array to avoid AES Julia JIT compilation)
for L in L_list
    N_warm = 10 # a small N for warmup
    speed_compare(L, N_warm, par, mode = :threads)
end

# Define all parameter ranges
par_ranges = Dict(
    20 => 0.0:0.01:1.0,
    40 => 0.0:0.01:1.0,
    80 => 0.0:0.01:1.0,
)

## Vary location-specific variable cost shocks σᵥ
all_Var_σᵥ = [];
for L in L_list
    # N = L
    println("L =$L, N = $N")
    Var_σᵥ = DataFrame(L = L, N = N, σᵥ = Float64[], t_milp = Float64[], t_aes = Float64[], opt_gap = Float64[], loc_gap = Float64[], sdt_milp = Float64[], sdt_aes = Float64[], avg_N_nodes = Float64[], avg_N_undet = Float64[], avg_N_open = Float64[], avg_N_close = Float64[], avg_N_plants = Float64[])
    for σ in par_ranges[L]
        println(" σᵥ = $σ")
        t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants = speed_compare(L, N, ModelParams(par; σᵥ = σ), mode = :threads)
        push!(Var_σᵥ, [L, N, σ, t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants])
    end
    push!(all_Var_σᵥ, Var_σᵥ)
end
# Save dataframes in a single JLD2 file
dfV = vcat(all_Var_σᵥ...)
file_path = "Data/JLD/MILP_CDCP_VarSigmaV_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
save(file_path, "v", dfV)

## Vary path cost shocks σᵤ
all_Var_σᵤ = [];
@time for L in L_list
    # N = L
    println("L = $L, N = $N")
    Var_σᵤ = DataFrame(L = L, N = N, σᵤ = Float64[], t_milp = Float64[], t_aes = Float64[], opt_gap = Float64[], loc_gap = Float64[], sdt_milp = Float64[], sdt_aes = Float64[], avg_N_nodes = Float64[], avg_N_undet = Float64[], avg_N_open = Float64[], avg_N_close = Float64[], avg_N_plants = Float64[])
    for σ in par_ranges[L]
        println(" σᵤ = $σ")
        t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants = speed_compare(L, N, ModelParams(par; σᵤ = σ), mode = :threads)
        push!(Var_σᵤ, [L, N, σ, t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants])
    end
    push!(all_Var_σᵤ, Var_σᵤ)
end 
# Save dataframes in a single JLD2 file
dfU = vcat(all_Var_σᵤ...)
file_path = "Data/JLD/MILP_CDCP_VarSigmaU_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
save(file_path, "u", dfU)

## Vary trade elasticity β
all_Var_β = [];
for L in L_list[3]
    # N = L
    println("L = $L, N = $N")
    Var_β = DataFrame(L = L, N = N, β = Float64[], t_milp = Float64[], t_aes = Float64[], opt_gap = Float64[], loc_gap = Float64[], sdt_milp = Float64[], sdt_aes = Float64[], avg_N_nodes = Float64[], avg_N_undet = Float64[], avg_N_open = Float64[], avg_N_close = Float64[], avg_N_plants = Float64[])
    # for β in par_ranges[L]
    for β in [0.81, 0.82, 0.83, 0.84, 0.89, 0.9, 0.91]
        println(" β = $β")
        t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants = speed_compare(L, N, ModelParams(par; β = β), mode = :threads)
        push!(Var_β, [L, N, β, t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants])
    end
    push!(all_Var_β, Var_β)
end
# Save dataframes in a single JLD2 file
dfB = vcat(all_Var_β...)
file_path = "Data/JLD/MILP_CDCP_VarBeta_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
save(file_path, "B", dfB)


## Vary location-specific fixed cost shocks σᵩ
# Modify parameter range for σᵩ since we want to focus on the critical range between AES and MILP performances
par_ranges = Dict(
    20 => 0.0:0.02:2.0,
    40 => 0.0:0.02:2.0,
    80 => 0.0:0.02:2.0
)

all_Var_σᵩ = [];
for L in L_list
    # N = L
    println("L = $L, N = $N")
    Var_σᵩ = DataFrame(L = L, N = N, σᵩ = Float64[], t_milp = Float64[], t_aes = Float64[], opt_gap = Float64[], loc_gap = Float64[], sdt_milp = Float64[], sdt_aes = Float64[], avg_N_nodes = Float64[], avg_N_undet = Float64[], avg_N_open = Float64[], avg_N_close = Float64[], avg_N_plants = Float64[])
    for σ in par_ranges[L]
        println(" σᵩ = $σ")
        t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants = speed_compare(L, N, ModelParams(par; σᵩ = σ), mode = :threads)
        push!(Var_σᵩ, [L, N, σ, t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants])
    end
    push!(all_Var_σᵩ, Var_σᵩ)
end
# Save dataframes in a single JLD2 file
dfF = vcat(all_Var_σᵩ...)
file_path = "Data/JLD/MILP_CDCP_VarSigmaF_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
save(file_path, "phi", dfF)



## Vary market demand shocks σₐ
# Modify parameter range for σₐ since we want to focus on the critical range between AES and MILP performances
par_ranges = Dict(
    20 => 0.0:0.01:1.0,
    40 => 0.0:0.01:1.0,
    80 => 0.0:0.01:1.0
)
all_Var_σₐ = [];
for L in L_list
    # N = L
    println("L = $L, N = $N")
    Var_σₐ = DataFrame(L = L, N = N, σₐ = Float64[], t_milp = Float64[], t_aes = Float64[], opt_gap = Float64[], loc_gap = Float64[], sdt_milp = Float64[], sdt_aes = Float64[], avg_N_nodes = Float64[], avg_N_undet = Float64[], avg_N_open = Float64[], avg_N_close = Float64[], avg_N_plants = Float64[])
    for σ in par_ranges[L]
        println(" σₐ = $σ")
        t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants = speed_compare(L, N, ModelParams(par; σₐ = σ), mode = :threads)
        push!(Var_σₐ, [L, N, σ, t_milp, t_aes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants])
    end
    push!(all_Var_σₐ, Var_σₐ)
end
# Save dataframes in a single JLD2 file
dfA = vcat(all_Var_σₐ...)
file_path = "Data/JLD/MILP_CDCP_VarSigmaA_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
save(file_path, "A", dfA)

# σₐ does not affect computational performance of AES or MILP given that market entry is fixed. So, we exclude it in the critical sigma table or Var_plots. We keep the code here for completeness and potential future use.

