# This file is the master file that runs all the Julia files for the methodology and computation paper
# IMPORTANT!! press Ctrl+D or trashcan after every code cell (defined by ##) to clear workspace
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
# Pkg.instantiate()    # install anything missing
# Pkg.status()         # sanity check

################################################################################
## Speed tests
################################################################################
"""
Speed test for AES vs. MILP vs. Brute Force with varying number of locations L and consumers N
     - single-stage, single-product problem
     - serve all or none market
     - fix one set of parameters
Included files: 
    - speed_setup_par.jl, 
    - speed_functions_par.jl (which calls Module_ufl_pbf.jl)
    - write_to_tex_function.jl
Output files: 
    - Data/JLD/AES_MILP_times_$(CPU).jld2
    - Data/JLD/brute_times_$(CPU).jld2
    - Plots/Others/speed_AES_MILP_BF_$(CPU).pdf
    - Tables/Others/speed_AES_MILP_BF_$(CPU).tex
"""
include("speed_AES_MILP_BF.jl")

"""
Speed test for AES vs. MILP with varying parameters
    - single-stage, single-product problem
    - serve all or none market
    - fix three sets of locations L and consumers N = L
Included files:
    - speed_setup_par.jl,
    - speed_functions_par.jl (which calls Module_ufl_pbf.jl)
    - write_to_tex_function.jl
Output files:
    - Data/JLD/MILP_CDCP_Var{SigmaV, SigmaU, SigmaF, Beta, SigmaA}_L_$(L_list[1])_$(L_list[end])_X_$(CPU).jld2
"""
include("speed_UFLP_Var.jl")

"""
Tabulate critical parameter values where AES becomes faster than MILP for each parameter type and L (N=L)
Included files:
    - write_to_tex_function.jl
Input files:
    - Data/JLD/MILP_CDCP_Var{SigmaV, SigmaU, SigmaF, Beta, SigmaA}_L_$(L_list[1])_$(L_list[end])_X.jld2
Output files:
    - Tables/Others/speed_critical_sigmas.tex
"""
include("speed_critical_Var.jl")


"""
Plot time comparisons for AES vs. MILP for each parameter type
    - choose one set of locations L and consumers N = L
Included files:
    - plot_stack_function.jl
Input files:
    - Data/JLD/MILP_CDCP_Var{SigmaV, SigmaU, SigmaF, Beta, SigmaA}_L_$(L_list[1])_$(L_list[end])_X.jld2
Output files:
    - Plots/Others/speed_UFLP_Var{SigmaV, SigmaU, SigmaF, Beta, SigmaA}_L$(L)_N$(N).pdf
"""
include("speed_UFLP_Var_plots.jl")

"""
Computational performance of MILP with multiple stages of production 
    - multi-stage, single-product problem
    - endogeneous market entry decisions
    - fix one set of parameters
Included files:
    - speed_setup_par.jl
    - Module_ufl_pbf.jl
    - write_to_tex_function.jl
Output files:
    - Data/JLD/speed_MILP_K2.jld2
    - Data/JLD/speed_MILP_K3.jld2
    - Tables/Others/speed_MILP_K2.tex
    - Tables/Others/speed_MILP_K3.tex
    - Plots/Others/speed_MILP_K2.pdf
    - Plots/Others/speed_MILP_K2_Nvars.pdf
    - Plots/Others/speed_MILP_K3_Nvars.pdf
"""
include("speed_MILP_K2K3.jl")

################################################################################
## Flat world Simulation
################################################################################
"""
Flatworld simulation for the UFLP with two levels of β
    - single-product problem
    - K = 1 or 2 or 3 stages of production
    - endogeneous market entry decisions (except for K = 1)
    - fix parameters other than β
    - fix one set of locations L and consumers N
Included files:
    - speed_setup_par.jl, 
    - speed_functions_par.jl
    - Module_ufl_pbf.jl
    - flatworld_plot.jl
Output files:
    - flatworld_K$(K)_Beta$(par.β).pdf
"""
include("flatworld_UFLP.jl") 

# Some other files in the folder 
# - get_gdp.jl 
# - check_SCD.jl 
# - added_bar_chart.jl: where is this used?