## Optimization Problem : MILP                                  
include("Module_ufl_pbf.jl")
using .Module_ufl_pbf: ufl_K1_fs0, ufl_K1, ufl_pbf_K2, ufl_pbf_K3

"""
## Optimization Problem : AES  
Computes objective function for the AES algorithm.
- J: Binary vector indicating which facilities are open (true) or closed (false)
- fc: Fixed costs K-dimensional array with L[k]-sized vector for each facility location
- mc: Marginal cost matrix (Mf x N x L...)
- A: Demand shifters (Mf x N)
- par: Model parameters (ModelParams struct)
- Returns the total profit for the given facility configuration J
"""
function obj(J, fc, mc, A, par::ModelParams)
    open = findall(J .== true) 
    if isempty(open)
        return 0.0
    end
    mc_chosen = minimum(view(mc, :, :, open), dims=3) # Mf by N by L
    vprofit = vprofit_fn(mc_chosen, A, par)
    return sum(vprofit) - sum(fc[1] .* J) 
end

"""
## Marginal profit of adding/dropping facility j
- j : Index of the facility to consider adding/dropping.
- J: Binary vector indicating which facilities are open (true) or closed (false).
- fc: Fixed costs K-dimensional array with L[k]-sized vector for each facility location
- mc: Marginal cost matrix (Mf x N x L...)
- A: Demand shifters (Mf x N)
- par: Model parameters (ModelParams struct).
- Returns the difference in the objective function from adding or dropping facility j.
"""
function D_j_obj(j, J, fc, mc, A, par::ModelParams)
    bool_j = J[j]
    J = setindex!(J, true, j)
    marg = obj(J, fc, mc, A, par)
    J = setindex!(J, false, j)
    marg -= obj(J, fc, mc, A, par)
    J = setindex!(J, bool_j, j)
    return marg
end

"""
## Identify determined and undetermined locations after the initial squeezing step in AES.
The fewer undetermined, the less branching needed
- L: Number of facility locations
- fc: Fixed costs K-dimensional array with L[k]-sized vector for each facility location
- mc: Marginal cost matrix (Mf x N x L...)
- A: Demand shifters (Mf x N)
- par: Model parameters (ModelParams struct)
- Returns the number of undetermined locations, list of always open locations, and list of always closed locations after squeezing.
"""
function undetermined_loc(L, fc, mc, A, par::ModelParams)
    open_list = [] 
    close_list = []
    open_list_o = 1:L 
    close_list_o = 1:L
    while open_list != open_list_o || close_list != close_list_o
        open_list_o = copy(open_list)
        close_list_o = copy(close_list)
        # start from J all ones
        # always open if in this most competitive case, marginal profit of opening >0
        marg = zeros(Float64, L)
        J = ones(Bool, L)
        J[open_list_o] .= true # determined always open in last iteration
        J[close_list_o] .= false # determined always close in last iteration
        # update marginals only for undetermined locations
        for j in setdiff(1:L, union(open_list_o, close_list_o))
            marg[j] = D_j_obj(j, J, fc, mc, A, par)
        end
        open_list = union(findall(marg .> 0.0), open_list_o)

        # start from J all zeros
        # always close if in this least competitive case, marginal profit of opening <0
        marg = zeros(Float64, L)
        J = zeros(Bool, L)
        J[open_list_o] .= true
        J[close_list_o] .= false
        for j in setdiff(1:L, union(open_list_o, close_list_o))
            marg[j] = D_j_obj(j, J, fc, mc, A, par)
        end
        close_list = union(findall(marg .< 0.0), close_list_o)
    end 
    N_undetermined = L - (length(open_list) + length(close_list))
    return N_undetermined, open_list, close_list
end

"""
## MILP vs. CDCP comparison in time and solutions for a single simulation with given seed index, then use speed_compare2 to run multiple simulations and aggregate results.
# Arguments
- s: seed index
- L: number of facility locations
- N: number of customers
- par: ModelParams struct
# Returns
- DataFrame with results for this seed
"""  
function speed_compare_1sim(s, L, N, par::ModelParams; silent = true, Mf = 1, tol = 1e-6)
    if silent == false
        println("Simulation run ", s)
    end
    seed = s + par.seed_start
    Random.seed!(seed)  
    fc, mc, _, A, _, _, _, _ = setup_par_K(L, N, par, Mf = Mf)
    vprofit = vprofit_fn(mc, A, par)
    results = ones(L, 2)
    Π_list = zeros(2)
    
    # MILP K = 1
    Π_list[1], config, active_y, t_milp, N_nodes, results[:,1] = ufl_K1_fs0(Mf, N, L, fc, vprofit) #  ufl_K1_fs0 default silent = true, no logging

    # AES
    # old syntax:
    # sub = falses(L)
    # sup = trues(L)
    # aux = falses(L)
    # t_aes = @elapsed results[:,2] = solve!((sub, sup, aux); scdca = true, obj = x -> obj(x, fc, mc, A, par), D_j_obj = (x, y) -> D_j_obj(x, y, fc, mc, A, par)) 
    # Π_list[2] = obj(results[:,2], fc, mc, A, par)

    # new syntax: (make sure to install the latest StaticArrays to avoid L limit)
    t_aes = @elapsed AES_solve = solve(Squeezing, J -> obj(J, fc, mc, A, par), L, true) # objective obeys SCD-C from above, submodular
    results[:,2] = AES_solve.x .== included
    Π_list[2] = AES_solve.value

    # total number of AES plants opened
    N_plants = sum(results[:,2],dims=1)

    # check if squeezing determined all locations
    N_undetermined, open_list, close_list = undetermined_loc(L, fc, mc, A, par)
    N_open_list = length(open_list)
    N_close_list = length(close_list)
    check_open = sum(isapprox.(results[open_list,2], 0.0)) # approx to handle float64, check_open = 0 is good, meaning all determined open locations are indeed open
    check_close = sum(isapprox.(results[close_list,2], 1.0)) # check_close = 0 is good, meaning all determined closed locations are indeed closed
    if silent == false
        println("Simulation $s: MILP time $(t_milp) seconds, CDCP time $(t_aes) seconds,  MILP nodes explored $(N_nodes), Undetermined locations after squeezing: $(N_undetermined), Always open locations: $(N_open_list), Always closed locations: $(N_close_list)")
        if N_nodes > 1
            @warn "MILP nodes explored: $(N_nodes)"
        end
        if N_undetermined > 0
            @warn "Undetermined locations after squeezing: $(N_undetermined)"
        end
        if check_open > 0 || check_close > 0
            @warn "Wrong opening: $check_open, Wrong closing: $check_close"
        end
    end

    # check if MILP and AES get the same profit, if not, print a warning
    diff = Π_list[1] - Π_list[2]
    opt_gap = abs(diff) < tol ? 0.0 : diff
    loc_gap = norm(results[:,1] .- results[:,2], 1)
    if loc_gap > tol || opt_gap != 0.0
        @warn "Simulation $s: MILP and AES have different solutions. Location gap: $loc_gap, Profit gap (MILP - AES): $opt_gap"
    end

    return DataFrame(seed = seed, t_milp = t_milp, t_aes = t_aes, opt_gap = opt_gap,loc_gap = loc_gap, N_nodes = N_nodes, N_undetermined = N_undetermined, N_open_list = N_open_list, N_close_list = N_close_list, N_plants = N_plants)
end

"""
## MILP vs. CDCP comparison in time and solutions for multiple simulations sequentially or in parallel and aggregate results.
# NOTE: this function requires ThreadsX.jl or Distributed.jl for parallel processing, make sure to set up the environment accordingly
# Arguments
- L: Number of facility locations
- N: Number of customers
- par: ModelParams struct
- mode: :serial, :distributed, or :threads for how to run multiple simulations
- silent: if false, print detailed results for each simulation
- aggregate: if true, return average results across simulations; if false, return a concatenated DataFrame with each row a simulation
- Mf: number of products per firm (default 1)
- tol: tolerance for checking solution differences (default 1e-10)
# Returns
- aggregated results over all simulations
- or the DataFrame concatenating all simulation results if aggregate=false
"""

function speed_compare(L, N, par::ModelParams; mode=:serial, silent = true, aggregate = true, Mf = 1, tol = 1e-6)
    Nsim = par.num_sims

    if mode == :serial
        # (1) Run simulations sequentially
        df_list = map(s -> speed_compare_1sim(s, L, N, par; silent = silent, Mf = Mf, tol = tol), 1:Nsim)
    elseif mode == :distributed
        # (2) Run simulations in parallel using Distributed
        df_list = pmap(s -> speed_compare_1sim(s, L, N, par; silent = silent, Mf = Mf, tol = tol), 1:Nsim)
    elseif mode == :threads
        # (3) Run simulations in parallel using ThreadsX
        df_list = ThreadsX.map(s -> speed_compare_1sim(s, L, N, par; silent = silent, Mf = Mf, tol = tol), 1:Nsim)
    else
        error("mode must be :serial, :distributed, or :threads")
    end

    df = vcat(df_list...)
    # if aggregate is false, return per-simulation results as a DataFrame
    if aggregate == false
        return df
    end

    avg_t_milp = mean(df.t_milp)
    avg_t_aes = mean(df.t_aes)
    avg_opt_gap = mean(df.opt_gap) # optimal profit difference
    avg_loc_gap = mean(df.loc_gap) # optimal location vector difference
    std_t_milp = std(df.t_milp)
    std_t_aes = std(df.t_aes)
    avg_N_nodes = mean(df.N_nodes)
    avg_N_undet = mean(df.N_undetermined)
    avg_N_open = mean(df.N_open_list)
    avg_N_close = mean(df.N_close_list)
    avg_N_plants = mean(df.N_plants)

    return avg_t_milp, avg_t_aes, avg_opt_gap, avg_loc_gap, std_t_milp, std_t_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants
end


"""
## Prepare for Brute Force
Enumerate all possible facility opening configurations for L locations.
- L: Number of facility locations
- Returns a vector of all possible binary vectors of length L (except the all-closed case).
"""                                  
function get_loc_set(L)
    arg = []
    for l in 1:L
        push!(arg, [false, true])
    end
    loc_set = vec(collect(Iterators.product(arg...)))
    loc_set = loc_set[2:end] # remove the case where all facilities are closed
    return loc_set
end

"""
## Optimization Problem : Brute Force
L<=20 is easy, L>25 is hard and also require large memory (~200Gb), L>30 will probably choke computer.
- L: Number of facility locations
- N: Number of customers
- par: ModelParams struct
- Returns average time and memory usage for brute force enumeration of all facility configurations, including the time and memory to generate the location set.
"""
function speed_brute(L, N, par::ModelParams)
    Nsim = par.num_sims
    loc_set, t_loc_set, m_loc_set = @timed get_loc_set(L);
    println("To enumerate all possible location configurations required ", round(t_loc_set,digits =4)," seconds and ", round(m_loc_set/10^9,digits=4), " GB of memory")
  
    t_bf = zeros(Nsim);
    m_bf = zeros(Nsim);
    for s in 1:Nsim
        Random.seed!(s + par.seed_start)
        fc, mc, _, A, _, _, _, _ = setup_par_K(L, N, par)
        results = ones(L)
        timed_tuple = @timed begin
            πs = []
            for J in loc_set 
                π = obj(J, fc, mc, A, par)
                push!(πs, π)
            end
            results = [loc_set[argmax(πs)]...]
        end

        # CDCP.jl syntax for Brute Force:
        # obj_wrapped = Objective(J -> obj(J, fc, mc, A, par), Vector{Bool}(falses(L)))
        # naive_solve = solve(Naive, obj_wrapped, L)
        # naive_solve.x
        # naive_solve.value

        t_bf[s] = timed_tuple.time
        m_bf[s] = timed_tuple.bytes/10^9 #in GB
    end
    avg_t_bf = mean(t_bf)+t_loc_set
    sd_t_bf = std(t_bf)
    avg_m_bf = mean(m_bf)+m_loc_set/10^9
    return avg_t_bf, sd_t_bf, avg_m_bf
end

"""
getCPU() returns the CPU model as a string excluding the brand ("Apple") and any leading numbers. For example, if the CPU model is "Apple M1 Pro", getCPU() will return "M1Pro".
- useful for labelling results,
- example: filename = "speed_check_$(getCPU())"
"""
getCPU() = join(split(Sys.cpu_info()[1].model)[2:end])

"""
Solve MILP and check computational performance for K = 1 or K=2 or K=3 for a single simulation with given seed index
"""
function solve_MILP_1sim(s, L, N, par::ModelParams; silent = true, Mf = 1)
    if silent == false
        println("Simulation run ", s)
    end
    seed = s + par.seed_start
    Random.seed!(seed)  
    fc, mc, fs, A, Xc, Yc, Xf, Yf = setup_par_K(L, N, par, Mf = Mf)
    vprofit = vprofit_fn(mc, A, par)

    # Solve MUFLP-PBF for K = 2 or K = 3 only:
    K = length(L)
    if K == 1 && par.μₙ == 0.0
        Π, config, active_y, t_milp, N_nodes, _ = ufl_K1_fs0(Mf, N, L, fc, vprofit)
    elseif K == 1 && par.μₙ != 0.0
        Π, config, active_y, t_milp, N_nodes, _ = ufl_K1(Mf, N, L, fc, fs, vprofit)
    elseif K == 2
        Π, config, active_y, t_milp, N_nodes = ufl_pbf_K2(Mf, N, L, K, fc, fs, vprofit)
    elseif K == 3
        Π, config, active_y, t_milp, N_nodes = ufl_pbf_K3(Mf, N, L, K, fc, fs, vprofit)
    end

    return (Π=Π, config=config, active_y=active_y, Xc=Xc, Yc=Yc, Xf=Xf, Yf=Yf, fc=fc, A=A, t_milp=t_milp, N_nodes=N_nodes)
end

"""
MILP speed test for multiple simulations, with options to run sequentially or in parallel using Distributed.jl or ThreadsX.jl
"""
function speed_MILP(L, N, par::ModelParams; mode=:serial, silent = true, Mf = 1)
    Nsim = par.num_sims

    if mode == :serial
        # (1) Run simulations sequentially
        results = map(s -> solve_MILP_1sim(s, L, N, par; silent = silent, Mf = Mf), 1:Nsim)
    elseif mode == :distributed
        # (2) Run simulations in parallel using Distributed
        results = pmap(s -> solve_MILP_1sim(s, L, N, par; silent = silent, Mf = Mf), 1:Nsim)
    elseif mode == :threads
        # (3) Run simulations in parallel using ThreadsX
        results = ThreadsX.map(s -> solve_MILP_1sim(s, L, N, par; silent = silent, Mf = Mf), 1:Nsim)
    else
        error("mode must be :serial, :distributed, or :threads")
    end

    time_list = [r.t_milp for r in results]
    avg_t_milp = mean(time_list)
    std_t_milp = std(time_list)
    node_list = [r.N_nodes for r in results]
    avg_N_nodes = mean(node_list)

    return avg_t_milp, std_t_milp, avg_N_nodes, time_list, node_list
end