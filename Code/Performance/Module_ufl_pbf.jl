module Module_ufl_pbf
export ufl_K1_fs0, ufl_K1, ufl_pbf_K2, ufl_pbf_K3

using JuMP
import Gurobi
const ENV = Gurobi.Env()

function ufl_K1_fs0(Mf, N, L, fc, vprofit; silent = true, logfile = "Gurobi_logs/ufl_pbf_log.txt", num_threads::Int=1)
    y_ = zeros(Int, L)
    x_ = zeros(Mf, N, L)
    
    ufl = Model(() -> Gurobi.Optimizer(ENV))
    set_optimizer_attribute(ufl, "LogToConsole", 0)
    if silent == true
        set_silent(ufl)
    else
        set_optimizer_attribute(ufl, "LogFile", logfile)
    end
    set_optimizer_attribute(ufl, "WorkLimit", 200.0)
    set_optimizer_attribute(ufl, "Presolve", 0)
    set_optimizer_attribute(ufl, "Threads", num_threads)
    # set_optimizer_attribute(ufl, "MIPGap", 1e-4) # when there is u shocks, MIPGap should be relaxed, default is 1e-4
    set_optimizer_attribute(ufl, "MIPGap", 0)
    set_optimizer_attribute(ufl, "OptimalityTol", 1e-9)
    set_optimizer_attribute(ufl, "Method", 1) 
    # Variables
    @variable(ufl, y[1:L], Bin); # define binary y (active plants)
    @variable(ufl, x[1:Mf, 1:N, 1:L] >= 0); # define binary x (sourcing assignments)

    # Flow-conservation constraints
    # Each client is served exactly once, i.e. each consumer must be served by exactly one facility
    # When fixed cost of market entry is zero, markets either all enter or all exit, so we don't need to define z
    @constraint(ufl, client_service[m in 1:Mf, n in 1:N],
        sum(x[m, n, l] for l in 1:L) <= 1 
    );

    # Activity constraints
    # A facility must be open to serve a client
    @constraint(ufl, open_facility[m in 1:Mf, n in 1:N, l in 1:L],
        x[m, n, l] <= y[l]
    )

    # Objective: profit maximization
    @objective(ufl, Max, sum(vprofit .* x) - sum(fc[1] .* y));

    t_milp = @elapsed optimize!(ufl)

    # Count nodes explored
    node_count = MOI.get(ufl, MOI.NodeCount())

    # Check if the solution is optimal
    if termination_status(ufl) != MOI.OPTIMAL
        @warn "No optimal solution. Status: $(termination_status(ufl))"
        node_count = NaN # set node count to NaN if no optimal solution is found
    end

    y_ = round.(Int, value.(y));
    x_ = round.(Int, value.(x));
    # total profit
    Π = objective_value(ufl) + 0.0; # turn -0.0 to 0.0

    # update active assignments
    active_MN = findall(x -> x > 0, sum(x_,dims=3))
    active_MN = hcat(getindex.(active_MN, 1), getindex.(active_MN, 2))
    # update active assignments based on market entry decision
    config = []
    for m in 1:Mf
        active_N = active_MN[findall(x -> x == m, active_MN[:,1]), 2]
        tmp = zeros(Int, N)
        for n in active_N 
            tmp[n] = findall(x -> x > 0, x_[m,n,:,:])[1][1]
        end
        push!(config, tmp)
    end

    # update active plants based on all entry or all exit
    active_y = [unique(filter(x -> x != 0, vcat(config...)))]

    return Π, config, active_y, t_milp, node_count, y_
end

function ufl_K1(Mf, N, L, fc, fs, vprofit; silent = true, logfile = "Gurobi_logs/ufl_pbf_log.txt", num_threads::Int=1)
    z_ = zeros(Int, Mf, N)
    y_ = zeros(Int, L)
    x_ = zeros(Mf, N, L)
    
    ufl = Model(() -> Gurobi.Optimizer(ENV))
    set_optimizer_attribute(ufl, "LogToConsole", 0)
    if silent == true
        set_silent(ufl)
    else
        set_optimizer_attribute(ufl, "LogFile", logfile)
    end
    set_optimizer_attribute(ufl, "WorkLimit", 200.0)
    set_optimizer_attribute(ufl, "Presolve", 0)
    set_optimizer_attribute(ufl, "Threads", num_threads)
    # set_optimizer_attribute(ufl, "MIPGap", 1e-4) # when there is u shocks, MIPGap should be relaxed, default is 1e-4
    set_optimizer_attribute(ufl, "MIPGap", 0)
    set_optimizer_attribute(ufl, "OptimalityTol", 1e-9)
    set_optimizer_attribute(ufl, "Method", 1) 
    # Variables
    @variable(ufl, z[1:Mf, 1:N], Bin); # define binary z (market entry decision)
    @variable(ufl, y[1:L], Bin); # define binary y (active plants)
    @variable(ufl, x[1:Mf, 1:N, 1:L] >= 0); # define binary x (sourcing assignments)

    # Flow-conservation constraints
    # Each client is served exactly once, i.e. each consumer must be served by exactly one facility
    @constraint(ufl, client_service[m in 1:Mf, n in 1:N],
        sum(x[m, n, l] for l in 1:L) <= z[m, n] 
    );

    # Activity constraints
    # A facility must be open to serve a client
    @constraint(ufl, open_facility[m in 1:Mf, n in 1:N, l in 1:L],
        x[m, n, l] <= y[l]
    )

    # Objective: profit maximization
    @objective(ufl, Max, sum(vprofit .* x) - sum(fc[1] .* y) - sum(fs' .* z));

    t_milp = @elapsed optimize!(ufl)

    # Count nodes explored
    node_count = MOI.get(ufl, MOI.NodeCount())

    # Check if the solution is optimal
    if termination_status(ufl) != MOI.OPTIMAL
        @warn "No optimal solution. Status: $(termination_status(ufl))"
        node_count = NaN # set node count to NaN if no optimal solution is found
    end

    z_ = round.(Int, value.(z));
    y_ = round.(Int, value.(y));
    x_ = round.(Int, value.(x));
    # total profit
    Π = objective_value(ufl) + 0.0; # turn -0.0 to 0.0

    # update active assignments
    active_MN = findall(x -> x > 0, z_)
    active_MN = hcat(getindex.(active_MN, 1), getindex.(active_MN, 2))
    # update active assignments based on market entry decision
    config = []
    for m in 1:Mf
        active_N = active_MN[findall(x -> x == m, active_MN[:,1]), 2]
        tmp = zeros(Int, N)
        for n in active_N 
            tmp[n] = findall(x -> x > 0, x_[m,n,:,:])[1][1]
        end
        push!(config, tmp)
    end

    # update active plants based on all entry or all exit
    active_y = [unique(filter(x -> x != 0, vcat(config...)))]

    return Π, config, active_y, t_milp, node_count, y_
end

function ufl_pbf_K2(Mf, N, L, K, fc, fs, vprofit; silent = true, logfile = "Gurobi_logs/ufl_pbf_log.txt", num_threads::Int=1)
    z_ = zeros(Int, Mf, N)
    y1_ = zeros(Int, L[1])
    y2_ = zeros(Int, L[2])
    x_ = zeros(Mf, N, L[1], L[2])

    ufl = Model(() -> Gurobi.Optimizer(ENV))
    set_optimizer_attribute(ufl, "LogToConsole", 0)
    if silent == true
        set_silent(ufl)
    else
        set_optimizer_attribute(ufl, "LogFile", logfile) 
    end
    set_optimizer_attribute(ufl, "WorkLimit", 200.0)
    set_optimizer_attribute(ufl, "Presolve", 0)
    set_optimizer_attribute(ufl, "Threads", num_threads)
    set_optimizer_attribute(ufl, "MIPGap", 0) # when there is u shocks, MIPGap should be relaxed, default is 1e-4
    set_optimizer_attribute(ufl, "OptimalityTol", 1e-9)
    set_optimizer_attribute(ufl, "Method", 1) 
    # Variables
    @variable(ufl, z[1:Mf, 1:N], Bin);
    @variable(ufl, y1[1:L[1]], Bin);
    @variable(ufl, y2[1:L[2]], Bin);
    @variable(ufl, x[1:Mf, 1:N, 1:L[1], 1:L[2]] >= 0);

    # Flow-conservation constraints
    # Each client is served exactly once, i.e. each consumer must be served by exactly one facility
    @constraint(ufl, client_service[m in 1:Mf, n in 1:N],
        sum(x[m, n, l1, l2] for l1 in 1:L[1], l2 in 1:L[2]) <= z[m, n] 
    )
    # Capacity constraints
    # A facility must be open to serve a client
    #assembly: if l2 is closed, no path can pass through it
    @constraint(ufl, open_facility2[m in 1:Mf, n in 1:N, l2 in 1:L[2]],
        sum(x[m, n, l1, l2] for l1 in 1:L[1]) <= y2[l2]
    )
    #cells
    @constraint(ufl, open_facility1[m in 1:Mf, n in 1:N, l1 in 1:L[1]],
        sum(x[m, n, l1, l2] for l2 in 1:L[2]) <= y1[l1]
    )
    
    # Objective
    @objective(ufl, Max, sum(vprofit .* x) - sum(fc[1] .* y1) - sum(fc[2] .* y2) - sum(fs' .* z));

    t_milp = @elapsed optimize!(ufl)

    # Count nodes explored
    node_count = MOI.get(ufl, MOI.NodeCount())

    # Check if the solution is optimal
    if termination_status(ufl) != MOI.OPTIMAL
        @warn "No optimal solution. Status: $(termination_status(ufl))"
        node_count = NaN # set node count to NaN if no optimal solution is found
    end

    z_ = round.(Int, value.(z));
    y1_ = round.(Int, value.(y1));
    y2_ = round.(Int, value.(y2));
    x_ = round.(Int, value.(x));
    # total profit
    Π = objective_value(ufl) + 0.0; # turn -0.0 to 0.0

    active_MN = findall(x -> x > 0, z_)
    active_MN = hcat(getindex.(active_MN, 1), getindex.(active_MN, 2))
    # update active assignments based on market entry decision
    config = []
    for m in 1:Mf
        active_N = active_MN[findall(x -> x == m, active_MN[:,1]), 2]
        tmp = zeros(Int, N, K)
        for n in active_N
            for k in 1:K    
                tmp[n,k] = findall(x -> x > 0, x_[m,n,:,:])[1][k]
            end
        end
        push!(config, tmp)
    end
    # update active plants based on market entry decision
    active_y = []
    for k in 1:K
        push!(active_y, unique(filter(x -> x != 0, vcat(config...)[:,k])))
    end

    return Π, config, active_y, t_milp, node_count
end

#ufl_pbf K=3
function ufl_pbf_K3(Mf, N, L, K, fc, fs, vprofit; silent = true, logfile = "Gurobi_logs/round/ufl_pbf_log.txt", num_threads::Int=1)
    z_ = zeros(Int, Mf, N)
    y1_ = zeros(Int, L[1])
    y2_ = zeros(Int, L[2])
    y3_ = zeros(Int, L[3])
    x_ = zeros(Mf, N, L[1], L[2], L[3])

    ufl = Model(() -> Gurobi.Optimizer(ENV))
    set_optimizer_attribute(ufl, "LogToConsole", 0)
    if silent == true
        set_silent(ufl)
    else
        set_optimizer_attribute(ufl, "LogFile", logfile)
    end
    set_optimizer_attribute(ufl, "WorkLimit", 600.0)
    set_optimizer_attribute(ufl, "Presolve", 0)
    set_optimizer_attribute(ufl, "Threads", num_threads)
    set_optimizer_attribute(ufl, "MIPGap", 0) # when there is u shocks, MIPGap should be relaxed, default is 1e-4
    set_optimizer_attribute(ufl, "OptimalityTol", 1e-9)
    set_optimizer_attribute(ufl, "Method", 1) 
    # Variables
    @variable(ufl, z[1:Mf, 1:N], Bin);
    @variable(ufl, y1[1:L[1]], Bin);
    @variable(ufl, y2[1:L[2]], Bin);
    @variable(ufl, y3[1:L[3]], Bin);
    @variable(ufl, x[1:Mf, 1:N, 1:L[1], 1:L[2], 1:L[3]] >= 0);

    # Flow-conservation constraints
    # Each client is served exactly once, i.e. each consumer must be served by exactly one facility
    @constraint(ufl, client_service[m in 1:Mf, n in 1:N],
        sum(x[m, n, l1, l2, l3] for l1 in 1:L[1], l2 in 1:L[2], l3 in 1:L[3]) <= z[m, n] 
    )
    # Capacity constraints
    # A facility must be open to serve a client
#k=3 (assembly): if l3 is closed, no path can pass through it
    @constraint(ufl, open_facility3[m in 1:Mf, n in 1:N, l3 in 1:L[3]],
        sum(x[m, n, l1, l2, l3] for l1 in 1:L[1], l2 in 1:L[2]) <= y3[l3]
    )
#k =2 (packs)
    @constraint(ufl, open_facility2[m in 1:Mf, n in 1:N, l2 in 1:L[2]],
        sum(x[m, n, l1, l2, l3] for l1 in 1:L[1], l3 in 1:L[3]) <= y2[l2]
    )
#k=1 (cells)
    @constraint(ufl, open_facility1[m in 1:Mf, n in 1:N, l1 in 1:L[1]],
        sum(x[m, n, l1, l2, l3] for l2 in 1:L[2], l3 in 1:L[3]) <= y1[l1]
    )
    
    # Objective
    @objective(ufl, Max, sum(vprofit .* x) - sum(fc[1] .* y1) - sum(fc[2] .* y2) - sum(fc[3] .* y3) - sum(fs' .* z));

    t_milp = @elapsed optimize!(ufl)

    # Count nodes explored
    node_count = MOI.get(ufl, MOI.NodeCount())

    # Check if the solution is optimal
    if termination_status(ufl) != MOI.OPTIMAL
        @warn "No optimal solution. Status: $(termination_status(ufl))"
        node_count = NaN # set node count to NaN if no optimal solution is found
    end

    z_ = round.(Int, value.(z));
    y1_ = round.(Int, value.(y1));
    y2_ = round.(Int, value.(y2));
    y3_ = round.(Int, value.(y3));
    x_ = round.(Int, value.(x));
    # total profit
    Π = objective_value(ufl) + 0.0; # turn -0.0 to 0.0

    active_MN = findall(x -> x > 0, z_)
    active_MN = hcat(getindex.(active_MN, 1), getindex.(active_MN, 2))
    # update active assignments based on market entry decision
    config = []
    for m in 1:Mf
        active_N = active_MN[findall(x -> x == m, active_MN[:,1]), 2]
        tmp = zeros(Int, N, K)
        for n in active_N
            for k in 1:K    
                tmp[n,k] = findall(x -> x > 0, x_[m,n,:,:,:])[1][k]
            end
        end
        push!(config, tmp)
    end
    # update active plants based on market entry decision
    active_y = []
    for k in 1:K
        push!(active_y, unique(filter(x -> x != 0, vcat(config...)[:,k])))
    end

    return Π, config, active_y, t_milp, node_count
end

end # module