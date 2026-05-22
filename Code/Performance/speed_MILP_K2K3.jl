# Time MILP by varying dimension of the problem
using JuMP
import Gurobi
const ENV = Gurobi.Env()
using DataFrames
import DataFrames: groupby
using Plots
using JLD2
using Distributed, ThreadsX

# Import auxiliary functions
include("speed_setup_par.jl")
include("speed_functions_par.jl")
include("write_to_tex_function.jl")

CPU = getCPU()

# Parameter setting
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 9.5, num_sims = 100) 

## K = 2
# Define the dimension of the problem
L_list = [10, 20, 30, 40, 60, 80, 100]
N_list = [20, 40]
M_list = [1, 10, 20]

# Warmup once per L to avoid overhead of first run
for L in L_list
    N_warm = 10 # a small N for warmup
    speed_MILP([L L], N_warm, par, mode = :threads) # default Mf = 1
end

df_K2 = DataFrame(L = Int[], N = Int[], M = Int[], avg_time = Float64[], std_time = Float64[], avg_nodes = Float64[], time_list = Vector{Float64}[], node_list = Vector{Int}[])
for N in N_list
    for M in M_list
        for L in L_list
            t, s, nodes, time_list, node_list = speed_MILP([L L], N, par, Mf = M, mode = :threads)
            println("L = $L, N = $N, M = $M, avg_time = $t, std_time = $s, avg_nodes = $nodes")
            push!(df_K2, [L, N, M, t, s, nodes, time_list, node_list])
            GC.gc()
        end
    end
end
# save to jld2 file
file_path = "Data/JLD/speed_MILP_K2_$(CPU).jld2"
save(file_path, "df_K2", df_K2)

## K = 3
# Define the dimension of the problem
L_list = [10, 20, 30]
N_list = [20, 40]
M_list = [1, 10, 20]

# Warmup once per L to avoid overhead of first run
for L in L_list
    N_warm = 10 # a small N for warmup
    speed_MILP([L L L], N_warm, par, mode = :threads) # default Mf = 1
end

df_K3 = DataFrame(L = Int[], N = Int[], M = Int[], avg_time = Float64[], std_time = Float64[], avg_N_nodes = Float64[], time_list = Vector{Float64}[], node_list = Vector{Int}[])
for N in N_list
    for M in M_list
        for L in L_list
            t, s, nodes, time_list, node_list = speed_MILP([L L L], N, par, Mf = M, mode = :threads)
            println("L = $L, N = $N, M = $M, avg_time = $t, std_time = $s, avg_N_nodes = $nodes")
            push!(df_K3, [L, N, M, t, s, nodes, time_list, node_list], promote=true)
        end
    end
end
# save to jld2 file
file_path = "Data/JLD/speed_MILP_K3_$(CPU).jld2"
save(file_path, "df_K3", df_K3)

## Output the results
@load "Data/JLD/speed_MILP_K2_$(CPU).jld2" df_K2
@load "Data/JLD/speed_MILP_K3_$(CPU).jld2" df_K3

# Output the results to a LaTeX table
groups = groupby(df_K2[:, [:N, :M, :L, :avg_time]], [:N, :M])
df_K2_tex = hcat([select(g, Not(:N, :M)) for g in groups]..., makeunique=true)
# drop L_1, L_2 columns 
select!(df_K2_tex, Not(r"L_\d+"))
write_to_tex(df_K2_tex, "Tables/Others/speed_MILP_K2.tex", rounding = [0, 2, 2, 2, 2, 2, 2])
groups = groupby(df_K3[:, [:N, :M, :L, :avg_time]], [:N, :M]) # :avg_nodes
df_K3_tex = hcat([select(g, Not(:N, :M)) for g in groups]..., makeunique=true)
select!(df_K3_tex, Not(r"L_\d+"))
write_to_tex(df_K3_tex, "Tables/Others/speed_MILP_K3.tex", rounding = [0, 2, 2, 2, 2, 2, 2])

# Output the results to a plot
maxT = maximum(df_K2.avg_time)
maxL = maximum(df_K2.L)
pmarg = [5 30]
p = plot(ylabel = "Time in seconds", xlabel = "Number of locations per stage", xlims=(0, maxL+pmarg[2]), ylims= (0,maxT+pmarg[1]), xticks = unique(df_K2.L));
mycolors = [:blue, :orange, :brown, :red, :purple, :green];
# plot K = 2
N_list = [20, 40]
M_list = [1, 10, 20]
i = 0
for n in N_list
    for m in M_list
        i = i + 1
        plot!(df_K2[df_K2.N .== n .&& df_K2.M .== m, :L], df_K2[df_K2.N .== n .&& df_K2.M .== m, :avg_time],
        legend = false,
        ms = 3,
        mc = mycolors[i], lc = mycolors[i],
        shape = :circle)
        annotate!(maxL+pmarg[1], maximum(df_K2.avg_time[df_K2.N .== n .&& df_K2.M .== m]), text("N=$(n), M=$(m)", :left, mycolors[i], 11))
    end
end
display(p)
savefig(p,"Plots/Others/speed_MILP_K2.pdf")

# Check time as a function of Nvar or Nconstraints for K = 2
K = 2
# number of variables
df_K2.Nvar = df_K2.N .* df_K2.M .* (df_K2.L).^K + df_K2.L * K + df_K2.N .* df_K2.M
# number of paths (in thousands)
df_K2.Npath_ths = (df_K2.N .* df_K2.M .* (df_K2.L).^K) / 10^3
# number of constraints
df_K2.Nconst = df_K2.N .* df_K2.M + df_K2.N .* df_K2.M .* df_K2.L * K
# plot
ticks_p = [2, 10, 50, 250, 1000, 3000, 8000]
ticks_t = [0.01, 0.1, 1, 10, 50]
p = scatter(df_K2.Npath_ths, df_K2.avg_time, ylabel = "Time in seconds (log scale)", xlabel = "Number of paths (in ths, log scale)", yscale=:ln, xscale=:ln, legend = false, xticks=(ticks_p,ticks_p),yticks=(ticks_t,ticks_t)) # looks linear 
df_K2.time_linpred = df_K2.Npath_ths .* (df_K2.avg_time[1] / df_K2.Npath_ths[1])
df_K2_sorted = sort(df_K2, :Npath_ths)
plot!(df_K2_sorted.Npath_ths, df_K2_sorted.time_linpred, linestyle = :dash, color = :red)
annotate!(550,1.0,text("linear prediction \\n (time prop. to paths)", :left, :red, 11))
display(p)
savefig(p,"Plots/Others/speed_MILP_K2_Nvars.pdf")

# using GLM
# linear_model = lm(@formula(log(avg_time) ~ log(Npath_ths)), df_K2)

# Check time as a function of Nvar or Nconstraints for K = 3
K = 3
# number of variables
df_K3.Nvar = df_K3.N .* df_K3.M .* (df_K3.L).^K + df_K3.L * K + df_K3.N .* df_K3.M
# number of paths (in thousands)
df_K3.Npath_ths = (df_K3.N .* df_K3.M .* (df_K3.L).^K) / 10^3
# number of constraints
df_K3.Nconst = df_K3.N .* df_K3.M + df_K3.N .* df_K3.M .* df_K3.L * K
sort!(df_K3, :Npath_ths)
# Plot
p = scatter(df_K3.Npath_ths, df_K3.avg_time, ylabel = "Time in seconds (log scale)", xlabel = "Number of paths (in ths, log scale)", yscale=:ln, xscale=:ln, legend = false) # looks linear
df_K3.time_linpred = df_K3.Npath_ths .* (df_K3.avg_time[1] / df_K3.Npath_ths[1])
plot!(df_K3.Npath_ths, df_K3.time_linpred, linestyle = :dash, color = :red)
annotate!(minimum(df_K3.Npath_ths)+1, minimum(df_K3.avg_time),text("linear prediction", :left, :red, 11))
display(p)
savefig(p,"Plots/Others/speed_MILP_K3_Nvars.pdf")

## Tabulate node count
@load "Data/JLD/speed_MILP_K2_$(CPU).jld2" df_K2
@load "Data/JLD/speed_MILP_K3_$(CPU).jld2" df_K3
using StatsBase
K2_nodes = vcat(df_K2.node_list...)
K2_freq = countmap(K2_nodes)
println("K = 2 node count frequency: ", K2_freq)
K3_nodes = vcat(df_K3.node_list...) # the last setup L = 30, N = 40, M = 20 has NaN nodes
K3_freq = countmap(K3_nodes)
println("K = 3 node count frequency: ", K3_freq)

# log the node count frequency to a file
open("Gurobi_logs/node_freq.txt", "a") do io
    println(io, "CPU: ", CPU)
    println(io, "K = 2 node count frequency: ", sort(collect(K2_freq)))
    println(io, "K = 3 node count frequency: ", sort(collect(K3_freq)))
    println(io, "Note: K=3 last setup (L=30, N=40, M=20) has NaN nodes when there is no optimal solution found within the time limit.")
end
