# This code is for single-stage uncapacitated facility location problem using MILP to be compared with Arkolakis, Eckert and Shi (2023), and also Brute Force
# for a particular set of parameter settings and varying L and N.

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
include("write_to_tex_function.jl")

CPU = getCPU()

# Set the parameters for the speed test K = 1
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 4.5, μₙ = 0.0, num_sims = 300)

## AES vs MILP performances

# set L_list and N_list:
L_list = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
N_list = [20, 40]

# Warmup once per L (the size of location array to avoid AES Julia JIT compilation, CDCP uses StaticArrays)
for L in L_list
    N = 10 # a small N for warmup
    speed_compare(L, N, ModelParams(par; σᵥ = 0.5), mode = :threads)
end

# Actual speed test
df = DataFrame(L = Int[], N = Int[], avg_time_milp = Float64[], avg_time_aes = Float64[], opt_gap_milp_aes = Float64[], loc_gap_milp_aes = Float64[], std_time_milp = Float64[], std_time_aes = Float64[], avg_N_nodes = Float64[], avg_N_undet = Float64[], avg_N_open = Float64[], avg_N_close = Float64[], avg_N_plants = Float64[])
@time for N in N_list
    for L in L_list
        tmilp, taes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants = speed_compare(L, N, par, mode = :threads)
        push!(df, [L, N, tmilp, taes, opt_gap, loc_gap, sdt_milp, sdt_aes, avg_N_nodes, avg_N_undet, avg_N_open, avg_N_close, avg_N_plants])
        println("L = $L, N = $N, time_milp = $(round(tmilp; digits = 3)), time_aes = $(round(taes; digits = 3)), opt gap: milp-aes = $opt_gap, loc gap: milp-aes = $loc_gap, MILP nodes = $(avg_N_nodes), CDCP undet = $(avg_N_undet), always open = $(avg_N_open), always closed = $(avg_N_close), plants opened = $(avg_N_plants)")
        GC.gc()
    end
end 
# check output for opt gaps or branch and bounding, we want 0,0,1.0
maximum(df.loc_gap_milp_aes)
maximum(abs.(df.opt_gap_milp_aes))
maximum(df.avg_N_nodes)
# save to jld2 file
file_path = "Data/JLD/AES_MILP_times_sigmaV$(par.σᵥ)_$(CPU).jld2"
save(file_path, "df", df)

## Brute Force performance test
L_list = [10, 12, 14, 16, 18, 20]
N_list = [20]

dfbf = DataFrame(L = Int[], N = Int[], avg_time_bf = Float64[], sd_time_bf = Float64[], avg_mem_bf = Float64[])
@time for N in N_list
    for L in L_list
        tbf, sdtbf, mbf = speed_brute(L, N, par)
        push!(dfbf, [L, N, tbf, sdtbf, mbf])
        println("L = $L, N = $N, time_bf = $(round(tbf; digits = 4)), sdt_bf = $(round(sdtbf; digits = 4)), m_bf = $(round(mbf; digits = 3))")
        GC.gc()
    end
end
# save to jld2 file
file_path = "Data/JLD/brute_times_sigmaV$(par.σᵥ)_$(CPU).jld2"
save(file_path, "dfbf", dfbf)

## Output the results to a plot
dfbf = JLD2.load("Data/JLD/brute_times_sigmaV$(par.σᵥ)_$(CPU).jld2", "dfbf")
df = JLD2.load("Data/JLD/AES_MILP_times_sigmaV$(par.σᵥ)_$(CPU).jld2", "df")
ticks_p = unique([10, 20, unique(df.L)...])
ticks_t = [0.001, 0.01, 0.1, 1, 10]
maxL = maximum(ticks_p)
minL = minimum(ticks_p)
# maxL = maximum(df.L)
pmarg = [5 100]
include_bf = true
if include_bf == true 
    # maxy = 1.0
    maxy = maximum(ticks_t)
else 
    maxy = maximum(vcat(df.avg_time_milp, df.avg_time_aes)) + 0.1
end
p = scatter(ylabel = "Time in seconds (log scale)", xlabel = "Number of locations (log scale)", yscale=:ln, xscale=:ln, legend = false, xticks=(ticks_p,ticks_p),yticks=(ticks_t,ticks_t), xlims=(minimum(ticks_p)-1, maxL+pmarg[2]))
myshapes = [:circle, :square, :diamond, :utriangle];

# plot AES & MILP
i = 0
for n in unique(df.N)
    i = i+1
    plot!(df[df.N .== n,:L], df[df.N .== n,:avg_time_milp],
    legend = false, 
    ms=3, 
    mc=:blue,lc = :blue,
    shape=myshapes[i])
    plot!(df[df.N .== n,:L], df[df.N .== n,:avg_time_aes],
    legend = false, 
    ms=3, 
    mc=:orange,lc = :orange,
    shape=myshapes[i])
    #
    annotate!(maxL+pmarg[1], maximum(df.avg_time_milp[df.N .== n]),text("N = $n",:left,:blue, 11))
    annotate!(maxL+pmarg[1], maximum(df.avg_time_aes[df.N .== n]),text("N = $n",:left,:orange, 11))
end

# plot Brute force
if include_bf == true
    i = 0
    for n in unique(dfbf.N)
        i = i + 1
        plot!(dfbf[dfbf.N.==n, :L], dfbf[dfbf.N.==n, :avg_time_bf],
            legend=false,
            ms=3,
            mc=:red, lc=:red,
            shape=myshapes[i])

        annotate!(22, 0.95 * maxy, text("N = $n", :left, :red, 11))
    end
end
    # legend
if include_bf == true    
    annotate!(maxL + pmarg[2], maxy, text("BF", :right, :red, 11))
    annotate!(maxL + pmarg[2], 0.5 * maxy, text("AES", :right, :orange, 11))
    annotate!(maxL + pmarg[2], 0.25 * maxy, text("MILP", :right, :blue, 11))
else
    L_label = L_list[end]
    y_milp = df[(df.N .== maximum(N_list)) .&& (df.L .== L_label), :avg_time_milp]
    y_aes = df[(df.N .== maximum(N_list)) .&& (df.L .== L_label), :avg_time_aes]
    annotate!(L_label-10,y_aes , text("AES", :right, :orange, 11))
    annotate!(L_label-10,y_milp , text("MILP", :right, :blue, 11))
end
# add parameters as title on top
plot!(p,
    title = "σᵥ=$(par.σᵥ), σᵤ=$(par.σᵤ), σᵩ=$(par.σᵩ), μᵩ=$(par.μᵩ), μₐ=$(par.μₐ), β=$(par.β), sims=$(par.num_sims)",
    titlelocation = :center,
    titlefontsize = 11,
    titlefontcolor = :black,
    top_margin=8Plots.mm,
    size=(600, 500)
)
display(p)
savefig(p,"Plots/Others/speed_AES_MILP_BF_sigmaV$(par.σᵥ)_$(CPU).pdf")

# Output the results to a LaTeX table
import DataFrames: groupby
groups = groupby(df[:, [:N, :L, :avg_time_milp, :avg_time_aes, :loc_gap_milp_aes]], :N)
df_tex = hcat([select(g, Not(:N)) for g in groups]..., makeunique=true)
# drop L_1, L_2 columns 
select!(df_tex, Not(r"L_\d+"))

write_to_tex(df_tex, "Tables/Others/speed_AES_MILP_BF_sigmaV$(par.σᵥ)_$(CPU).tex", rounding = [0, 4, 4, 1, 4, 4, 1])
