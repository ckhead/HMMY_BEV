using Plots
using JLD2
using DataFrames

include("speed_setup_par.jl")
include("plot_stack_function.jl")
getCPU() = join(split(Sys.cpu_info()[1].model)[2:end])
CPU = getCPU()

# par = ModelParams(σᵥ=0.25, β=0.25, σᵩ=1.0, μᵩ=-1.0, μₐ = 2.0, num_sims = 100) 
# par = ModelParams(σᵥ=0.48, β=0.25, σᵩ=1.0, μᵩ=-2.0, num_sims = 100) 
par = ModelParams(β = 0.25, σᵩ = 1.5, μₐ = 4.5, μₙ = 0.0, num_sims = 300)

# Load input files 
L_list = [20, 40, 80]
L = 40; N= 40;
σᵤ_grid = round.(collect(0.0:0.1:1), digits=1)
σᵥ_grid = round.(collect(0.0:0.1:1), digits=1)
σᵩ_grid = round.(collect(0.0:0.2:2), digits=1)
β_grid = round.(collect(0.0:0.1:1), digits=1)

file_path = "Data/JLD/MILP_CDCP_VarSigmaV_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
dfV =  load(file_path, "v")
Var_σᵥ = dfV[(dfV.L .== L) .& (dfV.N .== N) .&
        (round.(dfV.σᵥ, digits=2) .∈ Ref(σᵥ_grid)), :]
rename!(Var_σᵥ, :avg_N_nodes => :nodes, :avg_N_undet => :undet, :avg_N_open => :open, :avg_N_close => :closed, :avg_N_plants => :plants)

file_path = "Data/JLD/MILP_CDCP_VarSigmaU_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
dfU =  load(file_path, "u")
Var_σᵤ = dfU[(dfU.L .== L) .& (dfU.N .== N) .&
        (round.(dfU.σᵤ, digits=2) .∈ Ref(σᵤ_grid)), :]
rename!(Var_σᵤ, :avg_N_nodes => :nodes, :avg_N_undet => :undet, :avg_N_open => :open, :avg_N_close => :closed, :avg_N_plants => :plants)

file_path = "Data/JLD/MILP_CDCP_VarSigmaF_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
dfF =  load(file_path, "phi")
Var_σᵩ = dfF[(dfF.L .== L) .& (dfF.N .== N) .&
        (round.(dfF.σᵩ, digits=2) .∈ Ref(σᵩ_grid)), :]
rename!(Var_σᵩ, :avg_N_nodes => :nodes, :avg_N_undet => :undet, :avg_N_open => :open, :avg_N_close => :closed, :avg_N_plants => :plants)


file_path = "Data/JLD/MILP_CDCP_VarBeta_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
dfB =  load(file_path, "B")
Var_β = dfB[(dfB.L .== L) .& (dfB.N .== N) .&
        (round.(dfB.β, digits=2) .∈ Ref(β_grid)), :]
rename!(Var_β, :avg_N_nodes => :nodes, :avg_N_undet => :undet, :avg_N_open => :open, :avg_N_close => :closed, :avg_N_plants => :plants)


## Facility variable cost std. dev. σᵥ
pmarg = [1.35, 0.6, 0.8, 1.05, 1.01]
Var_σᵥ.compositions = tuple.(Var_σᵥ.closed/L, Var_σᵥ.undet/L, Var_σᵥ.open/L)
maxy = maximum(vcat(Var_σᵥ.t_aes, Var_σᵥ.t_milp))
ylim = maxy * pmarg[1]
xlim = maximum(Var_σᵥ.σᵥ) 
xstep = Var_σᵥ.σᵥ[2] - Var_σᵥ.σᵥ[1]

plot_stack(Var_σᵥ.σᵥ, Var_σᵥ.t_aes, Var_σᵥ.t_milp, Var_σᵥ.compositions;
    xlabel="Std. dev. of variable facility costs (σᵥ)",color=:orange,
    ylabel="Time in seconds", # yscale=:log10,
    c_color=:white,
    u_color=:gray,
    o_color=:black,
    ylims = (0, ylim), xlims = (-xstep, xlim+xstep))
plot!(Var_σᵥ.σᵥ, Var_σᵥ.t_milp, label = "MILP", color = :blue,marker=:diamond, ms=4)
annotate!(xlim, Var_σᵥ.t_aes[end]+ylim/20, text("AES", :left, :orange, 11))
annotate!(xlim, Var_σᵥ.t_milp[end]+ylim/20, text("MILP", :left, :blue, 11))
# add number of nodes annotation for MILP
annotate!(xlim*0.6, Var_σᵥ.t_milp[end]*pmarg[2], text("#nodes", :center, :blue, 11))
annotate!(Var_σᵥ.σᵥ, Var_σᵥ.t_milp*pmarg[3], [text(s, :center, :blue, 8) for s in string.(round.(Var_σᵥ.nodes; digits=1))])
# add legend for open, undetermined, closed on top
annotate!(0.0,  ylim, text("\u25A0 Open", :left, :black, 11))
annotate!(xlim/2.5, ylim, text("\u25A0 Undetermined", :left, RGB(0.5, 0.5, 0.5), 11))
annotate!(xlim/1.2,  ylim, text("\u25A1 Closed",:left, :black, 11))
# add number of plants below the blocks bar 
annotate!(Var_σᵥ.σᵥ, maxy*pmarg[4], [text(s, :center, :black, 11) for s in string.(round.(Var_σᵥ.plants; digits=1))])
annotate!(xlim, maxy/pmarg[5], text("#plants", :center, :black, 11))
# add parameters as title on top
plot!(
    title = "σᵤ=$(par.σᵤ), β=$(par.β), σᵩ=$(par.σᵩ), μᵩ=$(par.μᵩ), μₐ=$(par.μₐ), sims=$(par.num_sims)",
    titlelocation = :center,
    titlefontsize = 11,
    titlefontcolor = :black
)

savefig("Plots/Others/speed_UFLP_VarSigmaV_L$(L)_N$(N)_$(CPU).pdf")    

## Path cost std. dev. σᵤ
pmarg = [1.35, 0.6, 0.8, 1.05, 1.01]
Var_σᵤ.compositions = tuple.(Var_σᵤ.closed/L, Var_σᵤ.undet/L, Var_σᵤ.open/L)
maxy = maximum(vcat(Var_σᵤ.t_aes, Var_σᵤ.t_milp))
ylim = maxy * pmarg[1]
xlim = maximum(Var_σᵤ.σᵤ)
xstep = Var_σᵤ.σᵤ[2] - Var_σᵤ.σᵤ[1]

plot_stack(Var_σᵤ.σᵤ, Var_σᵤ.t_aes, Var_σᵤ.t_milp, Var_σᵤ.compositions;
    xlabel="Std. dev. of path costs (σᵤ)",color=:orange,
    ylabel="Time in seconds", # yscale=:log10,
    c_color=:white,
    u_color=:gray,
    o_color=:black,
    ylims = (0, ylim), xlims = (-xstep, xlim+xstep))
plot!(Var_σᵤ.σᵤ, Var_σᵤ.t_milp, label = "MILP", color = :blue,marker=:diamond, ms=4)
annotate!(xlim, Var_σᵤ.t_aes[end]+ylim/20, text("AES", :left, :orange, 11))
annotate!(xlim,Var_σᵤ.t_milp[end]+ylim/20, text("MILP", :left, :blue, 11))
annotate!(xlim*0.6, Var_σᵤ.t_milp[end]*pmarg[2], text("#nodes", :center, :blue, 11))
annotate!(Var_σᵤ.σᵤ, Var_σᵤ.t_milp*pmarg[3], [text(s, :center, :blue, 8) for s in string.(round.(Var_σᵤ.nodes; digits=1))])
annotate!(0.0,  ylim, text("\u25A0 Open", :left, :black, 11))
annotate!(xlim/2.5, ylim, text("\u25A0 Undetermined", :left, RGB(0.5, 0.5, 0.5), 11))
annotate!(xlim/1.2,  ylim, text("\u25A1 Closed", :left, :black, 11))
# add number of plants below the blocks bar 
annotate!(Var_σᵤ.σᵤ, maxy*pmarg[4], [text(s, :center, :black, 11) for s in string.(round.(Var_σᵤ.plants; digits=1))])
annotate!(xlim, maxy/pmarg[5], text("#plants", :center, :black, 11))
# add parameters as title on top
plot!(
    title = "σᵥ=$(par.σᵥ), β=$(par.β), σᵩ=$(par.σᵩ), μᵩ=$(par.μᵩ), μₐ=$(par.μₐ), sims=$(par.num_sims)",
    titlelocation = :center,
    titlefontsize = 11,
    titlefontcolor = :black
)

savefig("Plots/Others/speed_UFLP_VarSigmaU_L$(L)_N$(N)_$(CPU).pdf")    


## Facility fixed cost std. dev. σᵩ
pmarg = [1.32, 0.55, 0.8, 1.05, 1.01]
Var_σᵩ.compositions = tuple.(Var_σᵩ.closed/L, Var_σᵩ.undet/L, Var_σᵩ.open/L)
maxy = maximum(vcat(Var_σᵩ.t_aes, Var_σᵩ.t_milp))
ylim = maxy * pmarg[1]
xlim = maximum(Var_σᵩ.σᵩ)
xstep = Var_σᵩ.σᵩ[2] - Var_σᵩ.σᵩ[1]

plot_stack(Var_σᵩ.σᵩ, Var_σᵩ.t_aes, Var_σᵩ.t_milp, Var_σᵩ.compositions;
    xlabel="Std. dev. of fixed costs (σᵩ)",color=:orange,
    ylabel="Time in seconds", # yscale=:log10,
    c_color=:white,
    u_color=:gray,
    o_color=:black,
    ylims = (0, ylim), xlims = (-xstep, xlim+xstep))
plot!(Var_σᵩ.σᵩ, Var_σᵩ.t_milp, label = "MILP", color = :blue,marker=:diamond, ms=4)
annotate!(xlim, Var_σᵩ.t_aes[end]+ylim/20, text("AES", :left, :orange, 11))
annotate!(xlim,Var_σᵩ.t_milp[end]+ylim/20, text("MILP", :left, :blue, 11))
annotate!(xlim*0.6, Var_σᵩ.t_milp[end]*pmarg[2], text("#nodes", :center, :blue, 11))
annotate!(Var_σᵩ.σᵩ, Var_σᵩ.t_milp*pmarg[3], [text(s, :center, :blue, 8) for s in string.(round.(Var_σᵩ.nodes; digits=1))])
annotate!(0.0,  ylim, text("\u25A0 Open", :left, :black, 11))
annotate!(xlim/2.5, ylim, text("\u25A0 Undetermined", :left, RGB(0.5, 0.5, 0.5), 11))
annotate!(xlim/1.2,  ylim, text("\u25A1 Closed", :left, :black, 11))
# add number of plants below the blocks bar 
annotate!(Var_σᵩ.σᵩ, maxy*pmarg[4], [text(s, :center, :black, 11) for s in string.(round.(Var_σᵩ.plants; digits=1))])
annotate!(xlim, maxy/pmarg[5], text("#plants", :center, :black, 11))
# add parameters as title on top
plot!(
    title = "σᵥ=$(par.σᵥ), σᵤ=$(par.σᵤ), β=$(par.β), μᵩ=$(par.μᵩ), μₐ=$(par.μₐ), sims=$(par.num_sims)",
    titlelocation = :center,
    titlefontsize = 11,
    titlefontcolor = :black
)

savefig("Plots/Others/speed_UFLP_VarSigmaPhi_L$(L)_N$(N)_$(CPU).pdf")  

## Distance cost elasticity β
pmarg = [1.32, 0.6, 0.8, 1.05, 1.01]
Var_β.compositions = tuple.(Var_β.closed/L, Var_β.undet/L, Var_β.open/L)
maxy = maximum(vcat(Var_β.t_aes, Var_β.t_milp))
ylim = maxy * pmarg[1]
xlim = maximum(Var_β.β) 
xstep = Var_β.β[2] - Var_β.β[1]

plot_stack(Var_β.β, Var_β.t_aes, Var_β.t_milp, Var_β.compositions;
    xlabel="Distance cost elasticity (β)",color=:orange,
    ylabel="Time in seconds", # yscale=:log10,
    c_color=:white,
    u_color=:gray,
    o_color=:black,
    ylims = (0, ylim), xlims = (-xstep, xlim+xstep))
plot!(Var_β.β, Var_β.t_milp, label = "MILP", color = :blue,marker=:diamond, ms=4)
annotate!(xlim, Var_β.t_aes[end]+ylim/20, text("AES", :left, :orange, 11))
annotate!(xlim, Var_β.t_milp[end]+ylim/20, text("MILP", :left, :blue, 11))
annotate!(xlim*0.6, Var_β.t_milp[end]*pmarg[2], text("#nodes", :center, :blue, 11))
annotate!(Var_β.β, Var_β.t_milp*pmarg[3], [text(s, :center, :blue, 8) for s in string.(round.(Var_β.nodes; digits=1))])
annotate!(0.0,  ylim, text("\u25A0 Open",         :left, :black,             11))
annotate!(xlim/2.5, ylim, text("\u25A0 Undetermined", :left, RGB(0.5, 0.5, 0.5), 11))
annotate!(xlim/1.2,  ylim, text("\u25A1 Closed",       :left, :black,            11))
# add number of plants below the blocks bar 
annotate!(Var_β.β, maxy*pmarg[4], [text(s, :center, :black, 11) for s in string.(round.(Var_β.plants; digits=1))])
annotate!(xlim, maxy/pmarg[5], text("#plants", :center, :black, 11))
# add parameters as title on top
plot!(
    title = "σᵥ=$(par.σᵥ), σᵤ=$(par.σᵤ), σᵩ=$(par.σᵩ), μᵩ=$(par.μᵩ), μₐ=$(par.μₐ), sims=$(par.num_sims)",
    titlelocation = :center,
    titlefontsize = 11,
    titlefontcolor = :black
)

savefig("Plots/Others/speed_UFLP_VarBeta_L$(L)_N$(N)_$(CPU).pdf")
