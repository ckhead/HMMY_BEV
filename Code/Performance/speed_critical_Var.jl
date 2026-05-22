using DataFrames
using JLD2
#load the dataframes
# Import auxiliary function for writing LaTeX tables
include("write_to_tex_function.jl")
getCPU() = join(split(Sys.cpu_info()[1].model)[2:end])
CPU = getCPU()

L_list = [20, 40, 80]

file_path = "Data/JLD/MILP_CDCP_VarSigmaU_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
dfU =  load(file_path, "u")
file_path = "Data/JLD/MILP_CDCP_VarSigmaF_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
dfF =  load(file_path, "phi")
file_path = "Data/JLD/MILP_CDCP_VarSigmaV_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
dfV =  load(file_path, "v")
file_path = "Data/JLD/MILP_CDCP_VarBeta_L_$(L_list[1])_$(L_list[end])_$(CPU).jld2"
dfB =  load(file_path, "B")

rename!(dfU,:σᵤ => :sigma);
rename!(dfV,:σᵥ => :sigma);
rename!(dfF,:σᵩ => :sigma);
rename!(dfB,:β => :sigma);
dfU.sigma_type .= "Path costs (\$\\sigma^u \\in \\{0, 1\\} \$)"
dfV.sigma_type .= "Variable facility costs (\$\\sigma^v \\in \\{0, 1\\} \$)"
dfF.sigma_type .= "Fixed costs (\$\\sigma^\\phi \\in \\{0, 2\\} \$)"
dfB.sigma_type .= "Distance cost elasticity (\$\\beta \\in \\{0, 1\\} \$)"

df = vcat(dfV, dfU) #, dfF, dfB)

# Find critical sigma with interpolation
grouped = DataFrames.groupby(df, [:L, :sigma_type])

critical_sigmas = combine(grouped) do subdf
    sort!(subdf, :sigma)
    
    # Find the transition point assuming AES decreases with sigma values
    diff_vals = subdf.t_aes .- subdf.t_milp
    idx = findfirst(diff_vals .< 0)
    
    if isnothing(idx)
        # MILP dominates until the end of parameter range
        # return DataFrame(critical_sigma = subdf.sigma[end]) 
        return DataFrame(critical_sigma = "-")
    elseif idx == 1
        # AES dominates from the first parameter value
        return DataFrame(critical_sigma = subdf.sigma[1])
    else
        # There is a crossover, and hence linear interpolation between points
        σ1, σ2 = subdf.sigma[idx-1], subdf.sigma[idx]
        d1, d2 = diff_vals[idx-1], diff_vals[idx]
        
        # Find where the line crosses zero
        critical_σ = σ1 - d1 * (σ2 - σ1) / (d2 - d1)
        
        return DataFrame(critical_sigma = critical_σ)
    end
end

# Reshape the table
result_table = unstack(critical_sigmas, :sigma_type, :L, :critical_sigma)

# Output the results to a LaTeX table
write_to_tex(result_table, "Tables/Others/speed_critical_sigmas_$(CPU).tex", rounding = [0, 2, 2, 2])