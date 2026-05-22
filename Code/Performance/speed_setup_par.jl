using Random
using LinearAlgebra
using Distributions 
using Parameters
"""
core parameters 
- arranged in a non-mutable keyword structure, 
- must be using Parameters.jl
- seed_start is the starting seed for random number generation
- num_sims is the number of simulations to run
- milp_threads is the number of threads to use in Gurobi for MILP, defaults to 1
- beta is the distance elasticity
- eta is the elasticity of substitution between varieties 
- eta >1 to ensure finite markups
- σᵤ is the standard deviation of the normal path shock to marginal costs
- σᵥ is the standard deviation of the log productivity shocks z
- we refer σᵤ=σᵥ=0 as the Euclidean case since location is the only source of randomness
- σᵩ is the standard deviation of  log fixed costs
- μᵩ is the mean of the log fixed costs
- set μᵩ=σᵩ = 0 for fixed costs = 1 for all locations
"""
@with_kw struct ModelParams
   seed_start::Int = 0
   num_sims::Int = 40
   milp_threads::Int = 1
   η::Float64 = 4.0 
   σᵤ::Float64 = 0.0
   # For the rest of params, assumed to be the same across stages in the multi-stage model
   β::Float64 = 0.25 
   σᵥ::Float64 = 0.0 
   μᵥ::Float64 = 1.0 # 0.0
   σᵩ::Float64 = 0.0 
   μᵩ::Float64 = 0.0 
   μₙ::Float64 = 1.0 # fixed cost of market entry
   # heterogeneity in demand shock A_n 
   σₐ::Float64 = 0.0
   μₐ::Float64 = 0.0
   @assert η > 1
   @assert β >= 0
   @assert σᵤ>=0
   @assert σᵥ>=0
   @assert σᵩ>=0
   @assert σₐ>=0
end

"""
Setup creates the data for a single-stage production model.
Inputs:
   - L: Number of potential facility locations
   - N: Number of consumers
   - par: Model parameters (ModelParams struct)
   - Mf: number of models per firm, defaults to 1
Outputs:
   - fc: Fixed costs vector (L x 1) for each facility location
   - mc: Marginal cost matrix (N x L) for consumers and facilities
   - fs: Fixed cost of market entry (N x 1)for each consumer
   - A: Demand shock (Mf x N) for each model and consumer
   - Xc, Yc: Consumers' locations
   - Xf, Yf: Facilities' potential locations
Additional features in a multi-stage model compared to single stage:
   - market entry decision z 
   - fixed cost of market entry, fs 
"""
function setup_par_K(L, N, par::ModelParams; Mf = 1)
   K = length(L) # number of stages
   @unpack η, β, σᵤ, σᵥ, μᵥ, σᵩ, μᵩ, μₙ, σₐ, μₐ = par
   # cost share 
   α = [1 / k for k in 1:K] # α: each element is α_kk
   γ = [prod(1 .- α[(k+1):K]) for k in 1:(K-1)] # assume CRS

   # Consumers' locations, X = lng , Y = lat
   Xc = rand(N);
   Yc = rand(N); 

   # Facilities' potential locations
   Xf = [];
   Yf = []; 
   # Variable costs
   dbn_v = Normal(μᵥ, σᵥ)
   v = []; 
   # Fixed costs
   dbn_ϕ = LogNormal(μᵩ, σᵩ)
   fc = []; 
   for k in K:-1:1 
      pushfirst!(Xf, rand(L[k])) # push in reverse order
      pushfirst!(Yf, rand(L[k]))
      pushfirst!(v, rand(dbn_v, L[k]))
      pushfirst!(fc, rand(dbn_ϕ, L[k]))
   end
   # normalize when there are multiple stages of production
   fc = fc ./ K

   # Fixed costs of sale / distribution
   # fs = rand(dbn_ϕ, N)
   fs = ones(N) .* μₙ # set fixed cost of market entry to be the same across consumers

   # normally distributed path shocks
   dbn_u = Normal(0, σᵤ)
   u = rand(dbn_u, Mf, N, L...);

   # Trade costs
   # Distance between level-(k+1) and level-k facilities
   home_d = 1.0 # home distance, 1.0 or can be set to 0 for no internal geography cost
   d_all = [];
   for k in 1:K
      if k < K
         d = [home_d + norm([Xf[k][l_u] - Xf[k+1][l_d], Yf[k][l_u] - Yf[k+1][l_d]], 2) for l_u in 1:L[k], l_d in 1:L[k+1]]
      elseif k == K
         d = [home_d + norm([Xf[k][l] - Xc[n], Yf[k][l] - Yc[n]], 2) for l in 1:L[k], n in 1:N] # last stage is distance to consumers
      end
      push!(d_all, d)
   end

   ## Path costs
   #  log cost for each stage k
   lc_all = Vector{Array{Float64}}(undef, K)
   for k in 1:K
      if k < K
         lc_all[k] = γ[k] .* (α[k] .* v[k] .+ β .* log.(d_all[k]))  # (L[k], L[k+1])
      elseif k == K
         lc_all[K] = α[K] .* v[K] .+ β .* log.(d_all[K])              # (L[K], N)
      end
   end

   #  total log cost as array of shape (N, L...)
   lmc = zeros(Float64, N, L...);
   for k in 1:K
      if k < K
         # add upstream costs k = 1,...,K-1
         shp = ones(Int, K + 1)     # dimensions are (N, L[1], ..., L[K])
         shp[k+1] = L[k]
         shp[k+2] = L[k+1]
         lmc .+= reshape(lc_all[k], Tuple(shp)) # reshape to fit into the dimensions
      elseif k == K
         # add final cost at stage K
         shp = ones(Int, K + 1) # dimensions are (N, L[1], ..., L[K])
         shp[1] = N
         shp[K+1] = L[K]
         lmc .+= reshape(lc_all[K]', Tuple(shp)) # reshape to fit into the dimensions
      end
   end
   # repeat lmc for every model m = 1,...,Mf, make it of shape (Mf, N, L...)
   lmc = repeat(reshape(lmc, 1, N, L...), Mf, ones(Int, K+1)...)

   # add path shock and exponentiate
   mc = exp.(lmc .+ u);

   # add heterogeneity in demand shock A_n
   dbn_a = LogNormal(μₐ, σₐ)
   A = rand(dbn_a, Mf, N)

   return fc, mc, fs, A, Xc, Yc, Xf, Yf
end

"""
## Compute variable profits given marginal costs and model parameters.
- Assumes CES demand with elasticity η
- mc: Marginal cost matrix (Mf x N x L...)
- A: Demand shifters (Mf x N)
- par: Model parameters (ModelParams struct)
- Returns a variable profits matrix (Mf x N x L...)
"""
function vprofit_fn(mc, A, par::ModelParams)
    @unpack η = par
    markup = η/(η-1)
    p = markup .* mc
    q = A .* p.^(-η)
    vprofit = q .* (p .- mc)
    return vprofit
end