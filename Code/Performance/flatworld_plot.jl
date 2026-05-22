function plot_solution(L, N, config, active_y, Xc, Yc, Xf, Yf, fc, A)
    K = length(L) # number of stages

    # Plot dictionary
    if K == 1
        color = [:red, :blue]
        label = ["Assembly", "Consumers"]
        shape = [:square, :circle]
    elseif K == 2
        color = [:green, :red, :blue]
        label = ["Battery", "Assembly", "Consumers"]
        shape = [:utriangle, :square, :circle]
    elseif K == 3
        color = [:black, :green, :red, :blue]
        label = ["Mineral", "Battery", "Assembly", "Consumers"]
        shape = [:diamond, :utriangle, :square, :circle]
    end
    lb = 3; ub = 8 # lower bound and upper bound for marker size

    # Display potential clients [destination markets]
    # rescale A to be between lb and ub for better visualization
    ms_n = [
        maximum(v) == minimum(v) ? fill((lb+ub)/2, length(v)) :
        lb .+ (ub-lb) .* (v .- minimum(v)) ./ (maximum(v) - minimum(v))
        for v in A
    ]
    p = scatter(Xc, Yc, label=label[K+1], markercolor=:white, markersize=ms_n, markershape=shape[K+1], markerstrokecolor=color[K+1], markerstrokewidth=0.5, xticks=0:0.25:1, yticks=0:0.25:1, xlims=(0,1), ylims=(0,1))

    # Display potential locations [production plants]
    # rescale fc to be between lb and ub for better visualization
    ms_l = [
        maximum(v) == minimum(v) ? fill((lb+ub)/2, length(v)) :
        lb .+ (ub-lb) .* (v .- minimum(v)) ./ (maximum(v) - minimum(v))
        for v in fc
    ]
    for k in K:-1:1 # to make legend in the right order
        scatter!(Xf[k], Yf[k], label=label[k], markercolor=:white, markersize=ms_l[k], markershape=shape[k],
        markerstrokecolor=color[k], markerstrokewidth=0.5)
    end

    # Show markets with positive sales 
    nonzero_q = sum(config[1],dims = 2) .> 0;
    active_N = sum(nonzero_q); # the number of markets entered
    mc = [(nonzero_q[n] ? color[K+1] : :white) for n in 1:N]
    scatter!(Xc, Yc, 
            markershape=shape[K+1], markercolor=mc, markersize=ms_n,
            markerstrokecolor=color[K+1], markerstrokewidth=0.5,
            label=nothing
        )

    # Show open facility
    for k in 1:K
        y_ = zeros(Int, L[k])
        y_[active_y[k]] .= 1
        mc = [(y_[l] == 1 ? color[k] : :white) for l in 1:L[k]]
        scatter!(Xf[k], Yf[k], 
            markershape=shape[k], markercolor=mc, markersize=ms_l[k],
            markerstrokecolor=color[k], markerstrokewidth=0.5,
            label=nothing
        )
    end
    # the number of locations produced
    active_L = [length(active_y[k]) for k in 1:K]

    # Show client-facility assignment
    for n in 1:N
        if nonzero_q[n] > 0
            for k in 1:K-1
                plot!([Xf[k][config[1][n,k]], Xf[k+1][config[1][n,k+1]]], [Yf[k][config[1][n,k]], Yf[k+1][config[1][n,k+1]]], color=color[k], label=nothing) 
            end
            plot!([Xc[n], Xf[K][config[1][n,K]]], [Yc[n], Yf[K][config[1][n,K]]], color=color[K], label=nothing)
        end
    end

    # Some graph formatting, legend, title, size
    scatter!(legend=:outerbottom, legendcolumns=K+1)
    plot!(title="Nsold=$active_N, Lproduced=$active_L",title_position=:left,subplot=1) 

    plot!(size=(800,800)) # adjust the size of the plot

    return p
end