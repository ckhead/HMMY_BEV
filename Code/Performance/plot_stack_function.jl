using Plots
"""
Plot y vs x with optional stacked composition bars above each data point.

# Arguments
- `x_vector`: A vector with equal increments (e.g., 0:0.1:1 or collect(0:0.1:1))
- `y_vector`: A vector of y values corresponding to x_vector
- `compositions`: Optional vector of (c, u, o) tuples where c + u + o = 1 for each x value.
                  When provided, draws stacked bars above the plot showing the composition.


# Keyword Arguments
- `title=""`: Plot title
- `xlabel="x"`: Label for x-axis
- `ylabel="y"`: Label for y-axis
- `linewidth=2`: Width of the plot line
- `marker=:circle`: Marker style for data points
- `markersize=4`: Size of markers
- `c_color=:white`: Color for the bottom (c) section of composition bars
- `u_color=:gray`: Color for the middle (u) section of composition bars
- `o_color=:black`: Color for the top (o) section of composition bars
- `yscale=:linear`: Scale for y-axis, either :linear or :log10
- `kwargs...`: Additional keyword arguments passed to `plot()`

# Returns
- A `Plots.Plot` object
 
# Example 
```julia
x = 0:0.1:1
y = [6, 0.5, 0.1, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 3e-4, 2e-4, 1e-4] # Example y values spanning several orders of magnitude

# Compositions: vector of (c, u, o) tuples for each x, from bottom to top
compositions = [
    (0.5, 0.3, 0.2),  # x = 0.0
    (0.3, 0.4, 0.3),  # x = 0.1
    (0.2, 0.5, 0.3),  # x = 0.2
    (0.4, 0.2, 0.4),  # x = 0.3
    (0.6, 0.2, 0.2),  # x = 0.4
    (0.3, 0.3, 0.4),  # x = 1.5
    (0.2, 0.4, 0.4),  # x = 0.6
    (0.4, 0.3, 0.3),  # x = 0.7
    (0.5, 0.2, 0.3),  # x = 0.8
    (0.3, 0.5, 0.2),  # x = 0.9
    (0.4, 0.4, 0.2)   # x = 1.0
]

p = plot_stack(x, y, compositions;
    title="Plot with Composition Bars",
    xlabel="x",
    ylabel="y",
    c_color=:white,
    u_color=:gray,
    o_color=:black)


p = plot_stack(x, y, compositions;
    title="Plot with Composition Bars",
    xlabel="x",
    ylabel="y", yscale=:log10,
    c_color=:white,
    u_color=:gray,
    o_color=:black)
```
"""
function plot_stack(x_vector, y1_vector, y2_vector, compositions=nothing; 
                 title="", 
                 xlabel="x", 
                 ylabel="y",
                 linewidth=2,
                 marker=:circle,
                 markersize=4,
                 c_color=:white,
                 u_color=:gray,
                 o_color=:black,
                 yscale=:identity,
                 kwargs...)
    
    # Convert to array if it's a range
    x = collect(x_vector)
    y1 = collect(y1_vector)
    y2 = collect(y2_vector)
    y = [y1; y2]
    
    # Calculate the increment (difference between consecutive x values)
    dx = x[2] - x[1]
    
    # Calculate axis limits
    x_min = minimum(x) - 0.8 * dx
    x_max = maximum(x) + 0.8 * dx
    y_min = minimum(y) * 0.8
    y_data_max = maximum(y)
    
    if yscale == :identity
        y_max = y_min + (y_data_max - y_min) * 1.25
    elseif yscale == :log10
        y_max = y_min * 10^(log10(y_data_max / y_min) * 1.25)
    end
    
    # Create the plot
    p = plot(x, y1;
             xlims=(x_min, x_max),
             ylims=(y_min, y_max),
             yscale=yscale,
             size=(600, 600),
             title=title,
             xlabel=xlabel,
             ylabel=ylabel,
             linewidth=linewidth,
             marker=marker,
             markersize=markersize,
             legend=false,
             kwargs...)
    
    # Add composition bars if provided
    if compositions !== nothing
        if yscale == :identity
            bar_bottom = y_min + (y_data_max - y_min) * 1.1
            bar_top = y_min + (y_data_max - y_min) * 1.3
            bar_height = bar_top - bar_bottom
        elseif yscale == :log10
            # Work in log space
            log_y_min = log10(y_min)
            log_y_data_max = log10(y_data_max)
            log_range = log_y_data_max - log_y_min
            
            log_bar_bottom = log_y_min + log_range * 1.05
            log_bar_top = log_y_min + log_range * 1.25
            log_bar_height = log_bar_top - log_bar_bottom
        end
        
        for (i, xi) in enumerate(x)
            c, u, o = compositions[i]
            
            # Bar boundaries
            x_left = xi - 0.4 * dx
            x_right = xi + 0.4 * dx
            
            if yscale == :identity
                # Linear: work in normal space
                c_bottom = bar_bottom
                c_top = bar_bottom + c * bar_height
                
                u_bottom = c_top
                u_top = u_bottom + u * bar_height
                
                o_bottom = u_top
                o_top = bar_top
            else
                # Log: work in log space, then convert back
                log_c_bottom = log_bar_bottom
                log_c_top = log_bar_bottom + c * log_bar_height
                
                log_u_bottom = log_c_top
                log_u_top = log_u_bottom + u * log_bar_height
                
                log_o_bottom = log_u_top
                log_o_top = log_bar_top
                
                # Convert back to linear for plotting
                c_bottom = 10^log_c_bottom
                c_top = 10^log_c_top
                u_bottom = 10^log_u_bottom
                u_top = 10^log_u_top
                o_bottom = 10^log_o_bottom
                o_top = 10^log_o_top
            end
            
            # Draw rectangles using Shape
            plot!(p, Shape([x_left, x_right, x_right, x_left], 
                          [c_bottom, c_bottom, c_top, c_top]),
                  fillcolor=c_color, linecolor=:black, linewidth=0.5)
            
            plot!(p, Shape([x_left, x_right, x_right, x_left], 
                          [u_bottom, u_bottom, u_top, u_top]),
                  fillcolor=u_color, linecolor=:black, linewidth=0.5)
            
            plot!(p, Shape([x_left, x_right, x_right, x_left], 
                          [o_bottom, o_bottom, o_top, o_top]),
                  fillcolor=o_color, linecolor=:black, linewidth=0.5)
        end
    end
    
    return p
end


