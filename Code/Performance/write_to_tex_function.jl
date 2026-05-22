using Printf
using DataFrames

"""
    write_to_tex(data, [filepath]; rounding, col_select)

Write a DataFrame or Matrix to a LaTeX-formatted table (rows separated by `&` and `\\\\`). Packages "DataFrames" and "Printf" are required.

# Arguments
- `data::Any`: a `DataFrame` or `Matrix` (matrices are converted to DataFrames automatically).
- `filepath::Union{Nothing, String}=nothing`: output file path. If `nothing` (default), prints to console.

# Keyword Arguments
- `rounding::Union{Nothing, Vector{Int}}=nothing`: number of decimal places for each column.
  Must match the number of columns *after* `col_select` is applied. `0` rounds to integer.
  Defaults to `[0, 0, ...]` (all integers) if omitted.
- `col_select::Union{Nothing, Vector{Symbol}, Vector{Int}}=nothing`: subset of columns to include,
  specified by name (`Symbol`) or position (`Int`).

# Examples
```julia
data = DataFrame(
    A = 1:5,
    B = [1.1234, 2.1234, 3.1234, 4.1234, 5.1234],
    C = [1000.55, 2000.75, 3000.95, 4000.15, 5000.35],
    D = [1, 2, 3, 4, 5],
    E = ["a", "b", "c", "d", "e"]
)
# Print to console with 2 decimal places for columns A and B, 1 for C
write_to_tex(data, rounding=[2, 2, 1], col_select=[:A, :B, :C])

# Write to file
write_to_tex(data, "output.tex", rounding=[0, 2, 1])

my_par_est = DataFrame(
       parameter = ["\$\\alpha\$","\$\\beta\$","\$\\gamma\$"],
       description = [" ", " " ," "],
       estimate = [0.32,0.3131,2]
)
write_to_tex(my_par_est, "my_par_est.tex", rounding=[0, 0, 2])

"""
function write_to_tex(data::Any, filepath::Union{Nothing, String}=nothing; rounding::Union{Nothing, Vector{Int}}=nothing, col_select::Union{Nothing, Vector{Symbol}, Vector{Int}} = nothing)
    if isa(data, Matrix)
        data = DataFrame(data, :auto)
    end

    if !isnothing(col_select)
        # Select the desired columns if specified
        data = select(data, col_select)
    end

    if isnothing(rounding)
        # Default to no rounding if not specified
        rounding = fill(0, ncol(data))
    end

    # Ensure the rounding vector matches the number of columns
    if length(rounding) != ncol(data)
        error("The length of the rounding vector must match the number of columns in the data.")
    end

    function format_rows(io)
        for row in eachrow(data)
            formatted_row = [
                if isa(row[i], Number)
                   rounding[i] == 0 ? round(Int, row[i]) : @sprintf("%.*f", rounding[i], row[i])
                   # rounding[i] == 0 ? round(Int, row[i]) : round(row[i], digits=rounding[i])
                else
                    row[i]
                end
                for i in 1:ncol(data)
            ]
            println(io, join(formatted_row, " & "), " \\\\")
        end
    end

    if isnothing(filepath)
        format_rows(stdout)
    else
        open(filepath, "w") do io
            format_rows(io)
        end
    end
end
