# Functions for grid generation and manipulation

# # Load functions
# include("system.jl")

# Construct a Cartesian grid
function gridGen(xL, xR, imax)
    # Build the grid
    dx = (xR - xL) / imax;
    x = collect(xL:dx:xR);
    return x
end
