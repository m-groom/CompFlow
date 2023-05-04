# Functions for grid generation and manipulation

# Read the domain information from a file and build the grid
function makeGrid(filename)
    # Initialise boundary condition array
    BCs = zeros(3, 2)
    # Read the domain information from a file
    report("Reading the domain information from file $(filename)")
    file = open(filename, "r")
    xL = parse(Float64, split(readline(file), "#")[1])
    xR = parse(Float64, split(readline(file), "#")[1])
    imax = parse(Int, split(readline(file), "#")[1])
    BCs[1, 1] = parse(Int, split(readline(file), "#")[1])
    BCs[1, 2] = parse(Int, split(readline(file), "#")[1])
    close(file)
    # Build the grid
    report("Building the grid")
    x = gridGen(xL, xR, imax)
    report("Number of cells: $(length(x)-1)")

    return x, BCs
end

# Construct a Cartesian grid
function gridGen(xL, xR, imax)
    # Build the grid
    dx = (xR - xL) / imax
    x = collect(xL:dx:xR)
    return x
end
