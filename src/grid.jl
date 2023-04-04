# Functions for grid generation and manipulation

# Read the domain information from a file and build the grid
function makeGrid(filename)
    # Read the domain information from a file
    file = open(filename, "r");
    xL = parse(Float64,split(readline(file), "#")[1]);
    xR = parse(Float64,split(readline(file), "#")[1]);
    imax = parse(Int,split(readline(file), "#")[1]);
    close(file)
    # Build the grid
    x = gridGen(xL, xR, imax);

    return x, imax
end

# Construct a Cartesian grid
function gridGen(xL, xR, imax)
    # Build the grid
    dx = (xR - xL) / imax;
    x = collect(xL:dx:xR);
    return x
end
