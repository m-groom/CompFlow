# Functions for computing boundary conditions

# Function for computing slopes at a boundary
function slopeBoundaryCondition(Q, BC, side)
    if (BC == 0) # Periodic
        error("Periodic boundary conditions not implemented yet")
    elseif (BC == 1) # Extrapolate
        return zeros(size(Q))
    elseif (BC == 2) # Reflective
        error("Reflective boundary conditions not implemented yet")
    elseif (BC == 3) # Wall
        error("Wall boundary conditions not implemented yet")
    elseif (BC == 4) # User defined
        error("User defined boundary conditions not implemented yet")
    else
        error("Invalid boundary condition: $(BC)")
    end
end