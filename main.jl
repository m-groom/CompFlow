# FV solver for the solution of the 1D Euler equations

# Load functions
include("riemann_solver.jl")
include("equation_of_state.jl")
include("system.jl")
include("reconstruction.jl")
include("plotting.jl")
include("initial_condition.jl")
include("timestepping.jl")
include("grid.jl")

# Define the domain. TODO: read this from an input file
xL = 0;
xR = 1;
imax = 100; # number of control volumes
Nmax = 10000; # Maximum number of time steps
CFL = 0.9; # Courant-Friedrichs-Lewy number
γ = 1.4; # Ratio of specific heats
# Build the grid
x = gridGen(xL, xR, imax);

# Initial condition
test = 2; # test case to use (0-6)
Q, tend = initialCondition(x, test, γ);

# Plot the initial condition. TODO: make this a general plotting function
fig1 = plotIC(x, Q, γ);

# Compute the approximate solution using the MUSCL-Hancock scheme
# TODO: save the solution at specified intervals
t = 0;
for n = 1:Nmax
    # Compute the time step
    Δt = getTimeStep(Q, x, γ, CFL);
    if (t + Δt > tend)
        Δt = tend - t;
    end
    # Stop criterion
    if (t >= tend)
        break
    end
    # Reconstruct the extrapolated values at the cell boundary
    QR, QL = reconstruct(Q, γ);
    # Evolve the extrapolated values at the cell boundary
    QR, QL = evolve(QR, QL, x, Δt, γ);
    # Compute the fluxes. TODO: turn this into a function called update
    Qnew = zeros(3, imax);
    for i = 1:imax
        if (i==1)
            # Dirichlet boundary condition. TODO: generalise this
            Qnew[:,i] = Q[:,i];
        elseif (i==imax)
            # Dirichlet boundary condition. TODO: generalise this
            Qnew[:,i] = Q[:,i];
        else
            Δx = x[i+1] - x[i];
            Fp = HLLC(QR[:,i], QL[:,i+1], γ)
            Fm = HLLC(QR[:,i-1], QL[:,i], γ)
            Qnew[:,i] = Q[:,i] - Δt/Δx * (Fp - Fm);
        end
    end 
    # Update the time and the solution
    global t = t + Δt;
    global Q = Qnew;
end

# Plot the solution. TODO: save the solution as well
fig2 = plotSolution(x, Q, γ, t, test);