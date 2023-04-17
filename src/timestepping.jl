# Functions for updating the solution at each timestep

# Load functions
include("system.jl")
include("riemann_solver.jl")
include("logging.jl")

# Read the solver settings from a file
function solverSettings(filename)
    # Read the solver settings from a file
    report("Reading the solver settings from file $(filename)");
    file = open(filename, "r");
    Nmax = parse(Int,split(readline(file), "#")[1]);
    Nout = parse(Int,split(readline(file), "#")[1]);
    CFL = parse(Float64,split(readline(file), "#")[1]);
    close(file)
    report("Maximum number of time steps: $(Nmax)")
    report("Number of time steps between outputs: $(Nout)")
    report("CFL number: $(CFL)")
    
    return Nmax, Nout, CFL
end

# Function for estimating the timestep size
function getTimeStep(Q, x, γ, CFL)
    # Number of cells
    imax = size(Q, 2); 
    # Estimate the maximum wave speed
    Smax = 0.0;
    Δx = x[end] - x[1];
    for i = 1:imax
        ai = speedOfSound(Q[:, i], γ);
        ui = Q[2, i] / Q[1, i];
        Smax = max(Smax, abs(ui) + ai);
        Δx = min(Δx, x[i+1] - x[i]);
    end
    
    return CFL * Δx / Smax

end

# Function for updating the solution
# TODO: generalise boundary conditions
function update(QR, QL, Q, x, Δt, γ)
    imax = length(x) - 1; # number of cells
    Qnew = zeros(3, imax);
    for i = 1:imax
        Δx = x[i+1] - x[i]; # grid spacing
        if (i==1)
            Qnew[:,i] = Q[:,i]; # Dirichlet boundary condition
        elseif (i==imax)
            Qnew[:,i] = Q[:,i]; # Dirichlet boundary condition
        else
            Fp = riemannSolver(QR[:,i], QL[:,i+1], γ) # Flux at x[i+1/2]
            Fm = riemannSolver(QR[:,i-1], QL[:,i], γ) # Flux at x[i-1/2]
            Qnew[:,i] = Q[:,i] - Δt/Δx * (Fp - Fm);
        end
    end 

    return Qnew

end