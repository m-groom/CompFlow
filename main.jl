# CompFlow: A finite-volume solver for the solution of the 1D Euler equations

# Load functions
include("src/riemann_solver.jl")
include("src/equation_of_state.jl")
include("src/system.jl")
include("src/reconstruction.jl")
include("src/output.jl")
include("src/initial_condition.jl")
include("src/timestepping.jl")
include("src/grid.jl")
include("src/logging.jl")

# Define the domain
report("Defining the computational domain...")
x = makeGrid("grid.par");
# Define the solver settings
report("Defining the solver settings...")
Nmax, Nout, CFL = solverSettings("solver.par");
# Print the solver settings
# Define the fluid properties
report("Defining the fluid properties...")
γ = fluidProperties("thermo.par");
# Set the initial time
t = 0.0;

# Initial condition
test = 1; # test case to use (0-6)
report("Computing the initial condition...")
report("Test case: $(test)")
Q, tend = initialCondition(x, test, γ);

# Plot the initial condition
report("Plotting the solution at time t = $(round(t, digits=6))")
fig1 = plotSolution(x, Q, γ, t, test);

# Save the initial condition
report("Saving the solution at time t = $(rpad(string(round(t, digits=6)), 8, "0"))")
writeSolution(x, Q, "solution_$(rpad(string(round(t, digits=6)), 8, "0")).vtr")

# Compute the approximate solution
startTime = report("Computing the approximate solution...", 1)
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
    report("Current time step: $(n)   Current time: $(rpad(string(round(t + Δt, digits=6)), 8, "0"))")
    # Reconstruct the extrapolated values at the cell boundary
    QR, QL = reconstruct(Q, γ);
    # Evolve the extrapolated values at the cell boundary
    QR, QL = evolve(QR, QL, x, Δt, γ);
    # Compute the solution at the next timestep
    Qnew = update(QR, QL, Q, x, Δt, γ);
    # Update the time and the solution
    global t = t + Δt;
    global Q = Qnew;
    # Write the solution at every Nout time steps
    if (n%Nout == 0)
        report("Saving the solution at time t = $(rpad(string(round(t, digits=6)), 8, "0"))")
        writeSolution(x, Q, "solution_$(rpad(string(round(t, digits=6)), 8, "0")).vtr")
    end
end
endTime = report("Simulation completed...", 1)
report("Elapsed time: $(endTime - startTime)")

# Plot the final solution
report("Plotting the solution at time t = $(round(t, digits=6))")
fig2 = plotSolution(x, Q, γ, t, test);

# Save the final solution
report("Saving the solution at time t = $(rpad(string(round(t, digits=6)), 8, "0"))")
writeSolution(x, Q, "solution_$(rpad(string(round(t, digits=6)), 8, "0")).vtr")