# CompFlow: A finite-volume solver for the solution of the 1D Euler equations

# Load modules
import PyPlot as plt
import Dates
using WriteVTK
# Load functions
include("src/equation_of_state.jl")
include("src/file_io.jl")
include("src/grid.jl")
include("src/logging.jl")
include("src/reconstruction.jl")
include("src/riemann_solver.jl")
include("src/system.jl")
include("src/timestepping.jl")
include("src/user_defined.jl")
include("initial_condition.jl")

# Define the domain
report("Defining the computational domain...")
x = makeGrid("grid.par");
# Define the solver settings
report("Defining the solver settings...")
Nmax, Nout, CFL = solverSettings("solver.par");
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
t = computeSolution!(Q, x, γ, CFL, Nmax, Nout, t, tend);
endTime = report("Simulation completed...", 1)
report("Elapsed time: $(endTime - startTime)")

# Plot the final solution
report("Plotting the solution at time t = $(round(t, digits=6))")
fig2 = plotSolution(x, Q, γ, t, test);

# Save the final solution
report("Saving the solution at time t = $(rpad(string(round(t, digits=6)), 8, "0"))")
writeSolution(x, Q, "solution_$(rpad(string(round(t, digits=6)), 8, "0")).vtr")