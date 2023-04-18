# Functions for plotting and saving the solver output

# Load modules
import PyPlot as plt
using WriteVTK
# Load functions
include("system.jl")
include("riemann_solver.jl")
include("logging.jl")
include("../initial_condition.jl")

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

# Read fluid properties from a file
function fluidProperties(filename)
    # Read the fluid properties from a file
    report("Reading the fluid properties from file $(filename)");
    file = open(filename, "r");
    γ = parse(Float64,split(readline(file), "#")[1]);
    close(file)
    report("Ratio of specific heats: $(γ)")
    
    return γ
end

# Save the solution to a VTK file
function writeSolution(x, Q, filename)
    # Define y and z vectors
    y = [0.0; 1.0];
    z = [0.0; 1.0];
    # Initialise arrays
    ρ = zeros(length(x)-1, length(y)-1, length(z)-1);
    ρu = zeros(length(x)-1, length(y)-1, length(z)-1);
    ρE = zeros(length(x)-1, length(y)-1, length(z)-1);
    # Fill arrays
    for i = 1:length(x)-1
        ρ[i, 1, 1] = Q[1, i];
        ρu[i, 1, 1] = Q[2, i];
        ρE[i, 1, 1] = Q[3, i];
    end
    # Write to .vtr file
    vtk_grid(filename, x, y, z) do vtk
        vtk["Density"] = ρ;
        vtk["Momentum"] = ρu;
        vtk["Energy"] = ρE;
    end

end

# Plot the solution
function plotSolution(x, Q, γ, t, test)
    # Calculate cell centres
    xc = collect(x[1:end-1] .+ 0.5 * (x[2:end] .- x[1:end-1]));
    # Convert conserved to primitive variables
    W = zeros(3, length(xc));
    for i = 1:length(xc)
        W[:, i] = consToPrim(Q[:, i], γ);
    end
    # Calculate the exact solution
    xe, We = exactSolution(x, t, γ, test);
    xec = xe[1:end-1] + 0.5 * (xe[2:end]-xe[1:end-1]);
    # Plot the solution
    fig = plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(xc, W[1,:], "bo", fillstyle="none", zorder = 5)
    plt.plot(xec, We[1,:], "k-", zorder = 1);
    plt.title("Solution: t = $(t)")
    plt.grid(true)
    plt.xlabel("x")
    plt.xticks(round.(linspace(xc[1], xc[end], 6), digits=1))
    plt.ylabel("ρ")
    plt.yticks(round.(linspace(minimum(W[1,:]), maximum(W[1,:]), 5), digits=2))
    plt.subplot(3, 1, 2)
    plt.plot(xc, W[2,:], "bo", fillstyle="none", zorder = 5)
    plt.plot(xec, We[2,:], "k-", zorder = 1);
    plt.grid(true)
    plt.xlabel("x")
    plt.xticks(round.(linspace(xc[1], xc[end], 6), digits=1))
    plt.ylabel("u")
    plt.yticks(round.(linspace(minimum(W[2,:]), maximum(W[2,:]), 5), digits=2))
    plt.subplot(3, 1, 3)
    plt.plot(xc, W[3,:], "bo", fillstyle="none", zorder = 5)
    plt.plot(xec, We[3,:], "k-", zorder = 1);
    plt.grid(true)
    plt.xlabel("x")
    plt.xticks(round.(linspace(xc[1], xc[end], 6), digits=1))
    plt.ylabel("p")
    plt.yticks(round.(linspace(minimum(W[3,:]), maximum(W[3,:]), 5), digits=2))
    # Save the plot
    plt.savefig("solution_$(t).png")

    return fig
end

# Calculate the exact solution
function exactSolution(x, t, γ, test, Ne = 10001)
    # Compute grid for exact solution
    xL = x[1]; xR = x[end];
    xe = linspace(xL, xR, Ne);
    # Initialise arrays
    Qe = zeros(3,Ne-1);
    We = zeros(3,Ne-1);
    # Get left and right data
    Q1, Q2, x0, _ = riemannProblem(test, γ);
    if (t > 0)
        # Calculate solution from exact Riemann solver
        for i = 1:Ne-1
            # Calculate the cell centre
            xc = xe[i] + 0.5 * (xe[i+1] - xe[i]);
            ξ = (xc - x0) / t;
            Qe[:,i] = exactRiemannSolver(Q1, Q2, ξ, γ);
            # Convert to primitive variables
            We[:,i] = consToPrim(Qe[:,i], γ);
        end
    else
        # Return the initial condition
        Qe, _ = initialCondition(xe, test, γ);
        # Convert to primitive variables
        for i = 1:Ne-1
            We[:,i] = consToPrim(Qe[:,i], γ);
        end
    end

    return xe, We

end

# Linspace function
function linspace(start, stop, n)
    return collect(range(start, stop, length=n))
end