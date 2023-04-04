# Functions for plotting the solver output

# Load modules
using PyPlot
# Load functions
include("system.jl")
include("riemann_solver.jl")
include("initial_condition.jl")

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
    fig = figure()
    subplot(3, 1, 1)
    plot(xc, W[1,:], "bo", fillstyle="none", zorder = 5)
    plot(xec, We[1,:], "k-", zorder = 1);
    title("Solution: t = $(t)")
    grid(true)
    xlabel("x")
    plt.xticks(round.(linspace(xc[1], xc[end], 6), digits=1))
    ylabel("ρ")
    plt.yticks(round.(linspace(minimum(W[1,:]), maximum(W[1,:]), 5), digits=2))
    subplot(3, 1, 2)
    plot(xc, W[2,:], "bo", fillstyle="none", zorder = 5)
    plot(xec, We[2,:], "k-", zorder = 1);
    grid(true)
    xlabel("x")
    plt.xticks(round.(linspace(xc[1], xc[end], 6), digits=1))
    ylabel("u")
    plt.yticks(round.(linspace(minimum(W[2,:]), maximum(W[2,:]), 5), digits=2))
    subplot(3, 1, 3)
    plot(xc, W[3,:], "bo", fillstyle="none", zorder = 5)
    plot(xec, We[3,:], "k-", zorder = 1);
    grid(true)
    xlabel("x")
    plt.xticks(round.(linspace(xc[1], xc[end], 6), digits=1))
    ylabel("p")
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