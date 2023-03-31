# Functions for plotting the solver output

# Load modules
using PyPlot
# Load functions
include("system.jl")
include("riemann_solver.jl")
include("initial_condition.jl")

# Plot the initial condition
function plotIC(x, Q, γ)
    # Calculate cell centres
    xc = collect(x[1:end-1] .+ 0.5 * (x[2:end] .- x[1:end-1]));
    # Convert conserved to primitive variables
    W = zeros(3, length(xc));
    for i = 1:length(xc)
        W[:, i] = consToPrim(Q[:, i], γ);
    end
    # Plot the initial condition
    fig = figure(1)
    subplot(3, 1, 1)
    plot(xc, W[1,:], "bo", fillstyle="none")
    title("Initial condition")
    xlabel("x")
    ylabel("ρ")
    subplot(3, 1, 2)
    plot(xc, W[2,:], "bo", fillstyle="none")
    xlabel("x")
    ylabel("u")
    subplot(3, 1, 3)
    plot(xc, W[3,:], "bo", fillstyle="none")
    xlabel("x")
    ylabel("p")
    plt.savefig("IC.png")

    return fig
end

# Plot the solution
# TODO: make this into a general plotting function for each time interval (plus IC)
# TODO: only plot the exact solution if test > 0
function plotSolution(x, Q, γ, t, test = 0)
    # Calculate cell centres
    xc = collect(x[1:end-1] .+ 0.5 * (x[2:end] .- x[1:end-1]));
    # Convert conserved to primitive variables
    W = zeros(3, length(xc));
    for i = 1:length(xc)
        W[:, i] = consToPrim(Q[:, i], γ);
    end
    # Calculate the exact solution
    xe, We = exactSolution(x, t, γ, test);
    # Plot the solution
    fig = figure(2)
    subplot(3, 1, 1)
    plot(xc, W[1,:], "bo", fillstyle="none", zorder = 5)
    plot(xe, We[1,:], "k-", zorder = 1)
    title("Solution")
    xlabel("x")
    ylabel("ρ")
    subplot(3, 1, 2)
    plot(xc, W[2,:], "bo", fillstyle="none", zorder = 5)
    plot(xe, We[2,:], "k-", zorder = 1)
    xlabel("x")
    ylabel("u")
    subplot(3, 1, 3)
    plot(xc, W[3,:], "bo", fillstyle="none", zorder = 5)
    plot(xe, We[3,:], "k-", zorder = 1)
    xlabel("x")
    ylabel("p")
    plt.savefig("Solution.png")

    return fig
end

# Calculate the exact solution
function exactSolution(x, t, γ, test, Ne = 10000)
    # Compute grid for exact solution
    xL = x[1]; xR = x[end];
    xe = collect(range(xL, xR, length=Ne));
    # Initialise arrays
    Qe = zeros(3,Ne);
    We = zeros(3,Ne);
    # Get left and right data
    Q1, Q2, x0, _ = riemannProblem(test, γ);
    # Calculate solution from exact Riemann solver
    for i = 1:Ne
        ξ = (xe[i] - x0) / t;
        Qe[:,i] = exactRiemannSolver(Q1, Q2, ξ, γ);
        # Convert to primitive variables
        We[:,i] = consToPrim(Qe[:,i], γ);
    end

    return xe, We

end