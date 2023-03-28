# FV solver for the solution of the 1D Euler equations
# TODO: write an analogous semi-discrete scheme

using PyPlot
include("riemann_solver.jl")
include("equation_of_state.jl")
include("system.jl")
include("reconstruction.jl")

γ = 1.4;

# Define the domain
# TODO: read this from an input file
t = 0;
xL = 0;
xR = 1;
imax = 100; # number of control volumes
Nmax = 10000; # Maximum number of time steps
CFL = 0.5; # Courant-Friedrichs-Lewy number
# TODO: turn this into a function
dx = (xR - xL) / imax;
x = collect(xL:dx:xR);

# Define the problem
test = 6; # test case to use (0-6)

if (test == 0) # RP0 from Toro
    ρ1 = 1.0;
    ρ2 = 0.125;
    u1 = 0.5; # 0.0;
    u2 = 0.5; # 0.0;
    p1 = 1.0;
    p2 = 1.0;
    x0 = 0.5;
    tend = 0.2;
elseif (test == 1) # RP1 from Toro
    ρ1 = 1.0;
    ρ2 = 0.125;
    u1 = 0.75;
    u2 = 0.0;
    p1 = 1.0;
    p2 = 0.1;
    x0 = 0.3;
    tend = 0.2;
elseif (test == 2) # RP2 from Toro
    ρ1 = 1.0;
    ρ2 = 1.0;
    u1 = -2.0;
    u2 = 2.0;
    p1 = 0.4;
    p2 = 0.4;
    x0 = 0.5;
    tend = 0.15;
elseif (test == 3) # RP3 from Toro
    ρ1 = 1.0;
    ρ2 = 1.0;
    u1 = 0.0;
    u2 = 0.0;
    p1 = 1000.0;
    p2 = 0.01;
    x0 = 0.5;
    tend = 0.012;
elseif (test == 4) # RP4 from Toro
    ρ1 = 5.99924;
    ρ2 = 5.99242;
    u1 = 19.5975;
    u2 = -6.19633;
    p1 = 460.894;
    p2 = 46.0950;
    x0 = 0.4;
    tend = 0.035;
elseif (test == 5) # RP5 from Toro
    ρ1 = 1.0;
    ρ2 = 1.0;
    u1 = -19.59745;
    u2 = -19.59745;
    p1 = 1000.0;
    p2 = 0.01;
    x0 = 0.8;
    tend = 0.012;
elseif (test == 6) # RP6 from Toro
    ρ1 = 1.0;
    ρ2 = 1.0;
    u1 = 2.0;
    u2 = -2.0;
    p1 = 0.1;
    p2 = 0.1;
    x0 = 0.5;
    tend = 0.8;
else
    error("Invalid test case")
end

# Primitive variables
W1 = [ρ1; u1; p1];
W2 = [ρ2; u2; p2];

# Conserved variables
Q1 = primToCons(W1, γ);
Q2 = primToCons(W2, γ);

# Initial condition
Q = zeros(3, imax);
W = zeros(3, imax);
for i = 1:imax
    if x[i] < x0
        Q[:, i] = Q1;
        W[:, i] = W1;
    else
        Q[:, i] = Q2;
        W[:, i] = W2;
    end
end

# Plot the initial condition as a 3x1 subplot
# TODO: turn this into a function
xc = collect(xL+dx/2:dx:xR-dx/2);
figure(1)
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

# Compute the approximate solution using the MUSCL-Hancock scheme
for n = 1:Nmax
    # Compute the time step
    # TODO: turn this into a function
    Smax = 0.0;
    for i = 1:imax
        ai = speedOfSound(Q[:, i], γ);
        ui = Q[2, i] / Q[1, i];
        Smax = max(Smax, abs(ui) + ai);
    end
    dt = CFL * dx / Smax;
    if (t + dt > tend)
        dt = tend - t;
    end
    # Stop criterion
    if (t >= tend)
        break
    end
    # Reconstruct the extrapolated values at the cell boundary
    QR, QL = reconstruct(Q, γ);
    # Evolve the extrapolated values at the cell boundary
    QR, QL = evolve(QR, QL, x, dt, γ);
    # Compute the fluxes
    # TODO: turn this into a function called update
    Qnew = zeros(3, imax);
    for i = 1:imax
        if (i==1)
            # Dirichlet boundary condition
            Qnew[:,i] = Q1;
        elseif (i==imax)
            # Dirichlet boundary condition
            Qnew[:,i] = Q2;
        else
            Fp = HLL(QR[:,i], QL[:,i+1], γ)
            Fm = HLL(QR[:,i-1], QL[:,i], γ)
            Qnew[:,i] = Q[:,i] - dt/dx * (Fp - Fm);
        end
    end 
    # Update the time and the solution
    global t = t + dt;
    global Q = Qnew;
end

# Convert to primitive variables
for i = 1:imax
    W[:, i] = consToPrim(Q[:, i], γ);
end

# Compute the exact solution
# TODO: turn this into a function
Ne = 10000;
xe = collect(range(xL+dx/2, xR-dx/2, length=Ne));
Qe = zeros(3,Ne);
We = zeros(3,Ne);
for i = 1:Ne
    ξ = (xe[i] - x0) / t;
    Qe[:,i] = exactRiemannSolver(Q1, Q2, ξ, γ);
    We[:,i] = consToPrim(Qe[:,i], γ);
end

# Plot the solution as a 3x1 subplot
# TODO: turn this into a function
figure(2)
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