# FV solver for the solution of the 1D Euler equations

using PyPlot
include("riemann_solver.jl")
include("equation_of_state.jl")

γ = 1.4;

test = 6; # test case to use (0-6)

# Define the problem
if (test == 0) # RP0 from Toro
    ρ1 = 1.0;
    ρ2 = 0.125;
    u1 = 0.0;
    u2 = 0.0;
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
Q1 = conserved(W1, γ);
Q2 = conserved(W2, γ);

# Define the domain
t = 0;
xL = 0;
xR = 1;
imax = 100; # number of control volumes
dx = (xR - xL) / imax;
x = collect(xL+dx/2:dx:xR-dx/2);
Nmax = 10000; # Maximum number of time steps
CFL = 0.5; # Courant-Friedrichs-Lewy number

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
figure(1)
subplot(3, 1, 1)
plot(x, W[1,:], "bo", fillstyle="none")
title("Initial condition")
xlabel("x")
ylabel("ρ")
subplot(3, 1, 2)
plot(x, W[2,:], "bo", fillstyle="none")
xlabel("x")
ylabel("u")
subplot(3, 1, 3)
plot(x, W[3,:], "bo", fillstyle="none")
xlabel("x")
ylabel("p")
plt.savefig("IC.png")

# Compute the approxiamte solution using MUSCL-Hancock
# for n = 1:Nmax
#     # Compute the time step
#     amax = 0.0;
#     for i = 1:imax
#         eig = lambda(Q[:, i], g);
#         amax = max(amax, maximum(abs.(eig))); # Maximum eigenvalue
#     end
#     dt = CFL * dx / amax;
#     if (t + dt > tend)
#         dt = tend - t;
#     end
#     # Stop criterion
#     if (t >= tend)
#         break
#     end
#     # Piecewise linear reconstruction in space
#     slope = zeros(3,imax);
#     QR = zeros(3,imax);
#     QL = zeros(3,imax);
#     for i = 1:imax
#         if (i==1)
#             slope[:,i] = zeros(3,1); 
#         elseif (i==imax)
#             slope[:,i] = zeros(3,1); 
#         else
#             slope[:,i] = minmod.(Q[:,i] - Q[:,i-1], Q[:,i+1] - Q[:,i]);
#             # slope[:,i] = superbee.(Q[:,i] - Q[:,i-1], Q[:,i+1] - Q[:,i]);
#             # slope[:,i] = eno.(Q[:,i] - Q[:,i-1], Q[:,i+1] - Q[:,i]);
#         end
#         # Compute the extrapolated values at the cell boundary
#         QR[:,i] = Q[:,i] + 0.5 * slope[:,i];
#         QL[:,i] = Q[:,i] - 0.5 * slope[:,i];
#         # Compute the time derivative (Cauchy-Kovalevskaya): Q_t = -F_x
#         Qt = - (F(QR[:,i], g) - F(QL[:,i], g)) / dx;
#         # Update the extrapolated values
#         QR[:,i] = QR[:,i] + 0.5 * dt * Qt;
#         QL[:,i] = QL[:,i] + 0.5 * dt * Qt;
#     end 
#     # Compute the fluxes
#     Qnew = zeros(3, imax);
#     for i = 1:imax
#         if (i==1)
#             # Dirichlet boundary condition
#             Qnew[:,i] = Q1;
#             # # Reflective boundary condition
#             # QGod = ExactRiemannSolver(QR[:,i], QL[:,i+1], 0, g);
#             # Fp = F(QGod, g);
#             # QBC = [QL[1,i]; -QL[2,i]; -QL[3,i]];
#             # QGod = ExactRiemannSolver(QBC, QL[:,i], 0, g);
#             # Fm = F(QGod, g);
#             # Qnew[:,i] = Q[:,i] - dt/dx * (Fp - Fm);
#         elseif (i==imax)
#             # Dirichlet boundary condition
#             Qnew[:,i] = Q2;
#             # # Reflective boundary condition
#             # QBC = [QR[1,i]; -QR[2,i]; -QR[3,i]];
#             # QGod = ExactRiemannSolver(QR[:,i], QBC, 0, g);
#             # Fp = F(QGod, g);
#             # QGod = ExactRiemannSolver(QR[:,i-1], QL[:,i], 0, g);
#             # Fm = F(QGod, g);
#             # Qnew[:,i] = Q[:,i] - dt/dx * (Fp - Fm);
#         else
#             QGod = ExactRiemannSolver(QR[:,i], QL[:,i+1], 0, g);
#             Fp = F(QGod, g);
#             QGod = ExactRiemannSolver(QR[:,i-1], QL[:,i], 0, g);
#             Fm = F(QGod, g);
#             Qnew[:,i] = Q[:,i] - dt/dx * (Fp - Fm);
#         end
#     end 
#     # Update the time and the solution
#     global t = t + dt;
#     global Q = Qnew;
# end

t = tend;

# Compute the exact solution
xe = collect(range(xL, xR, length=10000));
Qe = zeros(3,10000);
We = zeros(3,10000);
for i = 1:10000
    ξ = (xe[i] - x0) / t;
    Qe[:,i] = exactRiemannSolver(Q1, Q2, ξ, γ);
    We[:,i] = primitive(Qe[:,i], γ);
end

# Plot the solution as a 3x1 subplot
figure(2)
subplot(3, 1, 1)
# plot(x, Q[1,:], "bo", fillstyle="none", zorder = 5)
plot(xe, We[1,:], "k-", zorder = 1)
title("Solution")
xlabel("x")
ylabel("ρ")
subplot(3, 1, 2)
# plot(x, Q[2,:], "bo", fillstyle="none", zorder = 5)
plot(xe, We[2,:], "k-", zorder = 1)
xlabel("x")
ylabel("u")
subplot(3, 1, 3)
# plot(x, Q[3,:], "bo", fillstyle="none", zorder = 5)
plot(xe, We[3,:], "k-", zorder = 1)
xlabel("x")
ylabel("p")
plt.savefig("Solution.png")

