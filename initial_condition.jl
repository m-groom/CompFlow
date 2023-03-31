# Functions for calculating the initial condition

# Load functions
include("system.jl")

# Calculate the initial condition
function initialCondition(x, test, γ)
    # Initialise array
    imax = length(x)-1;
    Q = zeros(3, imax);
    # Get left and right data
    Q1, Q2, x0, tend = riemannProblem(test, γ);
    # Assign to Q
    for i = 1:imax
        if x[i] < x0
            Q[:, i] = Q1;
        else
            Q[:, i] = Q2;
        end
    end

    return Q, tend

end

# Define a Riemann problem
function riemannProblem(test, γ)
    # Test case to use
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

    return Q1, Q2, x0, tend

end