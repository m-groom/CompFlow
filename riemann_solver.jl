# Functions for computing the solution of the Riemann problem
include("equation_of_state.jl")
include("newton.jl")

# Compute the exact solution of the Riemann problem (from p152 of Toro)
function ExactRiemannSolver(QL, QR, ξ, γ)
    # Calculate primitive variables
    DL = QL[1];        # left density
    DR = QR[1];        # right density
    UL = QL[2]/QL[1];  # left velocity
    UR = QR[2]/QR[1];  # right velocity
    PL = eos(QL, γ);   # left pressure
    PR = eos(QR, γ);   # right pressure

    # Compute γ related constants
    G1 = (γ - 1.0) / (2.0 * γ);
    G2 = (γ + 1.0) / (2.0 * γ);
    G3 = 2.0 * γ / (γ - 1.0);
    G4 = 2.0 / (γ - 1.0);
    G5 = 2.0 / (γ + 1.0);
    G6 = (γ - 1.0) / (γ + 1.0);
    G7 = (γ - 1.0) / 2.0;
    G8 = γ - 1.0;
    # Package these into a single vector
    G = [G1; G2; G3; G4; G5; G6; G7; G8];

    # Compute the speeds of sound
    CL = sos(QL, γ);
    CR = sos(QR, γ);

    # Check pressure positivity condition
    if (G4 * (CL + CR) <= UR - UL)
        error("Pressure positivity condition violated")
    end

    # Compute exact solution for pressure and velocity in the star region
    WL = [DL; UL; PL]; # left primitive variables
    WR = [DR; UR; PR]; # right primitive variables
    PS, US = newton(WL, WR, CL, CR, G);




    
    # Return vector of conserved variables
    # return [h; h*u; h*psi]
end

