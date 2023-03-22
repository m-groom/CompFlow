# Compute the exact solution of the Riemann problem
include("equation_of_state.jl")
function ExactRiemannSolver(QL, QR, xi, γ = 5/3)
    # Calculate primitive variables
    DL = QL[1];            # left density
    DR = QR[1];            # right density
    UL = QL[2]/QL[1];      # left velocity
    UR = QR[2]/QR[1];      # right velocity
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

    # Compute the speeds of sound
    CL = sos(QL, γ);
    CR = sos(QR, γ);

    # Check pressure positivity condition
    if (G4 * (CL + CR) <= UR - UL)
        error("Pressure positivity condition violated")
    end



    
    # Return vector of conserved variables
    # return [h; h*u; h*psi]
end