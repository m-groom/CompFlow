# Compute the exact solution of the Riemann problem
include("equation_of_state.jl")
function ExactRiemannSolver(QL, QR, xi, gamma)
    # Calculate primitive variables
    DL = QL[1];            # left density
    DR = QR[1];            # right density
    UL = QL[2]/QL[1];      # left velocity
    UR = QR[2]/QR[1];      # right velocity
    PL = eos(QL, gamma);   # left pressure
    PR = eos(QR, gamma);   # right pressure

    # Compute gamma related constants
    G1 = (gamma - 1.0) / (2.0 * gamma);
    G2 = (gamma + 1.0) / (2.0 * gamma);
    G3 = 2.0 * gamma / (gamma - 1.0);
    G4 = 2.0 / (gamma - 1.0);
    G5 = 2.0 / (gamma + 1.0);
    G6 = (gamma - 1.0) / (gamma + 1.0);
    G7 = (gamma - 1.0) / 2.0;
    G8 = gamma - 1.0;

    # Compute the speeds of sound
    CL = sos(QL, gamma);
    CR = sos(QR, gamma);

    # Check pressure positivity condition
    if (G4 * (CL + CR) <= UR - UL)
        error("Pressure positivity condition violated")
    end



    
    # Return vector of conserved variables
    # return [h; h*u; h*psi]
end