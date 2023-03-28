# Functions for computing the solution of the Riemann problem
include("equation_of_state.jl")
include("newton.jl")
include("system.jl")

# Compute the exact solution of the Riemann problem (from p152 of Toro)
# TODO: fix error, currently not giving the correct solution
function exactRiemannSolver(QL, QR, ξ, γ)
    # Calculate primitive variables
    DL, UL, PL = consToPrim(QL, γ);
    DR, UR, PR = consToPrim(QR, γ);

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
    CL = speedOfSound(QL, γ);
    CR = speedOfSound(QR, γ);

    # Check pressure positivity condition
    if (G4 * (CL + CR) <= UR - UL)
        error("Pressure positivity condition violated")
    end

    # Compute exact solution for pressure and velocity in the star region
    WL = [DL; UL; PL]; # left primitive variables
    WR = [DR; UR; PR]; # right primitive variables
    PM, UM = newton(WL, WR, CL, CR, G);

    # Sample the solution
    if (ξ <= UM) # Left of the contact discontinuity
        if (PM <= PL) # Left rarefaction
            SHL = UL - CL;
            if (ξ <= SHL) # Left data state
                ρ = WL[1];
                u = WL[2];
                p = WL[3];
            else
                CML = CL*(PM/PL)^(G[1]);
                STL = UM - CML;
                if (ξ > STL) # Left star region
                    ρ = DL*(PM/PL)^(1.0/γ);
                    u = UM;
                    p = PM;
                else # Inside rarefaction fan
                    C = G[5]*(CL + G[7]*(UL - ξ));
                    ρ = DL*(C/CL)^(G[4]);
                    u = G[5]*(CL + G[7]*UL + ξ);
                    p = PL*(C/CL)^(G[3]);
                end
            end
        else # Left shock
            PML = PM/PL;
            SL = UL - CL*sqrt(G[2]*PML + G[1]);
            if (ξ <= SL) # Left data state
                ρ = WL[1];
                u = WL[2];
                p = WL[3];
            else # Left star region
                ρ = DL*(PML + G[6])/(PML*G[6] + 1.0);
                u = UM;
                p = PM;
            end
        end
    else # right of contact discontinuity
        if (PM > PR) # Right shock
            PMR = PM/PR;
            SR = UR + CR*sqrt(G[2]*PMR + G[1]);
            if (ξ >= SR) # Right data state
                ρ = WR[1];
                u = WR[2];
                p = WR[3];
            else # Right star region
                ρ = DR*(PMR + G[6])/(PMR*G[6] + 1.0);
                u = UM;
                p = PM;
            end
        else # Right rarefaction
            SHR = UR + CR;
            if (ξ >= SHR) # Right data state
                ρ = WR[1];
                u = WR[2];
                p = WR[3];
            else
                CMR = CR*(PM/PR)^(G[1]);
                STR = UM + CMR;
                if (ξ <= STR) # Right star region
                    ρ = DR*(PM/PR)^(1.0/γ);
                    u = UM;
                    p = PM;
                else # Inside rarefaction fan
                    C = G[5]*(CR - G[7]*(UR - ξ));
                    ρ = DR*(C/CR)^(G[4]);
                    u = G[5]*(-CR + G[7]*UR + ξ);
                    p = PR*(C/CR)^(G[3]);
                end
            end
        end
    end

    # Convert primitive to conserved variables
    Q = primToCons([ρ; u; p], γ);
    
    # Return vector of conserved variables
    return Q

end

# HLLC solver
function HLLC(QL, QR, γ)
    # Calculate primitive variables
    WL = consToPrim(QL, γ);
    WR = consToPrim(QR, γ);

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
    CL = speedOfSound(QL, γ);
    CR = speedOfSound(QR, γ);

    # Estimate the wave speeds
    SL, SM, SR = waveSpeeds(WL, WR, CL, CR, G);

    # # Rusanov solver
    # Smax = max(abs(SL), abs(SR));
    # flux = 0.5*(Fa(QL, γ) + Fa(QR, γ)) - 0.5*Smax*(QR - QL);

    # Compute the flux
    if (SL >= 0) # Supersonic flow to the right
        flux = Fa(QL, γ);
    elseif (SR <= 0) # Supersonic flow to the left
        flux = Fa(QR, γ);
    else
        # # HLL solver
        # flux = (SR*Fa(QL, γ) - SL*Fa(QR, γ) + SR*SL*(QR - QL))/(SR - SL);
        if (SM >= 0) # Subsonic flow to the right
            DL = QL[1]; EL = QL[3]/QL[1]; UL = WL[2]; PL = WL[3];
            QM = DL * (SL - UL) / (SL - SM) * [1.0; SM; EL + (SM-UL)*(SM + PL/(DL*(SL-UL)))];
            flux = Fa(QL, γ) + SL*(QM - QL);
        else # Subsonic flow to the left
            DR = QR[1]; ER = QR[3]/QR[1]; UR = WR[2]; PR = WR[3];
            QM = DR * (SR - UR) / (SR - SM) * [1.0; SM; ER + (SM-UR)*(SM + PR/(DR*(SR-UR)))];
            flux = Fa(QR, γ) + SR*(QM - QR);
        end
    end

    return flux

end

# Function for estimating wave speeds for the HLLC solver
function waveSpeeds(WL, WR, CL, CR, G)
    # Extract the primitive variables
    DL, UL, PL = WL;
    DR, UR, PR = WR;
    # Estimate the pressure and velocity in the star region
    PM, _ = ANRS(WL, WR, CL, CR, G);
    
    # Compute the left and right wave speeds
    if (PM <= PL)
        SL = UL - CL;
    else 
        SL = UL - CL*sqrt(1.0 + G[2]*(PM/PL - 1.0));
    end

    if (PM <= PR)
        SR = UR + CR;
    else
        SR = UR + CR*sqrt(1.0 + G[2]*(PM/PR - 1.0));
    end

    SM = (PR - PL + DL*UL*(SL-UL) - DR*UR*(SR-UR)) / (DL*(SL-UL) - DR*(SR-UR));

    return SL, SM, SR

end