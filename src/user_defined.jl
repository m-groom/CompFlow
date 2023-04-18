# User defined functions

# Newton iteration to compute p* and u* for the exact Riemann solver
# Note: currently only valid for fixed γ and ideal gas EOS
function newton(WL, WR, CL, CR, G)
    # Extract velocities
    UL = WL[2]; UR = WR[2]; 

    # Newton iteration to compute p* and u*
    PS, _ = ANRS(WL, WR, CL, CR, G);   # Initial guess for p*
    ϵ = 1e-12;                      # tolerance
    maxIter = 100;                  # maximum number of iterations

    # Compute pressure in the star region
    for i=1:maxIter
        g, dg = f(PS, WL, WR, CL, CR, G);  # function to solve and its derivative
        res = abs(g); 
        if(res < ϵ)        
            break # we have found the solution so we exit
        end
        Δp = -g / dg; # Newton step   
        # line search globalization 
        δ = 1; # scale factor 0<δ<=1 
        for ii=1:20   
            gk, _ = f(PS + δ * Δp, WL, WR, CL, CR, G);       
            if(abs(gk) >= res)
                # if the residual of the next iteration increases then the Newton step is reduced by 2
                δ = 0.5 * δ; 
            else 
                # residual does not increase => exit inner loop
                break
            end
        end
        # Update the solution. For δ=1, globalization method reduces to standard Newton. 
        PS = PS + δ * Δp; 
    end

    # Compute velocity in the star region
    fR, _ = fK(PS, WR, CR, G); fL, _ = fK(PS, WL, CL, G);
    US = 0.5 * ((UL + UR) + (fR - fL));

    return PS, US

end

# Function gk used in the root finding algorithm
function f(PS, WL, WR, CL, CR, G)
    # Extract velocities
    UL = WL[2]; UR = WR[2]; 
    # Evaluate fR, fL and their derivatives
    fR, dfR = fK(PS, WR, CR, G);
    fL, dfL = fK(PS, WL, CL, G);
    # Compute g and its derivative
    g = fR + fL + UR - UL;
    dg = dfR + dfL;

    return g, dg

end

# Function fK containing the Rankine-Hugoniot relation, Riemann invariants and Lax entropy condition
function fK(P, WK, CK, G)
    # Extract primitive variables
    DK = WK[1]; PK = WK[3];
    AK = G[5] / DK;
    BK = G[6] / DK;

    # Lax entropy condition
    if (P > PK) 
        # shock (Rankine-Hugoniot conditions)
        ϕ = (P - PK) * sqrt(AK / (P + BK));
        dϕ = (1.0 - 0.5*(P - PK)/(BK + P)) * sqrt(AK / (P + BK));
    else
        # rarefaction (Riemann invariants) 
        ϕ = G[4] * CK * ((P / PK)^(G[1]) - 1.0);
        dϕ = (1.0 / (DK*CK)) * (P / PK)^(-G[2]);
    end

    return ϕ, dϕ

end