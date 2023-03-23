# Functions for equation of state calculations

# Take the vector of conserved variables as input and returns the pressure
function eos(Q, γ, eosType = "ideal")
    ρ = Q[1]; u = Q[2]/Q[1]; ρE = Q[3];
    if (eosType == "ideal")
        p = (γ - 1.0)*(ρE - 0.5*ρ*u^2);
    else
        error("EOS type $eosType currently not supported")
    return p
end

# Take the vector of conserved variables as input and returns the speed of sound
function sos(Q, γ = 5/3, eosType = "ideal")
    ρ = Q[1]; u = Q[2]/Q[1]; ρE = Q[3];
    if (eosType == "ideal")
        p = (γ - 1.0)*(ρE - 0.5*ρ*u^2);
        c = sqrt(γ * p / ρ);
    else
        error("EOS type $eosType currently not supported")
    return c
end