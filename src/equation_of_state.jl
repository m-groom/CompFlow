# Functions for equation of state calculations

# Take the vector of conserved variables as input and returns the pressure
function pressure(Q, γ, eosType="ideal")
    ρ = Q[1]
    u = Q[2] / Q[1]
    ρE = Q[3]
    if (eosType == "ideal")
        p = (γ - 1.0) * (ρE - 0.5 * ρ * u^2)
    else
        error("EOS type $eosType currently not supported")
    end
    return p
end

# Take the vector of conserved variables as input and returns the speed of sound
function speedOfSound(Q, γ, eosType="ideal")
    ρ = Q[1]
    u = Q[2] / Q[1]
    ρE = Q[3]
    if (eosType == "ideal")
        p = (γ - 1.0) * (ρE - 0.5 * ρ * u^2)
        c = sqrt(γ * p / ρ)
    else
        error("EOS type $eosType currently not supported")
    end
    return c
end

# Take the vector of primitive variables as input and returns the internal energy
function internalEnergy(W, γ, eosType="ideal")
    ρ = W[1]
    p = W[3]
    if (eosType == "ideal")
        e = p / ((γ - 1) * ρ)
    else
        error("EOS type $eosType currently not supported")
    end
    return e
end