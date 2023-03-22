# Functions for equation of state calculations

# Takes the vector of conserved variables as input and returns the pressure
function eos(Q, gamma = 5/3, eosType = "ideal")
    rho = Q[1]; u = Q[2]/Q[1]; rhoE = Q[3];
    if (eosType == "ideal")
        p = (gamma - 1.0)*(rhoE - 0.5*rho*u^2);
    else
        error("EOS type currently not supported")
    return p
end

# Takes the vector of conserved variables as input and returns the speed of sound
function sos(Q, gamma = 5/3, eosType = "ideal")
    rho = Q[1]; u = Q[2]/Q[1]; rhoE = Q[3];
    if (eosType == "ideal")
        p = (gamma - 1.0)*(rhoE - 0.5*rho*u^2);
        c = sqrt(gamma * p / rho);
    else
        error("EOS type currently not supported")
    return c
end