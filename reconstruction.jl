# Functions for reconstructing slopes
# TODO: implement characteristic-based reconstruction
# TODO: use different reconstruction schemes for different characteristic variables

# Minmod slope limiter
function minmod(a, b)
    if (b >= 0)
        return max(0, min(a, b))
    else
        return min(0, max(a, b))
    end
    # return max(0, min(a, b))
end

# Superbee slope limiter
function superbee(a, b)
    if (b >= 0)
        return max(0, min(2.0*a, b), min(a, 2.0*b))
    else
        return min(0, max(2.0*a, b), max(a, 2.0*b))
    end
    # return max(0, min(2.0*a, b), min(a, 2.0*b))
end