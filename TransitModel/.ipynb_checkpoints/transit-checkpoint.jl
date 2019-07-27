include("limbdark/src/transit_poly_struct.jl")

# produces an array of impact parameters for a transit assuming constant velocity 
function b_array(t, t0, b0, rp, duration)
    b = zeros(Real, length(t))
    chordlength = 2*sqrt((1+rp)^2-b0^2)
    v = chordlength/duration
    
    for i in 1:length(t)
        if abs(t[i]-t0) > chordlength/(2*v)
            b[i] = 1+rp
        else
            x = v*(t[i]-t0)
            b[i] = sqrt(b0^2 + x^2)
        end
    end
    return b
end

# transit light curve at times t (center=t0, maximum impact parameter = b0, radius=rp, duration=d, limb darkening parameters = u)
function transit(t, t0, b0, rp, d, u)
    b = b_array(t, t0, b0, rp, d)
    trans = transit_init(rp, b0, u, true)
    flux = zeros(Real, length(b))
    gradflux = zeros(Real, length(b), 2 + length(u))
    for i in 1:length(b)
        trans.b = b[i]
        flux[i] = transit_poly!(trans)
        gradflux[i,:] = trans.dfdrbu
    end
    return flux, gradflux
end
