include("transit_orb.jl")

function transit(t, x)
    nsub = 5
    n = length(t)
    flux::Array{Number, 1} = zeros(n)
    for i=1:n
        ti = t[i]
        if i > 1
            dt = ti-t[i-1]
        else
            dt = 1e-10
        end
        flux[i] = transit_orb(ti, x, dt, nsub)
    end
    return flux
end