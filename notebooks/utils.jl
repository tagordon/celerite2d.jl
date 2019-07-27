using FITSIO
using Statistics
using DelimitedFiles

# read in SOHO data with long timescale variations filtered out 
function read_soho(path)
    red = FITS("$path/SPM_red_intensity_series.fits")
    blue = FITS("$path/SPM_blue_intensity_series.fits")
    green = FITS("$path/SPM_green_intensity_series.fits")

    tstart = read_header(blue[1])["JDSTART"]
    tend = read_header(blue[1])["JDEND"]

    red = read(red[1])
    blue = read(blue[1])
    green = read(green[1])

    time = (collect(LinRange(tstart, tend, length(blue))) .- tstart)
    return time, red, green, blue
end

function write_transit(t, r, g, b, path)
    open(path, "w") do io
        writedlm(io, [t r g b])
    end
end

function DefineCeleriteDist()
    
    extensions = quote

        ## Load needed packages and import methods to be extended
        using Distributions
        import Distributions: length, insupport, _logpdf
    
        ## Type declaration
        mutable struct CeleriteDist <: ContinuousMultivariateDistribution
            t::Array{Real, 1}
            c::Array{Real, 1}
            log_s::Real
            log_q::Real
            log_w::Real
            log_wn::Array{Real, 1}
            mu::Array{Real, 1}
            gp::celerite.Celerite
        
            function CeleriteDist(t, log_s, log_q, log_w; log_wn=[-12], c=[1], mu=zeros(Base.length(t)))
                gp = celerite.Celerite(celerite.SHOTerm(log_s, log_q, log_w), broadcast(*, c, c'))
                new(t, c, log_s, log_q, log_w, log_wn, mu, gp)
            end
        end

        ## The following method functions must be implemented

        ## Dimension of the distribution
        length(d::CeleriteDist) = length(d.t)*length(d.c)

        ## Logical indicating whether x is in the support
        function insupport(d::CeleriteDist, x::AbstractVector{T}) where {T<:Real}
            length(d) == length(x) && all(isfinite.(x))
        end
    
        function _logpdf(d::CeleriteDist, x::AbstractVector{T}) where {T<:Real}
            m = Base.length(d.c)
            n = length(d.t)
            wn_vec = zeros(n*m)
            for i in 1:m
                wn_vec[i:m:end] .= exp(d.log_wn[i])
            end
            try
                celerite.compute!(d.gp, d.t, wn_vec)
                z = zeros(length(d.t)*m)
                for i in 1:m
                    z[i:m:end] .= x[i:m:end] .- d.mu
                end
                return celerite.log_likelihood(d.gp, z)
            catch
                return -Inf 
            end
        end
    
        function _rand!(d::CeleriteDist, x::AbstractVector{T}) where {T<:Real}
            m = Base.length(d.c)
            n = length(d.t)
            wn_vec = zeros(n*m)
            for i in 1:m
                wn_vec[i:m:end] .= exp(d.log_wn[i])
            end
            celerite.compute!(d.gp, d.t, wn_vec)
            u = randn(Base.length(d.t)*Base.length(d.c))
            x[:] = celerite.simulate_gp(d.gp, u)
            for i in 1:m
                x[i:m:end] .+= d.mu
            end
        end
    
        mean(d::CeleriteDist) = convert(Vector, d.mu)
    
    end
    return extensions
end