using Distributed 
addprocs(length(Sys.cpu_info()))
using LinearAlgebra
@everywhere basepath = "/gscratch/astro/tagordon/celmcmc/"
@everywhere srcpath = basepath*"src/"
@everywhere datapath = basepath*"data/full/"
@everywhere outpath = basepath*"data/full/tmp"
@everywhere inputpath = basepath*"input/full/"
@everywhere include(srcpath*"GPdist/celerite.jl")
@everywhere include(srcpath*"TransitModel/transit.jl")

## Define a new multivariate Distribution type for Mamba.
## The definition must be placed within an unevaluated quote block.
@everywhere extensions = quote

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
        
        function CeleriteDist(t, log_s, log_q, log_w; log_wn=-12, c=[1], mu=zeros(Base.length(c), Base.length(t)))
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
		wn_vec[i:m:end] .= exp.(log_wn[i])
	end
        try
            celerite.compute!(d.gp, d.t, wn_vec)
            return celerite.log_likelihood(d.gp, x.-d.mu)
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
            x[i:m:end] .+= d.mu[i,:]
        end
    end
    
    mean(d::CeleriteDist) = convert(Array{Real, 2}, d.mu)
    
end

@everywhere include(srcpath*"Mamba.jl/src/Mamba.jl")
@everywhere using Main.Mamba
@everywhere eval(extensions)

@everywhere using DelimitedFiles
@everywhere codefile = inputpath*"codes.in"
@everywhere paramfile = inputpath*"params_poly.in"
@everywhere runcodes = readdlm(codefile, String)
@everywhere params = readdlm(paramfile)
@everywhere paramdict = Dict(runcodes[i] => params[i,:] for i in 1:length(runcodes))

# build model 

@everywhere TransitModel = Model(

    # defines the probability of y given mu and sig as a multivariate normal likelihood function 
    y = Stochastic(1,
        (t, c, log_s, log_w, log_q, log_wn, trans) -> CeleriteDist(t, log_s, log_q, log_w, log_wn=log_wn, mu=trans, c=[c[1], c[2]]),
        false
    ),
    
    trans = Logical(1,
        (t0, b0, rp, d, u, t) -> transit(t, t0, b0, rp, d, u)[1],
        false,
    ),
    
    # transit parameters
    t0 = Stochastic(
        () -> Uniform(-1, 1)
    ),
    b0 = Stochastic(
        () -> Uniform(0, 1)
    ),
   
    log_rp = Stochastic(
        () -> Uniform(-20, 0),
        false
    ),
    rp = Logical(
        (log_rp) -> exp(log_rp)
    ),
    d = Stochastic(
        () -> Uniform(0, 10)
    ),    
    q = Stochastic(1,
        () -> UnivariateDistribution[Uniform(0, 1), Uniform(0, 1)],
        false
    ),
    u = Logical(1,
        (q) -> [2*sqrt(q[1])*q[2], sqrt(q[1])*(1-2*q[2])]
    ),

    #noise parameters
    log_s = Stochastic(
        () -> Uniform(-30, 0)
    ),
    log_q = Stochastic(
        () -> Uniform(-10, 10)
    ),
    log_w = Stochastic(
        () -> Uniform(-10, 10)
    ),
    log_wn = Stochastic(1,
        () -> UnivariateDistribution[Uniform(-40, 10), Uniform(-40, 10)]
    ),
    c = Stochastic(1,
        () -> UnivariateDistribution[Uniform(0, 10), Uniform(0, 10)]
    )
)

# define the sampling scheme: No U-Turn Sampling for the three parameters 
@everywhere scheme = [AMWG([:t0, :b0, :log_rp, :d, :q, :log_s, :log_q, :log_w, :log_wn, :c], 0.1*ones(13))]

# define the initial values for the walkers 
@everywhere epsilon = 1e-8
@everywhere nchains = 4
if nprocs() > length(runcodes) nchains = Int(floor(nprocs()/length(runcodes))) end

@everywhere function f(code)
    p = paramdict[code]
    log_wn1, log_wn2, c1, c2, b0, rp, d, q1, q2, log_s, log_q, log_w, t0 = p
    transfile = datapath*"transit_poly_"*code
    fluxdata = readdlm(transfile)
                
    data = Dict{Symbol, Any}(
        :t => fluxdata[:,1],
        :y => vec(permutedims(fluxdata[:, 2:end]))
    )
    inits = [
        Dict{Symbol, Any}(
            :y => data[:y],
                        
            #transit parameters
            :t0 => t0 + epsilon*randn(),
            :b0 => b0 + epsilon*randn(),
            :log_rp => log(rp) + epsilon*randn(),
            :d => d + epsilon*randn(),
            :q => [q1, q2] + epsilon*randn(2),
            
            #noise parameters
            :log_s => log_s + epsilon*randn(),
            :log_q => log_q + epsilon*randn(),
            :log_w => log_w + epsilon*randn(),
            :log_wn => [log_wn1, log_wn2] + epsilon*randn(2),
            :c => [c1, c2] + epsilon*randn(2)
        )
        for i in 1:nchains
    ]
                
    setsamplers!(TransitModel, scheme)
    simfile = outpath*"simpoly"*code*".jls"
    sim = mcmc(TransitModel, data, inits, 10, burnin=0, chains=nchains)
    write(simfile, sim)
end
pmap(f, runcodes)
