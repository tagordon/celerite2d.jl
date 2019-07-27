import ftplib
import numpy as np
import os
import batman
import astropy.constants as ac
import julia
import corner
import starry
import sys

def get_sim(simfile):
    jl = julia.Julia()
    jl.eval("using Mamba")
    jl.eval("include(\"utils.jl\")")
    jl.eval("include(\"../src/celerite.jl\")")
    if not jl.eval("@isdefined CeleriteDist"):
        jl.eval("eval(DefineCeleriteDist())")
        #jl.eval(jl.DefineCeleriteDist())
    sim = jl.eval("sim = read(\"{0}\", ModelChains)".format(simfile))
    sim = jl.eval("sim.value")
    labels = jl.eval("sim.names")
    return sim, labels

def make_cornerplot(simfile, chains, burnin=0, parameters=slice(None, None), truths=None):
    sim, labels = get_sim(simfile)
    return corner.corner(sim[burnin:, parameters, chains], labels=labels, truths=truths)

def get_names(simfile, parameters=slice(None, None)):
    sim, labels = get_sim(simfile)
    return labels[parameters]

def get_samples(simfile):
    sim, _ = get_sim(simfile)
    return sim

def get_statistics(simfile):
    jl = julia.Julia()
    jl.eval("using Mamba")
    jl.eval("include(\"utils.jl\")")
    jl.eval("include(\"../src/celerite.jl\")")
    if not jl.eval("@isdefined CeleriteDist"):
        jl.eval("eval(DefineCeleriteDist())")
    sim = jl.eval("sim = read(\"{0}\", ModelChains)".format(simfile))
    return jl.eval("summarystats(sim)")
    
def get_autocor(simfile, lags):
    lags = iterable_to_string(lags)
    jl = julia.Julia()
    jl.eval("using Mamba")
    jl.eval("include(\"utils.jl\")")
    jl.eval("include(\"../src/celerite.jl\")")
    if not jl.eval("@isdefined CeleriteDist"):
        jl.eval("eval(DefineCeleriteDist())")
    sim = jl.eval("sim = read(\"{0}\", ModelChains)".format(simfile))
    return jl.eval("autocor(sim, lags={0})".format(lags))

def iterable_to_string(arr):
    st = "["
    for a in arr[:-1]:
        st += "{0},".format(a)
    st += "{0}]".format(arr[-1])
    return st

def estimate_autocor_time(simfile, chain, C=5):
    M = 1
    estimate = 1 + 2*get_autocor(simfile, [M]).value[:, :, chain][:, 0]
    while all(M < C*estimate):
        M += 1
        estimate += 2*get_autocor(simfile, [M]).value[:, :, chain][:, 0]
    return estimate

def plot_samples(simfile, nsamp, xlims=None, ylims=None):
    data = np.loadtxt(simfile)
    t = data[:, 0]
    pl.figure(figsize=(12, 8))
    pl.plot(t, transit(t, 10, 0.0, 0.07, 100, 0.06, [0.5, 0.5]), 'k', linewidth=5)
    for i in range(nsamp):
        i = np.random.randint(np.shape(sim_poly)[0])
        d, b0, t0, rp = sim_poly[i, :, 1][9], sim_poly[i, :, 1][10], sim_poly[i, :, 1][11], np.exp(sim_poly[i, :, 1][12])
        u1, u2 = sim_poly[i, :, 1][13], sim_poly[i, :, 1][8]
        u = [u1, u2]
        p = 100
        pl.plot(t, transit(t, t0, b0, d, p, rp, u), 'b', alpha=0.1)
        d, b0, t0, rp = sim_mono[i, :, 1][2], sim_mono[i, :, 1][6], sim_mono[i, :, 1][3], np.exp(sim_mono[i, :, 1][4])
        u1, u2 = sim_mono[i, :, 1][9], sim_mono[i, :, 1][8]
        u = [u1, u2]
        pl.plot(t, transit(t, t0, b0, d, p, rp, u), 'r', alpha=0.1)
    if xlims is not None:
        pl.xlim(*xlims)
    if ylims is not None:
        pl.ylim(*ylims)
                
def make_data(t, log_S0, log_w0, log_wn, rp, d, tin, c, t0=None, log_Q=1/np.sqrt(2), shape="square", seed=None):
    if t0 is None:
        t0 = np.median(t)
    jl = julia.Julia()
    jl.eval("using Random")
    jl.eval("include(\"utils.jl\")")
    jl.eval("include(\"../src/celerite.jl\")")
    if not jl.eval("@isdefined CeleriteDist"):
        jl.eval("eval(DefineCeleriteDist())")
    t_string = iterable_to_string(t)
    c_string = iterable_to_string(c)
    jl.eval("cd = CeleriteDist({0}, {1}, {2}, {3}, log_wn=ones(length({5}))*{4}, c={5})".format(t_string, log_S0, log_Q, log_w0, log_wn, c_string))
    jl.eval("x = zeros({0})".format(len(t)*len(c)))
    if seed is not None:
        jl.eval("Random.seed!({0})".format(seed))
    jl.eval("_rand!(cd, x)")
    r, g, b = jl.eval("x[1:3:end], x[2:3:end], x[3:3:end]")
    
    if shape is "square":
        trans = np.ones(len(t))
        intransit = (t>(t0-d/2))&(t<(t0+d/2))
        trans[intransit] = 1-(rp**2)
    if shape is "trapezoid":
        trans = np.ones(len(t))
        t1 = (t0-d/2-tin)
        t2 = (t0-d/2)
        t3 = (t0+d/2)
        t4 = (t0+d/2+tin)
        intransit = (t<=t3)&(t>=t2)
        ingress = (t<t2)&(t>t1)
        egress = (t<t4)&(t>t3)
        trans[intransit] = 1-(rp**2)
        trans[ingress] = 1 -(t[ingress]-t1)*(rp**2)/tin
        trans[egress] = (1 - (rp**2))+(t[egress]-t3)*(rp**2)/tin
    x = jl.eval("x")
    return [x[i::len(c)] for i in range(len(c))]

# generate an entire set of simulations that differ only by the steepness of 
# the wavelength dependence, having the same transit and variability aprameters 
# otherwise. All will sum to the same monochromatic light curve, and thus 
# will differ in white noise between simulations (the white noise is 
# identicle between bands). wn refers to the white noise of the reference simulation 
# with c=1
def make_dataset(t, log_S0, log_w0, depth_over_amp, depth_over_sigma, duration_times_w0, egress_over_duration, c, nbins, bs, t0=None, log_Q=1/np.sqrt(2), shape="square", seed=None):
    sims = np.zeros((np.shape(bs)[1]-1, nbins, len(t)))
    fluxes = make_data(t, log_S0=log_S0, log_w0=log_w0, log_wn=-10, rp=0, d=0, tin=0, c=1+c*np.arange(nbins), log_Q=log_Q, shape="trapezoid")
    depth = np.std(sum(fluxes)/nbins)*depth_over_amp
    wn = np.sqrt(nbins)*depth/depth_over_sigma
    d = duration_times_w0/np.exp(log_w0)
    tin = egress_over_duration*d
    sig = np.random.randn(nbins, len(t))*wn
    rp = np.sqrt(depth)
    if t0 is None:
        t0 = np.median(t)
    if shape is "square":
        trans = np.ones(len(t))
        intransit = (t>(t0-d/2))&(t<(t0+d/2))
        trans[intransit] = 1-(rp**2)
    if shape is "trapezoid":
        trans = np.ones(len(t))
        t1 = (t0-d/2-tin)
        t2 = (t0-d/2)
        t3 = (t0+d/2)
        t4 = (t0+d/2+tin)
        intransit = (t<=t3)&(t>=t2)
        ingress = (t<t2)&(t>t1)
        egress = (t<t4)&(t>t3)
        trans[intransit] = 1-(rp**2)
        trans[ingress] = 1 -(t[ingress]-t1)*(rp**2)/tin
        trans[egress] = (1 - (rp**2))+(t[egress]-t3)*(rp**2)/tin
    
    for i in range(np.shape(bs)[0]):
        bins = bs[i]
        x = np.zeros((len(bins)-1, len(t)))
        for j in range(1, len(bins)):
            scales = 1+c*np.arange(nbins)
            norm = bins[j]-bins[j-1]
            x[j-1,:] = sum(fluxes[bins[j-1]:bins[j]]) + norm*trans + np.sqrt(sum(sig[bins[j-1]:bins[j]]**2))
        sims[:, i, :] = x
    params = [wn, rp, d, tin, c, log_S0, log_Q, log_w0, t0]
    return sims, params, trans