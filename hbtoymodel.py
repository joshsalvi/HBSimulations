
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def hbtoymodel(Fc, k, noiselevel, Fextmax, fr, tvec):

    # Local Variables: xNoiseSTD, tau, xsto, fdet, Dtfac, noiselevel, Ftvec, noiseampl, fzero, fNoiseSTD, t2n, pulsestim, noisevec, Fc, Fext2, Dt, fr, sinusoisalstim, xzero, fsto, N, Xdet, NVec, fdW, x0, plotyn, ndW, a, b, Xsto, Fext, k, j, xdW, xdet, tvec, Fextmax, whitenoise, t1n
    # Function calls: plot, cos, randn, figure, sqrt, RandStream, length, hbtoymodel, zeros, linspace, ylabel, xlabel, pi, round, findnearest
    #%
    #% This function simulates the hair-buyndle model from PNAS 2012.
    #%
    #% [Xdet, Xsto, Fext] = hbtoymodel(Fc,k,noiselevel,Fextmax,fr,tvec)
    #%
    #% Xdet : deterministic result
    #% Xsto : stochastic result
    #%
    #%  
    #% tvec : tvec vector
    #% Fc,k : control parameters
    #% noiselevels : standard deviation of stochastic noise in x and y
    #% fr : frequency of oscillation on the unstable side of the bifurcation
    #% Fextmax : amplitude in force of sinusoidal stimulation.
    #%
    #% Note that stiffnesses are scaled by a factor of 100 in the manuscript.
    #%
    #% By modifying the code, you can also add a step function or other external
    #% forcing. 
    #%
    #% jsalvi@rockefeller.edu
    #%
    #%Stochasic HB model integration
    #%EM Euler-Maruyama method
    #%Ito integral
    a = 3.5
    #%b > 1 has unbounded solutions
    b = 0.5
    tau = 10.
    xzero = 0.01
    fzero = 0.
    x0 = 0.
    #%Decrease tvec step size by factor of Dtfac to ensure convergence
    Dtfac = 10.**2.
    Dt = matdiv(tvec[1]-tvec[0], Dtfac)
    N = np.round(matdiv(tvec[int(0)-1], Dt))
    #%Set the default random number stream
    RandStream.setGlobalStream[int(RandStream[int('mt19937ar')-1,int('seed')-1,0])-1]
    xdW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    fdW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    ndW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% white noise driving
    xdet = np.zeros(1., N)
    fdet = np.zeros(1., N)
    xsto = np.zeros(1., N)
    fsto = np.zeros(1., N)
    xdet[0] = xzero
    xsto[0] = xzero
    fdet[0] = fzero
    fsto[0] = fzero
    #%Not using FD theorem
    xNoiseSTD = noiselevel
    fNoiseSTD = noiselevel
    #% equal noise levels
    sinusoisalstim = 0.
    Ftvec = np.linspace(tvec[0], tvec[int(0)-1], N)
    Fext = np.zeros(1., N)
    if sinusoisalstim == 1.:
        Fext = np.dot(Fextmax, np.cos(np.dot(np.dot(2.*np.pi, fr), Ftvec)))
        #%Fext = abs(Fextmax*sawtooth(2*pi*fr*Ftvec)) - mean(abs(Fextmax*sawtooth(2*pi*fr*Ftvec)));
        #%Fext = Fextmax.*Ftvec;
    
    
    pulsestim = 0.
    #% for a pulse stimulus, fr should be [t1 t2], in which:
    #% t1: starting time of pulse
    #% t2: ending time of pulse
    if pulsestim == 1.:
        t1n = findnearest(Ftvec, fr[0])
        t2n = findnearest(Ftvec, fr[1])
        t1n = t1n[0]
        t2n = t2n[0]
        Fext[int(t1n)-1:t2n] = Fextmax
    
    
    whitenoise = 0.
    if whitenoise == 1.:
        Fext = Fextmax*plt.randn(1., length(Ftvec))
    
    
    for j in np.arange(2., (N)+1):
        noiseampl = 0.2
        noisevec[int(j)-1] = ndW[int(j)-1]
        #%Deterministic integral
        xdet[int(j)-1] = xdet[int((j-1.))-1]+np.dot(Dt, np.dot(-k, xdet[int((j-1.))-1]-x0)+np.dot(a, xdet[int((j-1.))-1]-fdet[int((j-1.))-1])-(xdet[int((j-1.))-1]-fdet[int((j-1.))-1])**3.+Fc+Fext[int(j)-1])+noiseampl*ndW[int(j)-1]
        fdet[int(j)-1] = fdet[int((j-1.))-1]+matdiv(np.dot(Dt, np.dot(b, xdet[int((j-1.))-1])-fdet[int((j-1.))-1]), tau)
        #%Stochastic integral
        xsto[int(j)-1] = xsto[int((j-1.))-1]+np.dot(Dt, np.dot(-k, xsto[int((j-1.))-1]-x0)+np.dot(a, xsto[int((j-1.))-1]-fsto[int((j-1.))-1])-(xsto[int((j-1.))-1]-fsto[int((j-1.))-1])**3.+Fc+Fext[int(j)-1])+np.dot(xNoiseSTD, xdW[int(j)-1])+noiseampl*ndW[int(j)-1]
        fsto[int(j)-1] = fsto[int((j-1.))-1]+matdiv(np.dot(Dt, np.dot(b, xsto[int((j-1.))-1])-fsto[int((j-1.))-1]), tau)+matdiv(np.dot(fNoiseSTD, fdW[int(j)-1]), tau)
        
    Xdet = np.zeros(2., (length(tvec)-1.))
    Xsto = np.zeros(2., (length(tvec)-1.))
    #%Return vectors at tvecs specified by tvec.
    Xdet[0,:] = xdet[0:N:Dtfac]
    Xdet[1,:] = fdet[0:N:Dtfac]
    Xsto[0,:] = xsto[0:N:Dtfac]
    Xsto[1,:] = fsto[0:N:Dtfac]
    Fext2[0,:] = Fext[0:N:Dtfac]
    NVec[0,:] = noisevec[0:N:Dtfac]
    plotyn = 0.
    if plotyn == 1.:
        #% close all
    plt.figure
    plt.plot(Ftvec[0:], xsto, 'r')
    plt.hold(on)
    plt.plot(Ftvec[0:], xdet, 'k')
    plt.xlabel('tvec', 'FontSize', 24.)
    plt.ylabel('x', 'FontSize', 24., 'Rotation', 0., 'HorizontalAlignment', 'right')
    #%}
    #%{
    plt.figure
    plt.plot(Ftvec[0:], fsto, 'g')
    plt.hold(on)
    plt.plot(Ftvec[0:], fdet, 'k')
    plt.xlabel('tvec', 'FontSize', 24.)
    plt.ylabel('f', 'FontSize', 24., 'Rotation', 0., 'HorizontalAlignment', 'right')
    #%}
    #% figure
    #% plot(xsto,fsto,'k');xlabel('xsto');ylabel('fsto');
    
    return [Xdet, Xsto, Fext2, NVec]