
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def FHNstoch(I, NoiseSTD, tvec, Fextmax, fr):

    # Local Variables: j, Dtfac, wsto, wdet, Ftime, Vsto, pulsestim, wNoiseSTD, sinusoidalstim, Vdet, Dt, vsto, vzero, fr, vdW, I, pulseend, N, vdet, vNoiseSTD, wdW, wzero, Fext, NoiseSTD, plotyn, pulsestart, tvec, Fextmax
    # Function calls: real, subplot, randn, figure, sawtooth, title, plot, sqrt, RandStream, FHNstoch, length, zeros, linspace, pi
    #%
    #% This function simulates the FitzHugh-Nagumo model:
    #%
    #% v_dot = v - v^3/3 - w + I
    #% w_dot = 0.08*(v + 0.7 - 0.8*w)
    #%
    #% where I is an injected current and the control parameter. 
    #% V is the membrane potential.
    #% W is a recovery variable.
    #%
    #% Here we simulate both the deterministic and stochastic cases for the
    #% supercritical Hopf bifurcation 
    #%
    #% [Vdet, Vsto, Fext] = FHNstoch(I,NoiseSTD,tvec,Fextmax,fr)
    #%
    #% Vdet : deterministic result (V(1,:)=v, V(2,:)=w)
    #% Vsto : stochastic result
    #%
    #%  
    #% tvec : time vector
    #% mu : control parameter
    #% xNoiseSTD : standard deviation of stochastic noise in x
    #% yNoiseSTD : standard deviation of stochastic noise in y
    #% fosc : frequency of oscillation on the unstable side of the bifurcation
    #% Fextmax : amplitude of periodic forcing
    #% fr : frequency of driving
    #%
    #% By modifying the code, you can also add a step function or external
    #% forcing. 
    #%
    #% jsalvi@rockefeller.edu
    #%
    #% Initial condition
    vzero = 0.1
    wzero = -0.1
    #% Add external forcing if desired
    sinusoidalstim = 1.
    pulsestim = 0.
    #% pulse or sinusoid?
    #%Fextmax = 1e8;        % amplitude of sinusoidal stim OR pulse
    #%fr = 1.01;             % frequency of stimulation
    pulsestart = 1.
    #% start of pulse
    pulseend = 2.
    #% end of pulse
    #% Decrease time step size by factor of Dtfac to ensure convergence
    Dtfac = 10.**2.
    Dt = matdiv(tvec[1]-tvec[0], Dtfac)
    N = matdiv(tvec[int(0)-1], Dt)
    #%Set the default random number stream
    RandStream.setGlobalStream[int(RandStream[int('mt19937ar')-1,int('seed')-1,0])-1]
    vdW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    wdW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    vNoiseSTD = NoiseSTD
    wNoiseSTD = NoiseSTD
    vdet = np.zeros(1., N)
    vdet[0] = vzero
    vsto = np.zeros(1., N)
    vsto[0] = vzero
    wdet = np.zeros(1., N)
    wdet[0] = wzero
    wsto = np.zeros(1., N)
    wsto[0] = wzero
    #% External forcing
    if sinusoidalstim == 1.:
        Ftime = np.linspace(tvec[0], tvec[int(0)-1], N)
        #%Fext = Fextmax*sin(2*pi*fr*Ftime);
        Fext = np.dot(Fextmax, sawtooth(np.dot(np.dot(2.*np.pi, fr), Ftime)))
    elif pulsestim == 1.:
        Ftime = np.linspace(tvec[0], tvec[int(0)-1], N)
        Fext = np.dot((Ftime<pulseend)-(Ftime<pulsestart), Fextmax)
        
    else:
        Fext = np.zeros(1., N)
        
    
    #% Euler-Murayama Method with Ito Integration
    for j in np.arange(2., (N)+1):
        #%Deterministic integral
        
    Vdet = np.zeros(2., (length(tvec)-1.))
    Vsto = np.zeros(2., (length(tvec)-1.))
    #%Return vectors at times specified by Time.
    Vdet[0,:] = vdet[0:N:Dtfac]
    Vdet[1,:] = wdet[0:N:Dtfac]
    Vsto[0,:] = vsto[0:N:Dtfac]
    Vsto[1,:] = wsto[0:N:Dtfac]
    Fext = Fext[0:N:Dtfac]
    #% Make a plot of the data?
    plotyn = 1.
    if plotyn == 1.:
        plt.figure
        plt.subplot(1., 2., 1.)
        plt.hold(on)
        plt.plot(tvec[1:], Vsto[0,:], 'r')
        plt.plot(tvec[1:], Vdet[0,:], 'k')
        plt.title('Black=deterministic; Red=stochastic; real part only')
        plt.subplot(1., 2., 2.)
        plt.hold(on)
        plt.plot(Vsto[0,:], Vsto[1,:], 'r')
        plt.plot(Vdet[0,:], Vdet[1,:], 'k')
        plt.title('Black=deterministic; Red=stochastic')
    
    
    return [Vdet, Vsto, Fext]