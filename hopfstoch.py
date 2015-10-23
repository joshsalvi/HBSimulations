
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def hopfstoch(mu, fosc, NoiseSTD, tvec, Fextmax, fr):

    # Local Variables: xNoiseSTD, j, xsto, Dtfac, fosc, Ftime, ysto, mu, yNoiseSTD, ydet, pulsestim, yzero, fr, sinusoidalstim, pulseend, xzero, N, Xdet, Dt, Xsto, Fext, NoiseSTD, plotyn, pulsestart, xdW, xdet, tvec, Fextmax, ydW
    # Function calls: real, subplot, randn, figure, sawtooth, imag, title, plot, sqrt, RandStream, length, zeros, linspace, hopfstoch, pi
    #%
    #% This function simulates the normal form of the supercritical Hopf
    #% bifurcation, given by two planar equations:
    #%
    #% x_dot = mu*x - omega*y - x*(x^2 + y^2)
    #% y_dot = omega*x + mu*y - y*(x^2 + y^2)
    #%
    #% where mu is the control parameter. For mu>0, the system will oscillate at
    #% an amplitude that grows with sqrt(mu). Alternatively, one may express the
    #% above equations in polar coordinates, making the amplitude relationship
    #% with respect to mu more apparent:
    #%
    #% rho_dot = rho*(mu + i*omega - rho^2)
    #%
    #% Here we simulate both the deterministic and stochastic cases for the
    #% supercritical Hopf bifurcation 
    #%
    #% [Xdet Xsto] = hopfstoch(mu,fosc,NoiseSTD,tvec,Fextmax,fr)
    #%
    #% Xdet : deterministic result
    #% Xsto : stochastic result
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
    xzero = 1.
    yzero = -1.
    #% Add external forcing if desired
    sinusoidalstim = 0.
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
    xdW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    ydW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    xNoiseSTD = NoiseSTD
    yNoiseSTD = NoiseSTD
    xdet = np.zeros(1., N)
    xdet[0] = xzero
    xsto = np.zeros(1., N)
    xsto[0] = xzero
    ydet = np.zeros(1., N)
    ydet[0] = yzero
    ysto = np.zeros(1., N)
    ysto[0] = yzero
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
        
    Xdet = np.zeros(2., (length(tvec)-1.))
    Xsto = np.zeros(2., (length(tvec)-1.))
    #%Return vectors at times specified by Time.
    Xdet[0,:] = xdet[0:N:Dtfac]
    Xdet[1,:] = ydet[0:N:Dtfac]
    Xsto[0,:] = xsto[0:N:Dtfac]
    Xsto[1,:] = ysto[0:N:Dtfac]
    Fext = Fext[0:N:Dtfac]
    #% Make a plot of the data?
    plotyn = 0.
    if plotyn == 1.:
        plt.figure
        plt.subplot(1., 2., 1.)
        plt.hold(on)
        plt.plot(tvec[1:], Xsto[0,:], 'r')
        plt.plot(tvec[1:], Xdet[0,:], 'k')
        plt.title('Black=deterministic; Red=stochastic; real part only')
        plt.subplot(1., 2., 2.)
        plt.hold(on)
        plt.plot(Xsto[0,:], Xsto[1,:], 'r')
        plt.plot(Xdet[0,:], Xdet[1,:], 'k')
        plt.title('Black=deterministic; Red=stochastic')
    
    
    return [Xdet, Xsto, Fext]