
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def cuspstoch(b1, b2, xNoiseSTD, tvec):

    # Local Variables: xNoiseSTD, j, xsto, Dtfac, Ftime, b1, b2, pulsestim, Dt, fr, sinusoidalstim, pulseend, xzero, N, Xdet, Xsto, Fext, dens, plotyn, pulsestart, xdet, bw, tvec, Fextmax, dW
    # Function calls: subplot, plot, cos, randn, kde2d, figure, title, imag, cuspstoch, sqrt, RandStream, length, zeros, linspace, real, pi, hilbert, imagesc, mean
    #%
    #% This function simulates the normal form of a cusp bifurcation:
    #%
    #% x_dot = b1 + b2*x - x^3
    #%
    #% where b2 is a control parameter and b1 > 0.
    #%
    #% Here we simulate both the deterministic and stochastic cases for the
    #% fold bifurcation . 
    #%
    #% [Xdet, Xsto, Fext] = cuspstoch(b1,b2,xNoiseSTD,tvec)
    #%
    #% Xdet : deterministic result
    #% Xsto : stochastic result
    #% Fext : external force
    #% 
    #% tvec : time vector
    #% mu : control parameter, injected current
    #% xNoiseSTD : standard deviation of stochastic noise in theta
    #%
    #% By modifying the code, you can also add a step function or external
    #% forcing. 
    #%
    #% jsalvi@rockefeller.edu
    #%
    #% Initial condition
    xzero = 1.
    #% Add external forcing if desired
    sinusoidalstim = 0.
    pulsestim = 0.
    #% pulse or sinusoid?
    Fextmax = 0.
    #% amplitude of sinusoidal stim OR pulse
    fr = 5.
    #% frequency of stimulation
    pulsestart = 400.
    #% start of pulse
    pulseend = 410.
    #% end of pulse
    #% Decrease time step size by factor of Dtfac to ensure convergence
    Dtfac = 10.**2.
    Dt = matdiv(tvec[1]-tvec[0], Dtfac)
    N = matdiv(tvec[int(0)-1], Dt)
    #%Set the default random number stream
    RandStream.setGlobalStream[int(RandStream[int('mt19937ar')-1,int('seed')-1,0])-1]
    dW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    xdet = np.zeros(1., N)
    xdet[0] = xzero
    xsto = np.zeros(1., N)
    xsto[0] = xzero
    #% External forcing
    if sinusoidalstim == 1.:
        Ftime = np.linspace(tvec[0], tvec[int(0)-1], N)
        Fext = np.dot(Fextmax, np.cos(np.dot(np.dot(2.*np.pi, fr), Ftime)))
    elif pulsestim == 1.:
        Ftime = np.linspace(tvec[0], tvec[int(0)-1], N)
        Fext = np.dot((Ftime<pulseend)-(Ftime<pulsestart), Fextmax)
        
    else:
        Fext = np.zeros(1., N)
        
    
    #% Euler-Murayama Method with Ito Integration
    for j in np.arange(2., (N)+1):
        #%Deterministic integral
        
    Xdet = np.zeros(1., (length(tvec)-1.))
    Xsto = np.zeros(1., (length(tvec)-1.))
    #%Return vectors at times specified by Time.
    Xdet[0,:] = xdet[0:N:Dtfac]
    Xsto[0,:] = xsto[0:N:Dtfac]
    Fext = Fext[0:N:Dtfac]
    Xdet = Xdet-np.mean(Xdet)
    Xsto = Xsto-np.mean(Xsto)
    #% Make a plot of the data?
    plotyn = 0.
    if plotyn == 1.:
        plt.figure
        plt.subplot(1., 3., 1.)
        plt.hold(on)
        plt.plot(tvec[1:], Xsto, 'r')
        plt.plot(tvec[1:], Xdet, 'k')
        plt.title('Black=deterministic; Red=stochastic')
        plt.subplot(1., 3., 2.)
        plt.hold(on)
        plt.plot(np.real(hilbert(Xdet)), np.imag(hilbert(Xdet)), 'k')
        plt.plot(np.real(hilbert(Xsto)), np.imag(hilbert(Xsto)), 'r')
        plt.title('Black=deterministic; Red=stochastic')
        plt.subplot(1., 3., 3.)
        [bw, dens] = kde2d(np.array(np.hstack((np.real(hilbert(Xsto)).conj().T, np.imag(hilbert(Xsto)).conj().T))))
        imagesc(dens)
        plt.title('Stochastic Phase Space Density')
        #%figure;
        #%plot(tvec(1:length(Fext)),Fext);
    
    
    return [Xdet, Xsto, Fext]