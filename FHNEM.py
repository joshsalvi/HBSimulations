
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def FHNEM(Time):

    # Local Variables: xNoiseSTD, tau, fdet, Dtfac, xem, fzero, fNoiseSTD, fem, Fc, Dt, plotcheck, xzero, N, Xdet, fdW, Time, a, b, Xem, k, j, xdW, xdet, gamma
    # Function calls: plot, randn, figure, FHNEM, sqrt, RandStream, zeros, numel, ylabel, xlabel, find, mod
    #%Stochasic HB model integration
    #%EM Euler-Maruyama method
    #%Ito integral
    a = 4.
    #%b > 1 has unbounded solutions
    b = 0.1
    tau = 20.
    Fc = 0.
    gamma = 1.
    k = 3.5
    xzero = 1.
    fzero = 1.
    #%Decrease time step size by factor of Dtfac to ensure convergence
    Dtfac = 10.
    Dt = matdiv(Time[1]-Time[0], Dtfac)
    N = np.dot(numel(Time), Dtfac)
    #%Set the default random number stream
    RandStream.setGlobalStream[int(RandStream[int('mt19937ar')-1,int('seed')-1,0])-1]
    xdW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    fdW = np.dot(np.sqrt(Dt), plt.randn(1., N))
    #% White noise increments
    xdet = np.zeros(1., N)
    fdet = np.zeros(1., N)
    xem = np.zeros(1., N)
    fem = np.zeros(1., N)
    xdet[0] = xzero
    xem[0] = xzero
    fdet[0] = fzero
    fem[0] = fzero
    #%Not using FD theorem
    xNoiseSTD = 0.1
    fNoiseSTD = 0.1
    for j in np.arange(2., (N)+1):
        #%Deterministic integral
        
    Xdet = np.zeros(2., matdiv(N, Dtfac))
    Xem = np.zeros(2., matdiv(N, Dtfac))
    #%Return vectors at times specified by Time.
    Xdet[0,:] = xdet[int(nonzero((np.mod(np.arange(1., (N)+1), Dtfac) == 0.)))-1]
    Xdet[1,:] = fdet[int(nonzero((np.mod(np.arange(1., (N)+1), Dtfac) == 0.)))-1]
    Xem[0,:] = xem[int(nonzero((np.mod(np.arange(1., (N)+1), Dtfac) == 0.)))-1]
    Xem[1,:] = fem[int(nonzero((np.mod(np.arange(1., (N)+1), Dtfac) == 0.)))-1]
    plotcheck = 1.
    if plotcheck == 1.:
        plt.figure
        plt.plot(Time, Xem[0,:], 'r')
        plt.hold(on)
        plt.plot(Time, Xdet[0,:], 'k')
        plt.xlabel('Time', 'FontSize', 24.)
        plt.ylabel('x', 'FontSize', 24., 'Rotation', 0., 'HorizontalAlignment', 'right')
        plt.figure
        plt.plot(Time, Xem[1,:], 'g')
        plt.hold(on)
        plt.plot(Time, Xdet[1,:], 'k')
        plt.xlabel('Time', 'FontSize', 24.)
        plt.ylabel('f', 'FontSize', 24., 'Rotation', 0., 'HorizontalAlignment', 'right')
    
    
    return [Xdet, Xem]