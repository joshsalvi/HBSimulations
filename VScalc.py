
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def VScalc(f):

    # Local Variables: phaseSig, f2, Fextphi, zd, zf, detsigfft, phi_Sig, NFFT, Sighilb, DetSig, phaseDetSig2, VS_Sig2, zs, DetSigphi, VS_DetSig, VS_Sig, Fmax, Sigfhilb, DetSighilb, Sigfphi, detAmp, sigAmp, sigfft, Detfphi, Detfhilb, sigSens, VS_DetSig2, F, Sigphi, detSens, N, phi_DetSig, phaseSig2, actual_delta_f, Fmin, Time, Detf, Fexthilb, ke, Fext, i, f, Sigf, q, zsf, Sig, zdf, phaseDetSig
    # Function calls: plot, angle, figure, title, ONHmodel, max, fft, smooth, circ_mean, find, length, abs, circ_r, linspace, numel, VScalc, toc, hilbert, tic
    #% This script uses ONHmodel.m and HBmodelEM.m to calculate vector strength
    #% at a given operating point
    tic
    N = 50.
    #% Number of points to plot
    Fmin = 0.
    #% Minimum force applied
    Fmax = 200.
    #% Maximum force applied
    #%f = 1/97;   % Frequency
    F = np.arange(Fmin, (Fmax)+(matdiv(Fmax-Fmin, N-1.)), matdiv(Fmax-Fmin, N-1.))
    ke = 100e-6
    for i in np.arange(1., (N)+1):
        [Sig, DetSig, Time, Fext, Detf, Sigf] = ONHmodel(F[int(i)-1], f)
        clear(zszdzf)
        zs = smooth(Sig[499:], (length(Sig[499:])/2.))
        zd = smooth(DetSig[499:], (length(DetSig[499:])/2.))
        zf = smooth(Fext[499:], (length(Fext[499:])/2.))
        zdf = smooth(Detf[499:], (length(Detf[499:])/2.))
        zsf = smooth(Sigf[499:], (length(Sigf[499:])/2.))
        Sig[499:] = Sig[499:]-zs[0:length(Sig[499:])].conj().T
        DetSig[499:] = DetSig[499:]-zd[0:length(Sig[499:])].conj().T
        Fext[499:] = Fext[499:]-zf[0:length(Sig[499:])].conj().T
        Detf[499:] = Detf[499:]-zdf[0:length(Sig[499:])].conj().T
        Sigf[499:] = Sigf[499:]-zsf[0:length(Sig[499:])].conj().T
        Sighilb = hilbert(Sig[499:])
        DetSighilb = hilbert(DetSig[499:])
        Fexthilb = hilbert(Fext[499:])
        Detfhilb = hilbert(Detf[499:])
        Sigfhilb = hilbert(Sigf[499:])
        Sigphi = np.angle(Sighilb)
        DetSigphi = np.angle(DetSighilb)
        Fextphi = np.angle(Fexthilb)
        Detfphi = np.angle(Detfhilb)
        Sigfphi = np.angle(Sigfhilb)
        phaseSig = Sigphi-Fextphi
        phaseDetSig = DetSigphi-Fextphi
        phaseSig2 = Sigphi-Detfphi
        phaseDetSig2 = DetSigphi-Sigfphi
        #%phaseSig = Sigphi;
        #%phaseDetSig = DetSigphi;    
        VS_Sig[int(i)-1] = circ_r(phaseSig)
        VS_DetSig[int(i)-1] = circ_r(phaseDetSig)
        VS_Sig2[int(i)-1] = circ_r(phaseSig2)
        VS_DetSig2[int(i)-1] = circ_r(phaseDetSig2)
        phi_Sig[int(i)-1] = circ_mean(phaseSig)
        phi_DetSig[int(i)-1] = circ_mean(phaseDetSig)
        NFFT = numel(Time)
        f2 = np.dot(1000./2., np.linspace(0., 1., (NFFT/2.+1.)))
        actual_delta_f = 1000./numel(Time)
        #%double sided
        detsigfft = matdiv(np.fft(DetSig, NFFT), numel(Time))
        q = nonzero(np.logical_and(f2 > 3., f2<100.))
        detAmp = matcompat.max(np.abs(detsigfft[int(q)-1]))
        detSens[int(i)-1] = detAmp/np.dot(ke, F[int(i)-1])
        actual_delta_f = 1000./numel(Time)
        #%double sided
        sigfft = matdiv(np.fft(Sig, NFFT), numel(Time))
        q = nonzero(np.logical_and(f2 > 3., f2<100.))
        sigAmp = matcompat.max(np.abs(sigfft[int(q)-1]))
        sigSens[int(i)-1] = sigAmp/np.dot(ke, F[int(i)-1])
        #%}
        
    plt.figure
    plt.plot(F, VS_Sig, 'r')
    plt.hold(on)
    plt.plot(F, VS_DetSig, 'k')
    plt.title('Vector Strength')
    #%{
    plt.figure
    plt.plot(F, phi_Sig, 'r')
    plt.hold(on)
    plt.plot(F, phi_DetSig, 'k')
    plt.title('Phase')
    #%}
    plt.figure
    plt.plot(detSens, VS_DetSig, 'k')
    plt.hold(on)
    plt.plot(sigSens, VS_Sig, 'r')
    plt.title('VS versus sensitivity')
    plt.figure
    plt.plot(F, VS_Sig2, 'r')
    plt.hold(on)
    plt.plot(F, VS_DetSig2, 'k')
    plt.title('Vector Strength 2')
    toc
    return [F, VS_DetSig, VS_Sig, phi_DetSig, phi_Sig, sigSens, detSens, VS_DetSig2, VS_Sig2]