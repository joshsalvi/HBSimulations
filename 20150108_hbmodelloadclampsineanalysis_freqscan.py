
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

Fs = t(2.)-t(1.)
T = 1./Fs
L = length(Xsto.cell_getattr(1., 1., 1., 1., 1.))
NFFT = np.dot(2.**6., 2.**nextpow2(L))
nw = 10.
XsegL = np.floor(matdiv(length(Xsto.cell_getattr(1., 1., 1.)), nw))
welchwin = np.round(XsegL)
NPSD = np.floor(matdiv(NFFT, nw))
noverlap = 0.
winfunc = np.hamming(welchwin)
freq = 0.005
Xsine = np.sin((np.dot(2.*np.pi, freq)*t))
[Xsinepsd, fsinepsd] = pwelch(Xsine, winfunc, noverlap, NPSD, Fs)
winpeaknorm = np.sqrt((matcompat.max(Xsinepsd)*np.dot(2., Fs)*XsegL*np.sum((np.abs(winfunc)**2.))/XsegL))/XsegL
sizeX = matcompat.size(Xsto)
for j in np.arange(1., (sizeX[0])+1):
    for k in np.arange(1., (sizeX[1])+1):
        for l in np.arange(1., (sizeX[2])+1):
            for m in np.arange(1., (sizeX[3])+1):
                for n in np.arange(1., (sizeX[4])+1):
                    trsto = Xsto.cell_getattr(j, k, l, m, n, np.arange(np.round((L/4.)), (0)+1))()-smooth(Xsto.cell_getattr(j, k, l, m, n, np.arange(np.round((L/4.)), (0)+1))(), (length(Xsto.cell_getattr(j, k, l, m, n, np.arange(np.round((L/4.)), (0)+1))())/4.)).conj().T
                    trdet = Xdet.cell_getattr(j, k, l, m, n, np.arange(np.round((L/4.)), (0)+1))()-smooth(Xdet.cell_getattr(j, k, l, m, n, np.arange(np.round((L/4.)), (0)+1))(), (length(Xdet.cell_getattr(j, k, l, m, n, np.arange(np.round((L/4.)), (0)+1))())/4.)).conj().T
                    winfunc
                    noverlap
                    NPSD
                    fstoind = findnearest(fstofft.cell_getattr(j, k, l, m, n), fr(n))
                    fstoind = fstoind[0]
                    winfunc
                    noverlap
                    NPSD
                    fdetind = findnearest(fdetfft.cell_getattr(j, k, l, m, n), fr(n))
                    fdetind = fdetind[0]
                    fscale = 1e3
                    Xstofft.cell[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = Xstofft.cell[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1]/fscale
                    Xdetfft.cell[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = Xdetfft.cell[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1]/fscale
                    #% Is the peak a minimum height? If not, set ampl/freq to zero.
                    Xstofftmaxind = 1.
                    if isempty(Xstofftmaxind) == 0.:
                        fftamplsto[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = Xstofft.cell[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1,int(fstoind[int(Xstofftmaxind[0])-1])-1]()
                        fftamplsto[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = np.sqrt((fscale*fftamplsto[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1]*np.dot(2., Fs)*XsegL*np.sum((np.abs(winfunc)**2.))/XsegL))/XsegL/winpeaknorm
                        fftfreqsto[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = fstofft.cell_getattr(j, k, l, m, n, fstoind[int(Xstofftmaxind[0])-1])()
                        fftRMSsto[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = np.std(trsto)
                    else:
                        fftfreqsto[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = 0.
                        fftamplsto[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = 0.
                        
                    
                    Xdetfftmaxind = 1.
                    if isempty(Xdetfftmaxind) == 0.:
                        fftampldet[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = Xdetfft.cell[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1,int(fdetind[int(Xdetfftmaxind[0])-1])-1]()
                        fftampldt[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = np.sqrt((fscale*fftampldet[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1]*np.dot(2., Fs)*XsegL*np.sum((np.abs(winfunc)**2.))/XsegL))/XsegL/winpeaknorm
                        fftfreqdet[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = fdetfft.cell_getattr(j, k, l, m, n, fdetind[int(Xdetfftmaxind[0])-1])()
                        fftRMSdet[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = np.std(trdet)
                    else:
                        fftfreqdet[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = 0.
                        fftampldet[int(j)-1,int(k)-1,int(l)-1,int(m)-1,int(n)-1] = 0.
                        
                    
                    
                
            
        
    