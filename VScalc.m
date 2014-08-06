function [F, VS_DetSig,VS_Sig,phi_DetSig,phi_Sig,sigSens,detSens,VS_DetSig2,VS_Sig2]=VScalc(f)

% This script uses ONHmodel.m and HBmodelEM.m to calculate vector strength
% at a given operating point
tic;
N = 50;    % Number of points to plot
Fmin = 0;   % Minimum force applied
Fmax = 200; % Maximum force applied
%f = 1/97;   % Frequency
F = Fmin:(Fmax-Fmin)/(N-1):Fmax;

ke = 100e-6;


for i = 1:N
    [Sig, DetSig, Time, Fext, Detf, Sigf] = ONHmodel(F(i),f);
    
    clear zs zd zf
    
    zs = smooth(Sig(500:end),length(Sig(500:end))/2);
    zd = smooth(DetSig(500:end),length(DetSig(500:end))/2);
    zf = smooth(Fext(500:end),length(Fext(500:end))/2);
    zdf = smooth(Detf(500:end),length(Detf(500:end))/2);
    zsf = smooth(Sigf(500:end),length(Sigf(500:end))/2);
    
    Sig(500:end) = Sig(500:end) - zs(1:length(Sig(500:end)))';
    DetSig(500:end) = DetSig(500:end) - zd(1:length(Sig(500:end)))';
    Fext(500:end) = Fext(500:end) - zf(1:length(Sig(500:end)))';
    Detf(500:end) = Detf(500:end) - zdf(1:length(Sig(500:end)))';
    Sigf(500:end) = Sigf(500:end) - zsf(1:length(Sig(500:end)))';
    
    Sighilb = hilbert(Sig(500:end));
    DetSighilb = hilbert(DetSig(500:end));
    Fexthilb = hilbert(Fext(500:end));
    Detfhilb = hilbert(Detf(500:end));
    Sigfhilb = hilbert(Sigf(500:end));
    
    Sigphi = angle(Sighilb);
    DetSigphi = angle(DetSighilb);
    Fextphi = angle(Fexthilb);
    Detfphi = angle(Detfhilb);
    Sigfphi = angle(Sigfhilb);
    
    
    phaseSig = Sigphi - Fextphi;
    phaseDetSig = DetSigphi - Fextphi;
    phaseSig2 = Sigphi - Detfphi;
    phaseDetSig2 = DetSigphi - Sigfphi;
    %phaseSig = Sigphi;
    %phaseDetSig = DetSigphi;    
    
    VS_Sig(i) = circ_r(phaseSig);
    VS_DetSig(i) = circ_r(phaseDetSig);
    VS_Sig2(i) = circ_r(phaseSig2);
    VS_DetSig2(i) = circ_r(phaseDetSig2);
    
    phi_Sig(i) = circ_mean(phaseSig);
    phi_DetSig(i) = circ_mean(phaseDetSig); 
    
    NFFT = numel(Time);
    f2 = (1000/2)*linspace(0,1,NFFT/2+1);

    actual_delta_f = 1000/numel(Time);
    %double sided
    detsigfft = fft(DetSig,NFFT)/numel(Time);
    
    q=find(f2>3 & f2<100);
    detAmp = max(abs(detsigfft(q)));
    detSens(i) = detAmp./(ke*F(i));


    actual_delta_f = 1000/numel(Time);
    %double sided
    sigfft = fft(Sig,NFFT)/numel(Time);
    
    q=find(f2>3 & f2<100);
    sigAmp = max(abs(sigfft(q)));
    sigSens(i) = sigAmp./(ke*F(i));
    %}
end

figure;
plot(F,VS_Sig,'r');hold on;plot(F,VS_DetSig,'k');
title('Vector Strength');
%{
figure;
plot(F,phi_Sig,'r');hold on;plot(F,phi_DetSig,'k');
title('Phase');
%}
figure;
plot(detSens,VS_DetSig,'k');hold on;plot(sigSens,VS_Sig,'r');
title('VS versus sensitivity');

figure;
plot(F,VS_Sig2,'r');hold on;plot(F,VS_DetSig2,'k');
title('Vector Strength 2')


toc;
