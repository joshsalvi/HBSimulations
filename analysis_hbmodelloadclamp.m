evmv=0;
kind=1;
if evmv==1
    Xdet=Xdet{kind,:,:,:,:};
    Xsto=Xsto{kind,:,:,:,:};
end
    
Xdet=squeeze(Xdet);
Xsto=squeeze(Xsto);

Fs=t(2)-t(1);
T = 1/Fs;
L = length(Xsto{1,1,1,1});
NFFT = (2^4)*2^nextpow2(L);
nw=10;
XsegL = floor(length(Xsto{1,1,1})/nw);
welchwin = round(XsegL);
NPSD = floor(NFFT/nw);
noverlap = 0;
winfunc = hamming(welchwin);
freq = 0.005;
Xsine = sin(2*pi*freq.*t);
[Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
sizeX=size(Xsto);

for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        for l = 1:sizeX(3)
            for m = 1:sizeX(4)
                
                trsto=Xsto{j,k,l,m}(round(L/4):end)-smooth(Xsto{j,k,l,m}(round(L/4):end),length(Xsto{j,k,l,m}(round(L/4):end))/4)';
                trdet=Xdet{j,k,l,m}(round(L/4):end)-smooth(Xdet{j,k,l,m}(round(L/4):end),length(Xdet{j,k,l,m}(round(L/4):end))/4)';
                [Xstofft{j,k,l,m}, fstofft{j,k,l,m}]= pwelch(trsto,winfunc,noverlap,NPSD,Fs);
                fstoind = find(fstofft{j,k,l,m} > 1e-10);
                [Xdetfft{j,k,l,m}, fdetfft{j,k,l,m}] = pwelch(trdet,winfunc,noverlap,NPSD,Fs);
                fdetind = find(fdetfft{j,k,l,m} > 1e-10);
                fscale = 1e3;
                Xstofft{j,k,l,m} = Xstofft{j,k,l,m}./fscale;
                Xdetfft{j,k,l,m} = Xdetfft{j,k,l,m}./fscale;
                
                % Is the peak a minimum height? If not, set ampl/freq to zero.
            
                Xstofftmaxind = find(Xstofft{j,k,l,m}(fstoind)==max(Xstofft{j,k,l,m}(fstoind)));
                fftamplsto(j,k,l,m) = Xstofft{j,k,l,m}(fstoind(Xstofftmaxind(1)));
                fftamplsto(j,k,l,m) = (sqrt(fscale.*fftamplsto(j,k,l,m).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;        
                fftfreqsto(j,k,l,m) = fstofft{j,k,l,m}(fstoind(Xstofftmaxind(1)));           
                if  fftamplsto(j,k,l,m) < 1e-3
                    fftfreqsto(j,k,l,m) = 0;
                    fftamplsto(j,k,l,m) = 0;
                end
                
                Xdetfftmaxind = find(Xdetfft{j,k,l,m}(fdetind)==max(Xdetfft{j,k,l,m}(fdetind)));
                fftampldet(j,k,l,m) = Xdetfft{j,k,l,m}(fdetind(Xdetfftmaxind(1)));
                fftampldt(j,k,l,m) = (sqrt(fscale.*fftampldet(j,k,l,m).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
                fftfreqdet(j,k,l,m) = fdetfft{j,k,l,m}(fdetind(Xdetfftmaxind(1)));            
                if  fftampldet(j,k,l,m) < 1e-3
                    fftfreqdet(j,k,l,m) = 0;
                    fftampldet(j,k,l,m) = 0;
                end
    end
    end
end
end
