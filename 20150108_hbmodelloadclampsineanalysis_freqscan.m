 Fs=t(2)-t(1);
T = 1/Fs;
L = length(Xsto{1,1,1,1,1});
NFFT = (2^6)*2^nextpow2(L);
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
                for n = 1:sizeX(5)
             
                trsto=Xsto{j,k,l,m,n}(round(L/4):end)-smooth(Xsto{j,k,l,m,n}(round(L/4):end),length(Xsto{j,k,l,m,n}(round(L/4):end))/4)';
                trdet=Xdet{j,k,l,m,n}(round(L/4):end)-smooth(Xdet{j,k,l,m,n}(round(L/4):end),length(Xdet{j,k,l,m,n}(round(L/4):end))/4)';
                [Xstofft{j,k,l,m,n}, fstofft{j,k,l,m,n}]= pwelch(trsto,winfunc,noverlap,NPSD,Fs);
                fstoind = findnearest(fstofft{j,k,l,m,n},fr(n));fstoind=fstoind(1);
                [Xdetfft{j,k,l,m,n}, fdetfft{j,k,l,m,n}] = pwelch(trdet,winfunc,noverlap,NPSD,Fs);
                fdetind = findnearest(fdetfft{j,k,l,m,n},fr(n));fdetind=fdetind(1);
                fscale = 1e3;
                Xstofft{j,k,l,m,n} = Xstofft{j,k,l,m,n}./fscale;
                Xdetfft{j,k,l,m,n} = Xdetfft{j,k,l,m,n}./fscale;
                
                % Is the peak a minimum height? If not, set ampl/freq to zero.
                
                Xstofftmaxind = 1;
                if isempty(Xstofftmaxind) ==0
                fftamplsto(j,k,l,m,n) = Xstofft{j,k,l,m,n}(fstoind(Xstofftmaxind(1)));
                fftamplsto(j,k,l,m,n) = (sqrt(fscale.*fftamplsto(j,k,l,m,n).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;        
                fftfreqsto(j,k,l,m,n) = fstofft{j,k,l,m,n}(fstoind(Xstofftmaxind(1)));           
                fftRMSsto(j,k,l,m,n) = std(trsto);
                else
                    fftfreqsto(j,k,l,m,n) = 0;
                    fftamplsto(j,k,l,m,n) = 0;
                end
                Xdetfftmaxind = 1;
                 if isempty(Xdetfftmaxind) ==0
                fftampldet(j,k,l,m,n) = Xdetfft{j,k,l,m,n}(fdetind(Xdetfftmaxind(1)));
                fftampldt(j,k,l,m,n) = (sqrt(fscale.*fftampldet(j,k,l,m,n).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
                fftfreqdet(j,k,l,m,n) = fdetfft{j,k,l,m,n}(fdetind(Xdetfftmaxind(1)));            
                fftRMSdet(j,k,l,m,n) = std(trdet);
                 else
                     fftfreqdet(j,k,l,m,n) = 0;
                    fftampldet(j,k,l,m,n) = 0;
                 end
                 
    end
    end
end
end
end
