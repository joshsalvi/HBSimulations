evmv=2; % 1: evmv only;  2: evmv with forcing
    
Xdet=squeeze(Xdet);
Xsto=squeeze(Xsto);
if evmv==0
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
                clear Xstofft Xdetfft fstofft fdetfft
                trsto=Xsto{j,k,l,m}(round(L/4):end)-smooth(Xsto{j,k,l,m}(round(L/4):end),length(Xsto{j,k,l,m}(round(L/4):end))/4)';
                trdet=Xdet{j,k,l,m}(round(L/4):end)-smooth(Xdet{j,k,l,m}(round(L/4):end),length(Xdet{j,k,l,m}(round(L/4):end))/4)';
                [Xstofft, fstofft]= pwelch(trsto,winfunc,noverlap,NPSD,Fs);
                fstoind = find(fstofft > 1e-10);
                [Xdetfft, fdetfft] = pwelch(trdet,winfunc,noverlap,NPSD,Fs);
                fdetind = find(fdetfft > 1e-10);
                fscale = 1e3;
                Xstofft = Xstofft./fscale;
                Xdetfft = Xdetfft./fscale;
                
                % Is the peak a minimum height? If not, set ampl/freq to zero.
            
                Xstofftmaxind = find(Xstofft(fstoind)==max(Xstofft(fstoind)));
                fftamplsto(j,k,l,m) = Xstofft(fstoind(Xstofftmaxind(1)));
                fftamplsto(j,k,l,m) = (sqrt(fscale.*fftamplsto(j,k,l,m).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;        
                fftfreqsto(j,k,l,m) = fstofft(fstoind(Xstofftmaxind(1)));           
                if  fftamplsto(j,k,l,m) < 1e-3
                    fftfreqsto(j,k,l,m) = 0;
                    fftamplsto(j,k,l,m) = 0;
                end
                
                Xdetfftmaxind = find(Xdetfft(fdetind)==max(Xdetfft(fdetind)));
                fftampldet(j,k,l,m) = Xdetfft(fdetind(Xdetfftmaxind(1)));
                fftampldt(j,k,l,m) = (sqrt(fscale.*fftampldet(j,k,l,m).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
                fftfreqdet(j,k,l,m) = fdetfft(fdetind(Xdetfftmaxind(1)));            
                if  fftampldet(j,k,l,m) < 1e-3
                    fftfreqdet(j,k,l,m) = 0;
                    fftampldet(j,k,l,m) = 0;
                end
    end
    end
end
end
elseif evmv==1
    Fs=t(2)-t(1);
T = 1/Fs;
L = length(Xsto{1,1,1,1,1});
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
                for n = 1:sizeX(5)
                trsto=Xsto{j,k,l,m,n}(round(L/4):end)-smooth(Xsto{j,k,l,m,n}(round(L/4):end),length(Xsto{j,k,l,m,n}(round(L/4):end))/4)';
                trdet=Xdet{j,k,l,m,n}(round(L/4):end)-smooth(Xdet{j,k,l,m,n}(round(L/4):end),length(Xdet{j,k,l,m,n}(round(L/4):end))/4)';
                [Xstofft{j,k,l,m,n}, fstofft{j,k,l,m,n}]= pwelch(trsto,winfunc,noverlap,NPSD,Fs);
                fstoind = find(fstofft{j,k,l,m,n} > 1e-10);
                [Xdetfft{j,k,l,m,n}, fdetfft{j,k,l,m,n}] = pwelch(trdet,winfunc,noverlap,NPSD,Fs);
                fdetind = find(fdetfft{j,k,l,m,n} > 1e-10);
                fscale = 1e3;
                Xstofft{j,k,l,m,n} = Xstofft{j,k,l,m,n}./fscale;
                Xdetfft{j,k,l,m,n} = Xdetfft{j,k,l,m,n}./fscale;
                
                % Is the peak a minimum height? If not, set ampl/freq to zero.
            
                Xstofftmaxind = find(Xstofft{j,k,l,m,n}(fstoind)==max(Xstofft{j,k,l,m,n}(fstoind)));
                if isempty(Xstofftmaxind) ==0
                fftamplsto(j,k,l,m,n) = Xstofft{j,k,l,m,n}(fstoind(Xstofftmaxind(1)));
                fftamplsto(j,k,l,m,n) = (sqrt(fscale.*fftamplsto(j,k,l,m,n).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;        
                fftfreqsto(j,k,l,m,n) = fstofft{j,k,l,m,n}(fstoind(Xstofftmaxind(1)));           
                if  fftamplsto(j,k,l,m,n) < 1e-3
                    fftfreqsto(j,k,l,m,n) = 0;
                    fftamplsto(j,k,l,m,n) = 0;
                end
                else
                    fftfreqsto(j,k,l,m,n) = 0;
                    fftamplsto(j,k,l,m,n) = 0;
                end
                Xdetfftmaxind = find(Xdetfft{j,k,l,m,n}(fdetind)==max(Xdetfft{j,k,l,m,n}(fdetind)));
                 if isempty(Xdetfftmaxind) ==0
                fftampldet(j,k,l,m,n) = Xdetfft{j,k,l,m,n}(fdetind(Xdetfftmaxind(1)));
                fftampldt(j,k,l,m,n) = (sqrt(fscale.*fftampldet(j,k,l,m,n).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
                fftfreqdet(j,k,l,m,n) = fdetfft{j,k,l,m,n}(fdetind(Xdetfftmaxind(1)));            
                if  fftampldet(j,k,l,m,n) < 1e-3
                    fftfreqdet(j,k,l,m,n) = 0;
                    fftampldet(j,k,l,m,n) = 0;
                end
                 else
                     fftfreqdet(j,k,l,m,n) = 0;
                    fftampldet(j,k,l,m,n) = 0;
                 end
                 
    end
    end
end
end
end
else
  Fs=t(2)-t(1);
T = 1/Fs;
L = length(Xsto{1,1,1,1,1,1});
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
                
                for n = 1:sizeX(5)
                    for o = 1:sizeX(6)
                trsto=Xsto{j,k,l,m,n,o}(round(L/4):end)-smooth(Xsto{j,k,l,m,n,o}(round(L/4):end),length(Xsto{j,k,l,m,n,o}(round(L/4):end))/4)';
                trdet=Xdet{j,k,l,m,n,o}(round(L/4):end)-smooth(Xdet{j,k,l,m,n,o}(round(L/4):end),length(Xdet{j,k,l,m,n,o}(round(L/4):end))/4)';
                [Xstofft, fstofft]= pwelch(trsto,winfunc,noverlap,NPSD,Fs);
                fstoind = find(fstofft > 1e-10);
                [Xdetfft, fdetfft] = pwelch(trdet,winfunc,noverlap,NPSD,Fs);
                fdetind = find(fdetfft > 1e-10);
                fscale = 1e3;
                Xstofft = Xstofft./fscale;
                Xdetfft = Xdetfft./fscale;
                
                % Is the peak a minimum height? If not, set ampl/freq to zero.
            
                Xstofftmaxind = find(Xstofft(fstoind)==max(Xstofft(fstoind)));
                if isempty(Xstofftmaxind) ==0
                fftamplsto(j,k,l,m,n,o) = Xstofft(fstoind(Xstofftmaxind(1)));
                fftamplsto(j,k,l,m,n,o) = (sqrt(fscale.*fftamplsto(j,k,l,m,n,o).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;        
                fftfreqsto(j,k,l,m,n,o) = fstofft(fstoind(Xstofftmaxind(1)));           
                if  fftamplsto(j,k,l,m,n,o) < 1e-3
                    fftfreqsto(j,k,l,m,n,o) = 0;
                    fftamplsto(j,k,l,m,n,o) = 0;
                end
                else
                    fftfreqsto(j,k,l,m,n,o) = 0;
                    fftamplsto(j,k,l,m,n,o) = 0;
                end
                Xdetfftmaxind = find(Xdetfft(fdetind)==max(Xdetfft(fdetind)));
                 if isempty(Xdetfftmaxind) ==0
                fftampldet(j,k,l,m,n,o) = Xdetfft(fdetind(Xdetfftmaxind(1)));
                fftampldt(j,k,l,m,n,o) = (sqrt(fscale.*fftampldet(j,k,l,m,n,o).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
                fftfreqdet(j,k,l,m,n,o) = fdetfft(fdetind(Xdetfftmaxind(1)));            
                if  fftampldet(j,k,l,m,n,o) < 1e-3
                    fftfreqdet(j,k,l,m,n,o) = 0;
                    fftampldet(j,k,l,m,n,o) = 0;
                end
                 else
                     fftfreqdet(j,k,l,m,n,o) = 0;
                    fftampldet(j,k,l,m,n,o) = 0;
                 end
                 
    end
    end
end
end
end    
end
end



%%
fftampldet2=fftampldet;fftfreqdet2=fftfreqdet;fftamplsto2=fftamplsto;fftfreqsto2=fftfreqsto;
%%
fftampldet=fftampldet2;fftfreqdet=fftfreqdet2;fftamplsto=fftamplsto2;fftfreqsto=fftfreqsto2;
%% Fc/kc

% Remove points in which a harmonic was picked up
fftampldet2=fftampldet;fftfreqdet2=fftfreqdet;fftamplsto2=fftamplsto;fftfreqsto2=fftfreqsto;
fftampldet(find(fftfreqdet>0.01))=0;
fftfreqdet(find(fftfreqdet>0.01))=0;
fftamplsto(find(fftfreqsto>0.01))=0;
fftfreqsto(find(fftfreqsto>0.01))=0;

load('blueredcm.mat');
indG = 3;
indksf = 2;

% Amplitude
figure(1);pcolor(kc,Fc,fftampldet(:,:,indG,indksf));shading flat;colormap(red1);
clim([-0.1 max(fftampldet(:,:,indG,indksf))]);
xlabel('stiffness');ylabel('force');
title('Deterministic - Amplitude');

% Amplitude
figure(2);pcolor(kc,Fc,fftfreqdet(:,:,indG,indksf));shading flat;colormap(blue1);
clim([-1.2e-6 max(fftfreqdet(:,:,indG,indksf))]);
xlabel('stiffness');ylabel('force');
title('Deterministic - Frequency');

%{
% Amplitude
figure(3);pcolor(kc,Fc,fftamplsto(:,:,indG,indksf));shading flat;colormap(red1);
title('Stochastic (noise=0.2) - Amplitude');

% Amplitude
figure(4);pcolor(kc,Fc,fftfreqsto(:,:,indG,indksf));shading flat;colormap(blue1);
title('Stochastic (noise=0.2) - Frequency');
%}

%% plot time traces for Fc/kc
indF = [16 16 16 16 24 8];
indk = [17 17 17 38 24 24];
indG = [1 1 3 3 3 3];
indksf = [1 2 2 2 2 2];
tstart=500;
g=findnearest(t,tstart);g=g(1);

figure;
for j = 1:6
    subplot(3,2,j);
    L=length(Xdet{indF(j),indk(j),indG(j),indksf(j)}(1,g:end));
    plot(t(1:L),Xsto{indF(j),indk(j),indG(j),indksf(j)}(1,g:end),'r','LineWidth',0.2);
    hold on;plot(t(1:L),Xdet{indF(j),indk(j),indG(j),indksf(j)}(1,g:end),'k','LineWidth',1);
    xlabel('time');ylabel('hb position');title(sprintf('%s%s %s%s %s%s %s%s','Fc=',num2str(Fc(indF(j))),'kc=',num2str(kc(indk(j))),'G=',num2str(G(indG(j))),'ksf=',num2str(ksf(indksf(j)))));
    axis([0 tstart -1.7 1.7])
end

%% plot time traces for Fe/kv
indF = [16 16 16 16 24 8];
indk = [17 30 36 38 24 24];
indG = [3 3 3 3 3 3];
indksf = [2 2 2 2 2 2];
tstart=500;
g=findnearest(t,tstart);g=g(1);

figure;
for j = 1:6
    subplot(3,2,j);
    L=length(Xdet{indF(j),indk(j),indG(j),indksf(j)}(1,g:end));
    plot(t(1:L),Xsto{indF(j),indk(j),indG(j),indksf(j)}(1,g:end),'r','LineWidth',0.2);
    hold on;plot(t(1:L),Xdet{indF(j),indk(j),indG(j),indksf(j)}(1,g:end),'k','LineWidth',1);
    xlabel('time');ylabel('hb position');title(sprintf('%s%s %s%s %s%s %s%s %s%s %s%s','Fc=',num2str(Fc(1)),'kc=',num2str(kc(1)),'Fe=',num2str(Fe(indF(j))),'kv=',num2str(kv(indk(j))),'G=',num2str(G(indG(j))),'ksf=',num2str(ksf(indksf(j)))));
    axis([0 tstart -1.7 1.7])
end

%% Fe/kv

% Remove points in which a harmonic was picked up
fftampldet2=fftampldet;fftfreqdet2=fftfreqdet;fftamplsto2=fftamplsto;fftfreqsto2=fftfreqsto;
fftampldet(find(fftfreqdet>0.01))=0;
fftfreqdet(find(fftfreqdet>0.01))=0;
fftamplsto(find(fftfreqsto>0.01))=0;
fftfreqsto(find(fftfreqsto>0.01))=0;

load('blueredcm.mat');
indG = 2;
indksf = 2;

% Amplitude
figure(1);pcolor(kv,Fe,fftampldet(:,:,indG,indksf));shading flat;colormap(red1);
caxis([-0.09 max(max(fftampldet(:,:,indG,indksf)))]);
xlabel('stiffness');ylabel('force');
title('Deterministic - Amplitude');

% Amplitude
figure(2);pcolor(kv,Fe,fftfreqdet(:,:,indG,indksf));shading flat;colormap(blue1);
caxis([-1e-6 max(max(fftfreqdet(:,:,indG,indksf)))]);
xlabel('stiffness');ylabel('force');
title('Deterministic - Frequency');

%{
% Amplitude
figure(3);pcolor(kc,Fc,fftamplsto(:,:,indG,indksf));shading flat;colormap(red1);
title('Stochastic (noise=0.2) - Amplitude');

% Amplitude
figure(4);pcolor(kc,Fc,fftfreqsto(:,:,indG,indksf));shading flat;colormap(blue1);
title('Stochastic (noise=0.2) - Frequency');
%}

%% plot time traces for ev/mv
indF = [1 1 1 1 1 1];
indk = ones(1,6).*3;
indG = [3 3 3 3 3 3];
indksf = [2 2 2 2 2 2];
indev = [16 31 1 16 16 16];
indmv = [16 16 16 31 1 16];
tstart=700;
g=findnearest(t,tstart);g=g(1);

figure;
for j = 1:6
    subplot(3,2,j);
    L=length(Xdet{indk(j),indev(j),indmv(j),indG(j),indksf(j)}(1,g:end));
    plot(t(1:L),Xsto{indk(j),indev(j),indmv(j),indG(j),indksf(j)}(1,g:end),'r','LineWidth',0.2);
    hold on;plot(t(1:L),Xdet{indk(j),indev(j),indmv(j),indG(j),indksf(j)}(1,g:end),'k','LineWidth',1);
    xlabel('time');ylabel('hb position');title(sprintf('%s%s %s%s %s%s %s%s %s%s %s%s','kv=',num2str(kv(indk(j))),'ev=',num2str(ev(indev(j))),'mv=',num2str(mv(indmv(j))),'G=',num2str(G(indG(j))),'ksf=',num2str(ksf(indksf(j)))));
    axis([0 t(end)-tstart -1.7 1.7])
end
