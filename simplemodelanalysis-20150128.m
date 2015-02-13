S= [0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5];
kc=3;
t=linspace(0,7e3,3e4);
fr=linspace(0.001,0.005,2e2);
Fextmax= [0 0.05 0.1 0.2];


Xdet=cell(length(S),length(fr),length(Fextmax));
Xsto=cell(length(S),length(fr),length(Fextmax));
Po1=cell(length(S),length(fr),length(Fextmax));
Cdet1=cell(length(S),length(fr),length(Fextmax));
meanPodet=zeros(length(S),length(fr),length(Fextmax));
stdPodet=zeros(length(S),length(fr),length(Fextmax));
meanCdet=zeros(length(S),length(fr),length(Fextmax));
stdPodet=zeros(length(S),length(fr),length(Fextmax));

for j = 1:length(S)
        disp(['j = ' num2str(j)]);
    for k = 1:length(fr);
        for l = 1:length(Fextmax)
        [Xdet{j,k,l},~,~,Po1,Cdet1]=hbsimplemodel(0,kc,S(j),0.1,Fextmax(l),fr(k),t);
                 meanPodet(j,k,l) = mean(Po1(round(length(Po1)/2):end));
                 stdPodet(j,k,l)=std(Po1(round(length(Po1)/2):end));
                 meanCdet(j,k,l)= mean(Cdet1(round(length(Cdet1)/2):end));
                 stdCdet(j,k,l)= std(Cdet1(round(length(Cdet1)/2):end));
                 clear Po1 Cdet1
        end
    end
end
disp('simulation complete');
save('/Users/joshsalvi/Downloads/11.mat','Xdet','fr','S','Fextmax','t','meanPodet','stdPodet','meanCdet','stdCdet','kc','-v7.3');
%}

Fs=t(2)-t(1);
T = 1/Fs;
L = length(Xdet{1,1,1});
NFFT = (2^4)*2^nextpow2(L);
nw=10;
XsegL = floor(length(Xdet{1,1,1})/nw);
welchwin = round(XsegL);
NPSD = floor(NFFT/nw);
noverlap = 0;
winfunc = hamming(welchwin);
freq = 0.005;
Xsine = sin(2*pi*freq.*t);
[Xsinepsd,fsinepsd] = pwelch(Xsine,winfunc,noverlap,NPSD,Fs);
winpeaknorm = sqrt(max(Xsinepsd).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;
sizeX=size(Xdet);
for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        for l = 1:sizeX(3)

                clear Xstofft Xdetfft fstofft fdetfft
                trdet=Xdet{j,k,l}(1,round(L/4):end)-smooth(Xdet{j,k,l}(1,round(L/4):end),length(Xdet{j,k,l}(1,round(L/4):end))/4)';
                vardet(j,k,l)=var(trdet);
                [Xdetfft, fdetfft] = pwelch(trdet,winfunc,noverlap,NPSD,Fs);
                fdetind = findnearest(fdetfft ,fr(k));fdetind=fdetind(1);
                fscale = 1e3;
                Xdetfft = Xdetfft./fscale;
                
                % Is the peak a minimum height? If not, set ampl/freq to zero.
            
             
                Xdetfftmaxind = fdetind;
                 if isempty(Xdetfftmaxind) ==0     
                fftampldet(j,k,l) = Xdetfft(fdetind);
                fftampldt(j,k,l) = (sqrt(fscale.*fftampldet(j,k,l).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL)./winpeaknorm;
                fftfreqdet(j,k,l) = fdetfft(fdetind);            
                if  fftampldet(j,k,l) < 1e-3
                    fftfreqdet(j,k,l) = 0;
                    fftampldet(j,k,l) = 0;
                end
                 else
                     fftfreqdet(j,k,l) = 0;
                    fftampldet(j,k,l) = 0;
                 end
    end
    end
end


for j = 1:length(S)
    S2{j}=num2str(S(j));
end

set(0,'DefaultAxesColorOrder',jet(sizeX(1)));
for j = 1:sizeX(3)
    figure(j);
    for k = 1:sizeX(1)
        plot(fr,((-fftampldt(k,:,j)+fftampldt(k,:,1))).*(1+sign((-fftampldt(k,:,j)+fftampldt(k,:,1)))));hold all;
        xlabel('Frequency (Hz)');ylabel('X-Xosc (nm)');
        title([' Fextmax = ' num2str(Fextmax(j)) ' k=3' ' Fc = 0']);
        legend(S2);
    end
end
for j = 1:sizeX(3)
    figure(j+sizeX(3));
    for k = 1:sizeX(1)
        plot(fr,meanPodet(k,:,j));hold all;
        xlabel('Frequency (Hz)');ylabel('Mean Po');
        title([' Fextmax = ' num2str(Fextmax(j)) ' k=3' ' Fc = 0']);
        legend(S2);
    end
end

for j = 1:sizeX(3)
    figure(j+2*sizeX(3));
    for k = 1:sizeX(1)
        plot(fr,meanCdet(k,:,j));hold all;
        xlabel('Frequency (Hz)');ylabel('Mean C');
        title([' Fextmax = ' num2str(Fextmax(j)) ' k=3' ' Fc = 0']);
        legend(S2);
    end
end
