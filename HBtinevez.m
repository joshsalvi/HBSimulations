function [Xdet, Xsto, Fext] = HBtinevez(fmax,S,NoiseSTD,tvec,Fextmax,stimtype,fr,plotyn)
%
% This function simulates a model of hair-bundle mechanics from Tinevez et
% al., 2007, Biophys J
%
% [Xdet, Xsto, Fext] =  HBtinevez(fmax,S,NoiseSTD,tvec,Fextmax,stimtype,fr,plotyn)
%
% --------------------------------------------------------
%                        OUTPUTS
% --------------------------------------------------------
% Xdet : deterministic solutions
%           Xdet(1,:) -> X  (m) 
%           Xdet(2,:) -> Xa (m)
%           Xdet(3,:) -> Po 
% Xsto : stochastic solutions
%           Xsto(1,:) -> X  (m)
%           Xsto(2,:) -> Xa (m)
%           Xsto(3,:) -> Po
% Fext : stimulus waveform  (N)
%
% --------------------------------------------------------
%                        INPUTS
% --------------------------------------------------------
% fmax : maximum force produced by motors        (N)
% S : strength of Ca2+ feedback
% NoiseSTD : std.dev. of additive noise
% Fextmax : force of external stimulus           (N)
% stimtype : stimulus type (0=none,1=sinusoidal,2=pulse)
% tvec : time vector                             (s)
% fr : if sinusoidal stim-> frequency            (Hz)
%      if force pulse-> [pulse_start pulse_end]  (s)
% plotyn : plot the data? (1=yes, 0=no)
%
% -----------------------
% Written by:
% Joshua D. Salvi
% jsalvi@rockefeller.edu
% -----------------------
%

% Parameters
Ns = 50;               % number of stereocilia
% fmax = 100e-12;         % max. force produced by motors (N)
% S = 2;                  % strength of Ca2+ feedback
kf = 0.4e-3;            % stimulus fiber stiffnes (N/m)
kes = 0.25e-3;             % extent spring stiffness (N/m)
ksp = 0.76e-3;             % stereociliary pivot stiffness (N/m)
kgs = 1e-3;           % gating spring stiffness (N/m)
D = 37.1e-9;             % (nm)
gam_x = .1e-6;         % drag coefficient in X (N·s/m)
gam_y = 30e-6;          % drag coefficient in Xa (N·s/m)
A = 1;                 % 1/(1+A): Po with severed TLs
z = kgs*D/Ns;           % gating force (N)
kT = 4.11e-21;          % kT energy (J)
xes = 0;                % position for which motors bear no tension (nm)

% Initial conditions
xzero = 1e-9;
yzero = -1e-9;

% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^2;
Dt = (tvec(2)-tvec(1))/Dtfac;
N = tvec(end)/Dt;

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments
ydW = sqrt(Dt)*randn(1,N); % White noise increments
xNoiseSTD = NoiseSTD;
yNoiseSTD = NoiseSTD;

% Initialize variables
xdet = zeros(1,N); xdet(1) = xzero;
xsto = zeros(1,N); xsto(1) = xzero;
ydet = zeros(1,N); ydet(1) = yzero;
ysto = zeros(1,N); ysto(1) = yzero;
Podet = zeros(1,N);
Posto = zeros(1,N);

% External forcing
if stimtype == 1
    Ftime = linspace(tvec(1),tvec(end),N);
    Fext = Fextmax*sin(2*pi*fr*Ftime);
elseif stimtype == 2
    Ftime = linspace(tvec(1),tvec(end),N);
    Fext = ((Ftime<fr(2))-(Ftime<fr(1)))*Fextmax;
else
    Fext = zeros(1,N);
end

% Euler-Murayama Method with Ito Integration
for j = 2:N
%Deterministic integral
Podet(j-1) = ( 1 + A*exp(-z*(xdet(j-1) - ydet(j-1))/(kT)) )^-1;
xdet(j) = xdet(j-1) + Dt*(-kgs/gam_x*(xdet(j-1) - ydet(j-1) - D*Podet(j-1)) - ksp/gam_x*xdet(j-1) + Fext(j)/gam_x);
ydet(j) = ydet(j-1) + Dt*(kgs/gam_y*(xdet(j-1) - ydet(j-1) - D*Podet(j-1)) - kes/gam_y*(ydet(j-1) - xes) - fmax/gam_y*(1-S*Podet(j-1))); 

%Stochastic integral
Posto(j-1) = ( 1 + A*exp(-z*(xsto(j-1) - ysto(j-1))/(kT)) )^-1;
xsto(j) = xsto(j-1) + Dt*(-kgs/gam_x*(xsto(j-1) - ysto(j-1) - D*Posto(j-1)) - ksp/gam_x*xsto(j-1) + Fext(j)/gam_x) + xNoiseSTD*xdW(j);
ysto(j) = ysto(j-1) + Dt*(kgs/gam_y*(xsto(j-1) - ysto(j-1) - D*Posto(j-1)) - kes/gam_y*(ysto(j-1) - xes) - fmax/gam_y*(1-S*Posto(j-1))) + yNoiseSTD*ydW(j); 

end

Xdet = zeros(2,length(tvec)-1);
Xsto = zeros(2,length(tvec)-1);

%Return vectors at times specified by Time.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = ydet(1:Dtfac:N);
Xdet(3,:) = Podet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Xsto(2,:) = ysto(1:Dtfac:N);
Xsto(3,:) = Posto(1:Dtfac:N);
Fext = Fext(1:Dtfac:N);

setfiguredefaults();close all;
if plotyn==1
    figure(1);    
    subplot(3,2,1);hold on;plot(tvec(1:end-1).*1e3,Xsto(1,:).*1e9,'r');plot(tvec(1:end-1).*1e3,Xdet(1,:).*1e9);legend([{'Stochastic'},{'Deterministic'}])
    xlabel('Time (ms)');ylabel('Bundle position (nm)');box on
    subplot(3,2,4);hold on;plot(Xsto(1,:).*1e9,Xsto(2,:).*1e9,'r');plot(Xdet(1,:).*1e9,Xdet(2,:).*1e9);
    xlabel('Bundle position (nm)');ylabel('Motor position (nm)');box on
    subplot(3,2,6);hold on;plot(Xsto(1,:).*1e9,Xsto(3,:),'r');plot(Xdet(1,:).*1e9,Xdet(3,:));
    xlabel('Bundle position (nm)');ylabel('Open probability');box on
    subplot(3,2,2);plot(tvec(1:end-1).*1e3,Xsto(1,:).*1e9,'r');hold on;plot(tvec(1:end-1).*1e3,Xdet(1,:).*1e9);xlabel('Time (ms)');ylabel('Bundle position (nm)');box on
    subplot(3,2,3);plot(tvec(1:end-1).*1e3,Xsto(2,:).*1e9,'r');hold on;plot(tvec(1:end-1).*1e3,Xdet(2,:).*1e9);xlabel('Time (ms)');ylabel('Motor position (nm)');box on
    subplot(3,2,5);plot(tvec(1:end-1).*1e3,Xsto(3,:),'r');hold on;plot(tvec(1:end-1).*1e3,Xdet(3,:));xlabel('Time (ms)');ylabel('Open probability');box on

end


end

function setfiguredefaults(N)
set(0,'defaultFigureColormap',jet)
if exist('N')==1
    set(0,'defaultAxesColorOrder',ametrine(N))
else
    set(0,'defaultAxesColorOrder',[0 0 0])
end
set(0,'DefaultTextFontName','Helvetica Neue')
set(0,'DefaultTextFontUnits','Points')
set(0,'defaultlinelinewidth',2)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontName','Helvetica Neue')
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultAxesGridLineStyle',':')
set(0,'defaultFigureColor','White')
set(0,'DefaultFigurePaperType','a4letter')
set(0,'defaultFigurePosition',[1200 100 700 700])
set(0,'defaultAxesColor',[1 1 1])
set(0,'defaultLineMarker','None')
end
