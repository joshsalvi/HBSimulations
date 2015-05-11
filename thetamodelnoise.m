function [Xdet, Xsto, Fext,rsto] = thetamodelnoise(I,thetaNoiseSTD,tvec)
%
% This function simulates the Ermentrout-Kopell canonical model, also known
% as the theta model. This is a one-dimensional model of neuron spiking,
% with spikes occuring at theta=pi. The form of the equation is also the
% normal form of the SNIC bifurcation:
%
% theta_dot = 1 - cos(theta) + (1+cos(theta))*I
%
% where I is the injected current and is a control parameter. There are no
% spikes for I<0 and the system passes a SNIC bifurcation at I=0, spiking
% for I>0. 
%
% Here we simulate both the deterministic and stochastic cases for the
% theta model. 
%
% [thetadet, thetasto, Fext] = thetamodel(I,thetaNoiseSTD,tvec)
%
% thetadet : deterministic result
% thetasto : stochastic result
% thetadet : deterministic angle
% thetasto : stochastic angle
% 
% tvec : time vector
% I : control parameter, injected current
% thetaNoiseSTD : standard deviation of stochastic noise in theta
%
% By modifying the code, you can also add a step function or external
% forcing. 
%
% jsalvi@rockefeller.edu
%

% Initial condition
thetazero = 1;

% Add external forcing if desired
sinusoidalstim = 0; pulsestim = 0;  % pulse or sinusoid?
Fextmax = 0;        % amplitude of sinusoidal stim OR pulse
fr = 5;             % frequency of stimulation
pulsestart = 400;     % start of pulse
pulseend = 410;       % end of pulse

% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^2;
Dt = (tvec(2)-tvec(1))/Dtfac;

N = tvec(end)/Dt;

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
tdW = sqrt(Dt)*randn(1,N); % White noise increments
rdW = sqrt(Dt)*randn(1,N); % White noise increments
tNoiseSTD = thetaNoiseSTD;
rNoiseSTD = thetaNoiseSTD;


thetadet = zeros(1,N); thetadet(1) = thetazero;
thetasto = zeros(1,N); thetasto(1) = thetazero;
rdet = zeros(1,N); rdet(1) = 1;
rsto = zeros(1,N); rsto(1) = 1;

% External forcing
if sinusoidalstim == 1
    Ftime = linspace(tvec(1),tvec(end),N);
    Fext = Fextmax*cos(2*pi*fr*Ftime);
elseif pulsestim == 1
    Ftime = linspace(tvec(1),tvec(end),N);
    Fext = ((Ftime<pulseend)-(Ftime<pulsestart))*Fextmax;
else
    Fext = zeros(1,N);
end

% Euler-Murayama Method with Ito Integration
for j = 2:N
%Deterministic integral
thetadet(j) = thetadet(j-1) + Dt*(1 - cos(thetadet(j-1)) + (1 + cos(thetadet(j-1)))*I + Fext(j));
rdet(j) = rdet(j-1);

%Stochastic integral
%thetasto(j) = thetasto(j-1) + Dt*(1 - cos(thetasto(j-1)) + (1 + cos(thetasto(j-1)))*I + Fext(j)) + thetaNoiseSTD*dW(j);
thetasto(j) = thetasto(j-1) + Dt*(1 - cos(thetasto(j-1)) + (1 + cos(thetasto(j-1)))*(I - thetaNoiseSTD^2/2*sin(thetasto(j-1))) + Fext(j)) + thetaNoiseSTD*(1 + cos(thetasto(j-1)))*tdW(j);
rsto(j) = rsto(j-1) + 0.1*rNoiseSTD*rdW(j);

end


thetadet1 = zeros(1,length(tvec)-1);
thetasto1 = zeros(1,length(tvec)-1);
rdet1 = zeros(1,length(tvec)-1);
rsto1 = zeros(1,length(tvec)-1);

%Return vectors at times specified by Time.
thetadet1(1,:) = thetadet(1:Dtfac:N);
rdet1(1,:) = rdet(1:Dtfac:N);
thetasto1(1,:) = thetasto(1:Dtfac:N);
rsto1(1,:) = rsto(1:Dtfac:N);
Fext = Fext(1:Dtfac:N);

% Convert these values into time traces (theta is the angle around a unit
% circle). Use thetadet and thetasto to find the phase plane of the signal
% if you so desire
Xdet = rdet1.*cos(thetadet1);
Xsto = rsto1.*cos(thetasto1);
Xdet=Xdet-mean(Xdet);
Xsto=Xsto-mean(Xsto);

% Make a plot of the data?
plotyn=0;

if plotyn==1
    figure;
    subplot(1,2,1);hold on;plot(tvec(2:end),Xsto,'r');plot(tvec(2:end),Xdet,'k');title('Black=deterministic; Red=stochastic');
    subplot(1,2,2);hold on;plot(Xsto,imag(hilbert(Xsto)),'r');plot(Xdet,imag(hilbert(Xdet)),'k');title('Black=deterministic; Red=stochastic');
    %figure;
    %plot(tvec(1:length(Fext)),Fext);
end


end
