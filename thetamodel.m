function [Xdet, Xsto, Fext] = thetamodel(I,thetaNoiseSTD,tvec)
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
% [Xdet Xsto thetadet thetasto] = thetamodel(mu,thetaNoiseSTD,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
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
dW = sqrt(Dt)*randn(1,N); % White noise increments

xdet = zeros(1,N); xdet(1) = thetazero;
xsto = zeros(1,N); xsto(1) = thetazero;

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
xdet(j) = xdet(j-1) + Dt*(1 - cos(xdet(j-1)) + (1 + cos(xdet(j-1)))*I + Fext(j));

%Stochastic integral
%xsto(j) = xsto(j-1) + Dt*(1 - cos(xsto(j-1)) + (1 + cos(xsto(j-1)))*I + Fext(j)) + thetaNoiseSTD*dW(j);
xsto(j) = xsto(j-1) + Dt*(1 - cos(xsto(j-1)) + (1 + cos(xsto(j-1)))*(I - thetaNoiseSTD^2/2*sin(xsto(j-1))) + Fext(j)) + thetaNoiseSTD*(1 + cos(xsto(j-1)))*dW(j);

end


thetadet = zeros(1,length(tvec)-1);
thetasto = zeros(1,length(tvec)-1);

%Return vectors at times specified by Time.
thetadet(1,:) = xdet(1:Dtfac:N);
thetasto(1,:) = xsto(1:Dtfac:N);
Fext = Fext(1:Dtfac:N);

% Convert these values into time traces (theta is the angle around a unit
% circle). Use thetadet and thetasto to find the phase plane of the signal
% if you so desire
Xdet = sin(thetadet);
Xsto = sin(thetasto);
Xdet=Xdet-mean(Xdet);
Xsto=Xsto-mean(Xsto);

% Make a plot of the data?
plotyn=0;

if plotyn==1
    figure;
    subplot(1,2,1);hold on;plot(tvec(2:end),Xsto,'r');plot(tvec(2:end),Xdet,'k');title('Black=deterministic; Red=stochastic');
    subplot(1,2,2);hold on;plot(exp(-1i*thetasto),'r');plot(exp(-1i*thetadet),'k');title('Black=deterministic; Red=stochastic');
    figure;
    plot(tvec(1:length(Fext)),Fext);
end


end
