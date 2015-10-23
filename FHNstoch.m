function [Vdet, Vsto, Fext] = FHNstoch(I,NoiseSTD,tvec,Fextmax,fr)
%
% This function simulates the FitzHugh-Nagumo model:
%
% v_dot = v - v^3/3 - w + I
% w_dot = 0.08*(v + 0.7 - 0.8*w)
%
% where I is an injected current and the control parameter. 
% V is the membrane potential.
% W is a recovery variable.
%
% Here we simulate both the deterministic and stochastic cases for the
% supercritical Hopf bifurcation 
%
% [Vdet, Vsto, Fext] = FHNstoch(I,NoiseSTD,tvec,Fextmax,fr)
%
% Vdet : deterministic result (V(1,:)=v, V(2,:)=w)
% Vsto : stochastic result
%
%  
% tvec : time vector
% mu : control parameter
% xNoiseSTD : standard deviation of stochastic noise in x
% yNoiseSTD : standard deviation of stochastic noise in y
% fosc : frequency of oscillation on the unstable side of the bifurcation
% Fextmax : amplitude of periodic forcing
% fr : frequency of driving
%
% By modifying the code, you can also add a step function or external
% forcing. 
%
% jsalvi@rockefeller.edu
%

% Initial condition
vzero = 0.1;
wzero = -0.1;

% Add external forcing if desired
sinusoidalstim = 1; pulsestim = 0;  % pulse or sinusoid?
%Fextmax = 1e8;        % amplitude of sinusoidal stim OR pulse
%fr = 1.01;             % frequency of stimulation
pulsestart = 1;     % start of pulse
pulseend = 2;       % end of pulse

% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^2;
Dt = (tvec(2)-tvec(1))/Dtfac;

N = tvec(end)/Dt;

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
vdW = sqrt(Dt)*randn(1,N); % White noise increments
wdW = sqrt(Dt)*randn(1,N); % White noise increments
vNoiseSTD = NoiseSTD;
wNoiseSTD = NoiseSTD;

vdet = zeros(1,N); vdet(1) = vzero;
vsto = zeros(1,N); vsto(1) = vzero;
wdet = zeros(1,N); wdet(1) = wzero;
wsto = zeros(1,N); wsto(1) = wzero;

% External forcing
if sinusoidalstim == 1
    Ftime = linspace(tvec(1),tvec(end),N);
    %Fext = Fextmax*sin(2*pi*fr*Ftime);
    Fext = Fextmax*sawtooth(2*pi*fr*Ftime);
elseif pulsestim == 1
    Ftime = linspace(tvec(1),tvec(end),N);
    Fext = ((Ftime<pulseend)-(Ftime<pulsestart))*Fextmax;
else
    Fext = zeros(1,N);
end

% Euler-Murayama Method with Ito Integration
for j = 2:N
%Deterministic integral
vdet(j) = vdet(j-1) + Dt*( vdet(j-1) - vdet(j-1)^3/3 - wdet(j-1) + I ) + real(Fext(j));
wdet(j) = wdet(j-1) + Dt*( 0.08*(vdet(j-1) + 0.7 - 0.8*wdet(j-1)) );

%Stochastic integral
vsto(j) = vsto(j-1) + Dt*( vsto(j-1) - vsto(j-1)^3/3 - wsto(j-1) + I ) + real(Fext(j)) + wNoiseSTD*wdW(j);
wsto(j) = wsto(j-1) + Dt*( 0.08*(vsto(j-1) + 0.7 - 0.8*wsto(j-1)) ) + vNoiseSTD*vdW(j);

end


Vdet = zeros(2,length(tvec)-1);
Vsto = zeros(2,length(tvec)-1);

%Return vectors at times specified by Time.
Vdet(1,:) = vdet(1:Dtfac:N);
Vdet(2,:) = wdet(1:Dtfac:N);
Vsto(1,:) = vsto(1:Dtfac:N);
Vsto(2,:) = wsto(1:Dtfac:N);
Fext = Fext(1:Dtfac:N);

% Make a plot of the data?
plotyn=1;

if plotyn==1
    figure;
    subplot(1,2,1);hold on;plot(tvec(2:end),Vsto(1,:),'r');plot(tvec(2:end),Vdet(1,:),'k');title('Black=deterministic; Red=stochastic; real part only');
    subplot(1,2,2);hold on;plot(Vsto(1,:),Vsto(2,:),'r');plot(Vdet(1,:),Vdet(2,:),'k');title('Black=deterministic; Red=stochastic');
end


end
