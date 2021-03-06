function [Xdet, Xsto, Fext] = cuspstoch(b1,b2,xNoiseSTD,tvec)
%
% This function simulates the normal form of a cusp bifurcation:
%
% x_dot = b1 + b2*x - x^3
%
% where b2 is a control parameter and b1 > 0.
%
% Here we simulate both the deterministic and stochastic cases for the
% fold bifurcation . 
%
% [Xdet, Xsto, Fext] = cuspstoch(b1,b2,xNoiseSTD,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
% Fext : external force
% 
% tvec : time vector
% mu : control parameter, injected current
% xNoiseSTD : standard deviation of stochastic noise in theta
%
% By modifying the code, you can also add a step function or external
% forcing. 
%
% jsalvi@rockefeller.edu
%

% Initial condition
xzero = 1;

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

xdet = zeros(1,N); xdet(1) = xzero;
xsto = zeros(1,N); xsto(1) = xzero;

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
xdet(j) = xdet(j-1) + Dt*(b1 + b2*xdet(j-1) - xdet(j-1)^3 + Fext(j));

%Stochastic integral
xsto(j) = xsto(j-1) + Dt*(b1 + b2*xsto(j-1) - xsto(j-1)^3 + Fext(j)) + xNoiseSTD*dW(j);
end


Xdet = zeros(1,length(tvec)-1);
Xsto = zeros(1,length(tvec)-1);

%Return vectors at times specified by Time.
Xdet(1,:) = xdet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Fext = Fext(1:Dtfac:N);

Xdet=Xdet-mean(Xdet);
Xsto=Xsto-mean(Xsto);

% Make a plot of the data?
plotyn=0;

if plotyn==1
    figure;
    subplot(1,3,1);hold on;plot(tvec(2:end),Xsto,'r');plot(tvec(2:end),Xdet,'k');title('Black=deterministic; Red=stochastic');
    subplot(1,3,2);hold on;plot(real(hilbert(Xdet)),imag(hilbert(Xdet)),'k');plot(real(hilbert(Xsto)),imag(hilbert(Xsto)),'r');title('Black=deterministic; Red=stochastic');
    subplot(1,3,3);[bw dens]=kde2d([real(hilbert(Xsto))' imag(hilbert(Xsto))']);imagesc(dens);title('Stochastic Phase Space Density');
    %figure;
    %plot(tvec(1:length(Fext)),Fext);
end


end
