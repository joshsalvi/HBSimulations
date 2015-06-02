function [Xdet, Xsto, Fext] = bistablestoch(r,xNoiseSTD,tvec)
%
% This function simulates the Langevin equation
% given by the following:
%
% x_dot = r*x - x^3 + noise
%
% where r is a control parameter. 
%
% Here we simulate both the deterministic and stochastic cases for the
% saddle node bifurcation.
%
% [Xdet, Xsto, Fext] = bistablestoch(r,xNoiseSTD,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
%
%  
% tvec : time vector
% r : control parameter
% xNoiseSTD : standard deviation of stochastic noise in x
%
% By modifying the code, you can also add a step function or external
% forcing. 
%
% jsalvi@rockefeller.edu
%

% Initial condition
xzero = -1e-14;

% Add external forcing if desired
sinusoidalstim = 0; pulsestim = 0;  % pulse or sinusoid?
Fextmax = 1;        % amplitude of sinusoidal stim OR pulse
fr = 5;             % frequency of stimulation
pulsestart = 1;     % start of pulse
pulseend = 2;       % end of pulse

% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^2;
Dt = (tvec(2)-tvec(1))/Dtfac;

N = tvec(end)/Dt;

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments

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
%xdet(j) = xdet(j-1) + Dt*(mu*xdet(j-1) - 2*pi*fosc*ydet(j-1) + xdet(j-1)*(xdet(j-1)^2 + ydet(j-1)^2) - 50*xdet(j-1)*(xdet(j-1)^4 + ydet(j-1)^4) + Fext(j));
%ydet(j) = ydet(j-1) + Dt*(2*pi*fosc*xdet(j-1) + mu*ydet(j-1) + ydet(j-1)*(xdet(j-1)^2 + ydet(j-1)^2) - 50*ydet(j-1)*(xdet(j-1)^4 + ydet(j-1)^4) + Fext(j));
xdet(j) = xdet(j-1) + Dt*(r*xdet(j-1) - xdet(j-1)^3 + (Fext(j)));

%Stochastic integral
%xsto(j) = xsto(j-1) + Dt*(mu*xsto(j-1) - 2*pi*fosc*ysto(j-1) + xsto(j-1)*(xsto(j-1)^2 + ysto(j-1)^2) - 50*xdet(j-1)*(xdet(j-1)^4 + ydet(j-1)^4) + Fext(j)) + xNoiseSTD*xdW(j);
%ysto(j) = ysto(j-1) + Dt*(2*pi*fosc*xsto(j-1) + mu*ysto(j-1) + ysto(j-1)*(xsto(j-1)^2 + ysto(j-1)^2) - 50*ydet(j-1)*(xdet(j-1)^4 + ydet(j-1)^4) + Fext(j)) + yNoiseSTD*ydW(j);
xsto(j) = xsto(j-1) + Dt*(r*xsto(j-1) - xsto(j-1)^3 + (Fext(j))) + xNoiseSTD*xdW(j);

end


Xdet = zeros(1,length(tvec)-1);
Xsto = zeros(1,length(tvec)-1);

%Return vectors at times specified by Time.
Xdet(1,:) = xdet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);

% Make a plot of the data?
plotyn=0;

if plotyn==1
    figure;
    subplot(1,2,1);hold on;plot(tvec(2:end),Xsto(1,:),'r');plot(tvec(2:end),Xdet(1,:),'k');title('Black=deterministic; Red=stochastic; real part only');
    subplot(1,2,2);hold on;plot(real(hilbert(Xsto(1,:))),imag(hilbert(Xsto(1,:))),'r');plot(real(hilbert(Xdet(1,:))),imag(hilbert(Xdet(1,:))),'k');title('Black=deterministic; Red=stochastic');
end
end

