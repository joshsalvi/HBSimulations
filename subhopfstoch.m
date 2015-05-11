function [Xdet, Xsto, Fext] = subhopfstoch(mu,omega,xNoiseSTD,yNoiseSTD,Fextmax,fr,tvec)
%
% This function simulates the normal form of the subcritical Hopf
% bifurcation, given by two polar equations:
%
% r_dot = mu*r + r^3 - r^5
% theta_dot = omega
%
% where mu is the control parameter. For mu>0, the system will oscillate at
% a high amplitude. 
%
% Here we simulate both the deterministic and stochastic cases for the
% subcritical Hopf bifurcation.
%
% [Xdet, Xsto, Fext] = subhopfstoch(mu,omega,xNoiseSTD,yNoiseSTD,Fextmax,fr,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
%
%  
% tvec : time vector
% mu : control parameter
% xNoiseSTD : standard deviation of stochastic noise in x (mag.)
% yNoiseSTD : standard deviation of stochastic noise in y (phase)
% omega : theta parameter
%
% By modifying the code, you can also add a step function or external
% forcing. 
%
% jsalvi@rockefeller.edu
%

% Initial condition
xzero = -1e-14;
yzero = 1e-14;

% Add external forcing if desired
sinusoidalstim = 1; pulsestim = 0;  % pulse or sinusoid?
%Fextmax = 0;        % amplitude of sinusoidal stim OR pulse
%fr = 5;             % frequency of stimulation
pulsestart = 500;     % start of pulse
pulseend = 500.01;       % end of pulse
fr=fr/(2*pi);

% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^2;
Dt = (tvec(2)-tvec(1))/Dtfac;

N = tvec(end)/Dt;

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments
ydW = sqrt(Dt)*randn(1,N); % White noise increments


xdet = zeros(1,N); xdet(1) = xzero;
xsto = zeros(1,N); xsto(1) = xzero;
ydet = zeros(1,N); ydet(1) = yzero;
ysto = zeros(1,N); ysto(1) = yzero;

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
xdet(j) = xdet(j-1) + Dt*(mu*xdet(j-1) + xdet(j-1)^3 - xdet(j-1)^5 + abs(Fext(j)));
ydet(j) = ydet(j-1) + Dt*(omega + xdet(j-1)^2 + angle(Fext(j)));

%Stochastic integral
%xsto(j) = xsto(j-1) + Dt*(mu*xsto(j-1) - 2*pi*fosc*ysto(j-1) + xsto(j-1)*(xsto(j-1)^2 + ysto(j-1)^2) - 50*xdet(j-1)*(xdet(j-1)^4 + ydet(j-1)^4) + Fext(j)) + xNoiseSTD*xdW(j);
%ysto(j) = ysto(j-1) + Dt*(2*pi*fosc*xsto(j-1) + mu*ysto(j-1) + ysto(j-1)*(xsto(j-1)^2 + ysto(j-1)^2) - 50*ydet(j-1)*(xdet(j-1)^4 + ydet(j-1)^4) + Fext(j)) + yNoiseSTD*ydW(j);
xsto(j) = xsto(j-1) + Dt*(mu*xsto(j-1) + xsto(j-1)^3 - xsto(j-1)^5 + abs(Fext(j))) + xNoiseSTD*xdW(j);
ysto(j) = ysto(j-1) + Dt*(omega + xsto(j-1)^2 + angle(Fext(j))) + yNoiseSTD*ydW(j);

end


Xdet = zeros(2,length(tvec)-1);
Xsto = zeros(2,length(tvec)-1);

%Return vectors at times specified by Time.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = ydet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Xsto(2,:) = ysto(1:Dtfac:N);
Fext = Fext(1:Dtfac:N);

% Make a plot of the data?
plotyn=0;

if plotyn==1
    figure;
    subplot(1,2,1);hold on;plot(tvec(2:end),Xsto(1,:).*sin(Xsto(2,:)),'r');plot(tvec(2:end),Xdet(1,:).*sin(Xdet(2,:)),'k');title('Black=deterministic; Red=stochastic');
    subplot(1,2,2);hold on;plot(Xsto(1,:).*cos(Xsto(2,:)) + 1i*Xsto(1,:).*sin(Xsto(2,:)),'r');plot(Xdet(1,:).*cos(Xdet(2,:)) + 1i*Xdet(1,:).*sin(Xdet(2,:)),'k');title('Black=deterministic; Red=stochastic');
end


end
