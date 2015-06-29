function [Xdet, Xsto, Fext] = hopfstoch(mu,fosc,NoiseSTD,tvec,Fextmax,fr)
%
% This function simulates the normal form of the supercritical Hopf
% bifurcation, given by two planar equations:
%
% x_dot = mu*x - omega*y - x*(x^2 + y^2)
% y_dot = omega*x + mu*y - y*(x^2 + y^2)
%
% where mu is the control parameter. For mu>0, the system will oscillate at
% an amplitude that grows with sqrt(mu). Alternatively, one may express the
% above equations in polar coordinates, making the amplitude relationship
% with respect to mu more apparent:
%
% rho_dot = rho*(mu + i*omega - rho^2)
%
% Here we simulate both the deterministic and stochastic cases for the
% supercritical Hopf bifurcation 
%
% [Xdet Xsto] = hopfstoch(mu,fosc,NoiseSTD,tvec,Fextmax,fr)
%
% Xdet : deterministic result
% Xsto : stochastic result
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
xzero = 1;
yzero = -1;

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
xdW = sqrt(Dt)*randn(1,N); % White noise increments
ydW = sqrt(Dt)*randn(1,N); % White noise increments
xNoiseSTD = NoiseSTD;
yNoiseSTD = NoiseSTD;

xdet = zeros(1,N); xdet(1) = xzero;
xsto = zeros(1,N); xsto(1) = xzero;
ydet = zeros(1,N); ydet(1) = yzero;
ysto = zeros(1,N); ysto(1) = yzero;

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
xdet(j) = xdet(j-1) + Dt*(mu*xdet(j-1) - 2*pi*fosc*ydet(j-1) - xdet(j-1)*(xdet(j-1)^2 + ydet(j-1)^2) + real(Fext(j)));
ydet(j) = ydet(j-1) + Dt*(2*pi*fosc*xdet(j-1) + mu*ydet(j-1) - ydet(j-1)*(xdet(j-1)^2 + ydet(j-1)^2) + imag(Fext(j)));

%Stochastic integral
xsto(j) = xsto(j-1) + Dt*(mu*xsto(j-1) - 2*pi*fosc*ysto(j-1) - xsto(j-1)*(xsto(j-1)^2 + ysto(j-1)^2) + real(Fext(j))) + xNoiseSTD*xdW(j);
ysto(j) = ysto(j-1) + Dt*(2*pi*fosc*xsto(j-1) + mu*ysto(j-1) - ysto(j-1)*(xsto(j-1)^2 + ysto(j-1)^2) + imag(Fext(j))) + yNoiseSTD*ydW(j);

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
    subplot(1,2,1);hold on;plot(tvec(2:end),Xsto(1,:),'r');plot(tvec(2:end),Xdet(1,:),'k');title('Black=deterministic; Red=stochastic; real part only');
    subplot(1,2,2);hold on;plot(Xsto(1,:),Xsto(2,:),'r');plot(Xdet(1,:),Xdet(2,:),'k');title('Black=deterministic; Red=stochastic');
end


end
