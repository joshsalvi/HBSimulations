function [Xdet, Xsto, Fext] = SNICmodelnoise(mu,NoiseSTD,tvec)
%
% This function simulates the normal form of the SNIC bifurcation
% in rectangular coordinates
%
% x_det = x/r^2*(1-r^2)-y/r^2*(1-cos(th)+mu*(1+cos(th)))
% y_det = y/r^2*(1-r^2)-x/r^2*(1-cos(th)+mu*(1+cos(th)))
% th = atan(y,x);
% r = sqrt(x^2+y^2);
%
% Here we simulate both the deterministic and stochastic cases for the
% SNIC bifurcation 
%
% [Xdet, Xsto, Fext] = SNICmodelnoise(mu,NoiseSTD,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
% 
% tvec : time vector
% mu : control parameter
% NoiseSTD : standard deviation of stochastic noise in real/imag (x/y)
%
% By modifying the code, you can also add a step function or external
% forcing. 
%
% jsalvi@rockefeller.edu
%

% Initial conditions
xzero = 1;
yzero = -1;

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
rdet(j-1) = sqrt(xdet(j-1)^2+ydet(j-1)^2);
xdet(j) = xdet(j-1) + Dt*( xdet(j-1)*(1 - rdet(j-1)^2) - ydet(j-1)*(1 - xdet(j-1)/rdet(j-1) + mu*(1 + xdet(j-1)/rdet(j-1)) ) );
ydet(j) = ydet(j-1) + Dt*( ydet(j-1)*(1 - rdet(j-1)^2) + xdet(j-1)*(1 - xdet(j-1)/rdet(j-1) + mu*(1 + xdet(j-1)/rdet(j-1)) ) );

%Stochastic integral
%thetasto(j) = thetasto(j-1) + Dt*(1 - cos(thetasto(j-1)) + (1 + cos(thetasto(j-1)))*I + Fext(j)) + thetaNoiseSTD*dW(j);
rsto(j-1) = sqrt(xsto(j-1)^2+ysto(j-1)^2);
xsto(j) = xsto(j-1) + Dt*( xsto(j-1)*(1 - rsto(j-1)^2) - ysto(j-1)*(1 - xsto(j-1)/rsto(j-1) + mu*(1 + xsto(j-1)/rsto(j-1)) ) ) + xNoiseSTD*xdW(j);
ysto(j) = ysto(j-1) + Dt*( ysto(j-1)*(1 - rsto(j-1)^2) + xsto(j-1)*(1 - xsto(j-1)/rsto(j-1) + mu*(1 + xsto(j-1)/rsto(j-1)) ) ) + yNoiseSTD*ydW(j);

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
    subplot(1,2,1);hold on;plot(tvec(2:end),Xsto(1,:),'r');plot(tvec(2:end),Xdet(1,:),'k');title('Black=deterministic; Red=stochastic');
    subplot(1,2,2);hold on;plot(Xsto(1,:),imag(hilbert(Xsto(1,:))),'r');plot(Xdet(1,:),imag(hilbert(Xdet(1,:))),'k');title('Black=deterministic; Red=stochastic');
    %figure;
    %plot(tvec(1:length(Fext)),Fext);
end


end
