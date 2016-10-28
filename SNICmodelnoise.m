function [Xdet, Xsto, Fext] = SNICmodelnoise(mu,NoiseSTD,tvec,fr,Fextmax)
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
%Fextmax = 1;        % amplitude of sinusoidal stim OR pulse



% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^0;
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

Ftvec = linspace(tvec(1),tvec(end),N);
Fext = zeros(1,N);
Fextdet = zeros(1,N);
Fextsto = zeros(1,N);
if sinusoidalstim == 1
    Fext = Fextmax*cos(2*pi*fr*Ftvec);
    Fextsto = Fext;
    Fextdet = Fext;
%Fext = abs(Fextmax*sawtooth(2*pi*fr*Ftvec)) - mean(abs(Fextmax*sawtooth(2*pi*fr*Ftvec)));
%Fext = Fextmax.*Ftvec;
end
pulsestim = 0;
% for a pulse stimulus, fr should be [t1 t2], in which:
% t1: starting time of pulse
% t2: ending time of pulse
if pulsestim==1
    t1n = findnearest(Ftvec,fr(1));t2n = findnearest(Ftvec,fr(2));t1n=t1n(1);t2n=t2n(1);
    Fext(t1n:t2n) = Fextmax;
    Fextsto = Fext;
    Fextdet = Fext;
end
hbmodelstim = 0;
Fc = fr;
k = Fextmax;
noiselevel = NoiseSTD * 2;

Fextsto = randn(1,N);
Fextdet = Fextsto;

if hbmodelstim == 1
    [Xdet, Xsto] = hbtoymodel(Fc,k,noiselevel,0.2,0.2,Ftvec);
    if max(abs(Xdet(1,:))) > 1
        Fextdet(2:end) = abs((Xdet(1,:)-mean(Xdet(1,:)))./std(Xdet(1,:)));
        Fextsto(2:end) = abs((Xsto(1,:)-mean(Xsto(1,:)))./std(Xsto(1,:)));
    else
        Fextdet(2:end) = abs((Xdet(1,:)-mean(Xdet(1,:))));
        Fextsto(2:end) = abs((Xsto(1,:)-mean(Xsto(1,:))));
    end
end
% Euler-Murayama Method with Ito Integration
for j = 2:N
%Deterministic integral
rdet(j-1) = sqrt(xdet(j-1)^2+ydet(j-1)^2);
xdet(j) = xdet(j-1) + Dt*( xdet(j-1)*(1 - rdet(j-1)^2) - ydet(j-1)*(1 - xdet(j-1)/rdet(j-1) + mu*(1 + xdet(j-1)/rdet(j-1)) ) + Fextdet(j));
ydet(j) = ydet(j-1) + Dt*( ydet(j-1)*(1 - rdet(j-1)^2) + xdet(j-1)*(1 - xdet(j-1)/rdet(j-1) + mu*(1 + xdet(j-1)/rdet(j-1)) ) );

%Stochastic integral
%thetasto(j) = thetasto(j-1) + Dt*(1 - cos(thetasto(j-1)) + (1 + cos(thetasto(j-1)))*I + Fext(j)) + thetaNoiseSTD*dW(j);
rsto(j-1) = sqrt(xsto(j-1)^2+ysto(j-1)^2);
xsto(j) = xsto(j-1) + Dt*( xsto(j-1)*(1 - rsto(j-1)^2) - ysto(j-1)*(1 - xsto(j-1)/rsto(j-1) + mu*(1 + xsto(j-1)/rsto(j-1)) ) + Fextsto(j)) + xNoiseSTD*xdW(j);
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
close all
if plotyn==1
    figure;
    subplot(1,2,1);hold on;plot(tvec(2:end),Xsto(1,:),'r');plot(tvec(2:end),Xdet(1,:),'k');title('Black=deterministic; Red=stochastic');
    subplot(1,2,2);hold on;plot(Xsto(1,:),imag(hilbert(Xsto(1,:))),'r');plot(Xdet(1,:),imag(hilbert(Xdet(1,:))),'k');title('Black=deterministic; Red=stochastic');
    %figure;
    %plot(tvec(1:length(Fext)),Fext);
end


end
