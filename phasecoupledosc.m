function [Xdet, Xsto, XdetX, XstoX] = phasecoupledosc(w1,w2,A,NoiseSTD,tvec)
%
% This code simulated two phase coupled oscillators.
%
% [Xdet, Xsto, XdetX, XstoX] = phasecoupledosc(w1,w2,A,NoiseSTD,tvec)
%
% w1,w2: angular frequencies of oscillators 1 and 2
% A: coupling coefficient
% NoiseSTD: standard deviation of the noise
% tvec: time vector
% Xdet,Xsto: deterministic and stochastic results for the phase of each
% oscillator
% XdetX,XstoX: cos(X) for each oscillator
%
% jsalvi@rockefeller.edu
%

% Initial condition
xzero = 1;
yzero = -1;

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


% Euler-Murayama Method with Ito Integration
for j = 2:N
%Deterministic integral
xdet(j) = xdet(j-1) + Dt*( w1 - A*sin(xdet(j-1) - ydet(j-1)) - (1-A)*sin(xdet(j-1) - 2*ydet(j-1)) );
ydet(j) = ydet(j-1) + Dt*( w2 + A*sin(xdet(j-1) - ydet(j-1)) + (1-A)*sin(xdet(j-1) - 2*ydet(j-1)) );

%Stochastic integral
xsto(j) = xsto(j-1) + Dt*( w1 - A*sin(xsto(j-1) - ysto(j-1)) - (1-A)*sin(xsto(j-1) - 2*ysto(j-1)) ) + xNoiseSTD*xdW(j);
ysto(j) = ysto(j-1) + Dt*( w2 + A*sin(xsto(j-1) - ysto(j-1)) + (1-A)*sin(xsto(j-1) - 2*ysto(j-1)) ) + yNoiseSTD*ydW(j);

end


Xdet = zeros(2,length(tvec)-1);
Xsto = zeros(2,length(tvec)-1);

%Return vectors at times specified by Time.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = ydet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Xsto(2,:) = ysto(1:Dtfac:N);

XdetX(1,:) = cos(Xdet(1,:));
XdetX(2,:) = cos(Xdet(2,:));
XstoX(1,:) = cos(Xsto(1,:));
XstoX(2,:) = cos(Xsto(2,:));

% Make a plot of the data?
plotyn=0;

if plotyn==1
    figure;
    subplot(1,2,1);hold on;plot(tvec(2:end),XstoX(1,:),'r');plot(tvec(2:end),XdetX(1,:),'k');title('Black=deterministic; Red=stochastic; real part only');
    subplot(1,2,2);hold on;plot(XstoX(1,:),XstoX(2,:),'r');plot(XdetX(1,:),XdetX(2,:),'k');title('Black=deterministic; Red=stochastic');
end


end
