function [Xdet, Xsto, Fext2] = hbtoymodeltuned(Fc,kzero,knoise,taua,noiselevel,Fextmax,fr,tvec)
%
% This function simulates the hair-bundle model but adds tuning to the
% parameter a.
%
% [Xdet, Xsto, Fext] = hbtoymodeltuned(Fc,kzero,knoise,noiselevel,Fextmax,fr,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
%
%  
% tvec : tvec vector
% kzero : starting stiffness value
% knoise : noise in k (STD)
% taua : rate of tuning
% Fc,k : control parameters
% noiselevels : standard deviation of stochastic noise in x and y
% fr : frequency of oscillation on the unstable side of the bifurcation
% Fextmax : amplitude in force of sinusoidal stimulation.
%
% Note that stiffnesses are scaled by a factor of 100 in the manuscript.
%
% By modifying the code, you can also add a step function or other external
% forcing. 
%
% jsalvi@rockefeller.edu
%

%Stochasic HB model integration

%EM Euler-Maruyama method
%Ito integral

azero = 3.5;
%b > 1 has unbounded solutions
b = 0.5;
tau = 10;
xzero = 1; 
fzero = 0;

%Decrease tvec step size by factor of Dtfac to ensure convergence
Dtfac = 10^0;
Dt = (tvec(2)-tvec(1))/Dtfac;

N = round(tvec(end)/Dt);

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments
fdW = sqrt(Dt)*randn(1,N); % White noise increments
adW = sqrt(Dt)*randn(1,N); % White noise increments

xdet = zeros(1,N);
fdet = zeros(1,N);
xsto = zeros(1,N);
fsto = zeros(1,N);
adet = zeros(1,N);
asto = zeros(1,N);

xdet(1) = xzero;
xsto(1) = xzero;

fdet(1) = fzero;
fsto(1) = fzero;

adet(1) = azero;
asto(1) = azero;

% Populate stiffness values
k = zeros(1,N);
k(1) = kzero;
for j = 2:N
    k(j) = k(1) + knoise*randn(1);
end




%Not using FD theorem
xNoiseSTD = noiselevel; fNoiseSTD = noiselevel; aNoiseSTD = noiselevel; % equal noise levels

Ftvec = linspace(tvec(1),tvec(end),N);
Fext = Fextmax*cos(2*pi*fr*Ftvec);

for j = 2:N
%Deterministic integral
xdet(j) = xdet(j-1) + Dt*(-k(j-1)*xdet(j-1) + adet(j-1)*(xdet(j-1)-fdet(j-1)) - (xdet(j-1)-fdet(j-1))^3 + Fc + Fext(j));
fdet(j) = fdet(j-1) + Dt*(b*xdet(j-1) - fdet(j-1))/tau;
adet(j) = adet(j-1) + Dt*(1/taua*(k(j-1) - adet(j-1)));

%Stochastic integral
xsto(j) = xsto(j-1) + Dt*(-k(j-1)*xsto(j-1) + asto(j-1)*(xsto(j-1)-fsto(j-1)) - (xsto(j-1)-fsto(j-1))^3 + Fc + Fext(j)) + xNoiseSTD*xdW(j);
fsto(j) = fsto(j-1) + Dt*(b*xsto(j-1) - fsto(j-1))/tau + fNoiseSTD*fdW(j)/tau;
asto(j) = asto(j-1) + Dt*(1/taua*(k(j-1) - asto(j-1))) + aNoiseSTD*adW(j)/taua;
end

Xdet = zeros(2,length(tvec)-1);
Xsto = zeros(2,length(tvec)-1);

%Return vectors at tvecs specified by tvec.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = fdet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Xsto(2,:) = fsto(1:Dtfac:N);
Fext2(1,:) = Fext(1:Dtfac:N);

plotyn = 0;
if plotyn == 1
figure
plot(Ftvec(1:end),xsto,'r');
hold on
plot(Ftvec(1:end),xdet,'k');
xlabel('tvec','FontSize',24) 
ylabel('x','FontSize',24,'Rotation',0,'HorizontalAlignment','right')

figure
plot(Ftvec(1:end),fsto,'g');
hold on
plot(Ftvec(1:end),fdet,'k');
xlabel('tvec','FontSize',24) 
ylabel('f','FontSize',24,'Rotation',0,'HorizontalAlignment','right')
end

end
