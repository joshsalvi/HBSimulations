function [Xdet, Xsto,Podet] = nadrowskimodel(S,fmax,noiselevel,tvec)
%
% This function simulates a hair bundle (Nadrowski et al 2004)
%
% [Xdet, Xsto, Fext] = hbtoymodel(Fc,k,noiselevel,Fextmax,fr,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
%
%  
% tvec : tvec vector
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
% close all;[Xdet, Xsto,podet] = nadrowskimodel(1e3,1e-5,0.1e-9,linspace(0,1e4,1e5));


%Stochasic HB model integration

%EM Euler-Maruyama method
%Ito integral

xzero = 0.1e-9;xazero=-0.05e-9; pozero=0.51;
czero = 0;

% Initial conditions
%ghb = 0.28e-6;
%ga = 10e-6;
ghb=1;ga=1;
kgs = 1e-3;
D = 40e-9;
ksp = 0.6e-3;
kes = 0.25e-3;
A = 1;
xes = 8e-12;
kT = 4e-21;

%Decrease tvec step size by factor of Dtfac to ensure convergence
Dtfac = 10^0;
Dt = (tvec(2)-tvec(1))/Dtfac;

N = round(tvec(end)/Dt);

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments
xadW = sqrt(Dt)*randn(1,N); % White noise increments
fdW = sqrt(Dt)*randn(1,N); % White noise increments

xdet = zeros(1,N);
cdet = zeros(1,N);
xsto = zeros(1,N);
csto = zeros(1,N);

xdet(1) = xzero;
xadet(1) = xazero;
xsto(1) = xzero;
xasto(1) = xazero;
Podet(1) = pozero;
Posto(1) = pozero;

cdet(1) = czero;
csto(1) = czero;

%Not using FD theorem
xNoiseSTD = noiselevel; xaNoiseSTD = noiselevel; cNoiseSTD = noiselevel; % equal noise levels

Ftvec = linspace(tvec(1),tvec(end),N);
Fext = 0*cos(2*pi*1*Ftvec);

for j = 2:N
%Deterministic integral
xdet(j) = xdet(j-1) + Dt*( -kgs/ghb*(xdet(j-1) - xadet(j-1) - D*Podet(j-1)) - ksp/ghb*xdet(j-1) + 1/ghb*Fext(j) );
xadet(j) = xadet(j-1) + Dt*( -kgs/ga*(xdet(j-1) - xadet(j-1) - D*Podet(j-1)) - kes/ga*(xadet(j-1) - xes) - fmax/ga*(1-S*Podet(j-1)) );
Podet(j) =  (1 + A*exp((-(xdet(j-1)-xadet(j-1)))/kT))^-1 ;

%Stochastic integral
xsto(j) = xsto(j-1) + Dt*( -kgs/ghb*(xsto(j-1) - xasto(j-1) - D*Posto(j-1)) - ksp/ghb*xsto(j-1) + 1/ghb*Fext(j) ) + xNoiseSTD*xdW(j);
xasto(j) = xasto(j-1) + Dt*( -kgs/ghb*(xsto(j-1) - xasto(j-1) - D*Posto(j-1)) - kes/ghb*(xasto(j-1) - xes) - fmax/ga*(1-S*Posto(j-1)) ) + xaNoiseSTD*xadW(j);
Posto(j) = (1 + A*exp((-(xsto(j-1)-xasto(j-1)))/kT))^-1 ;

end

Xdet = zeros(2,length(tvec)-1);
Xsto = zeros(2,length(tvec)-1);


%Return vectors at tvecs specified by tvec.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = xadet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Xsto(2,:) = xasto(1:Dtfac:N);
tstart = round(7*length(Xdet)/8);
tend = length(Xdet);
plotyn = 1;
if plotyn == 1
figure
plot(Ftvec(tstart:tend),xdet(tstart:tend),'r');
hold on
plot(Ftvec(tstart:tend),xsto(tstart:tend),'k');
xlabel('tvec','FontSize',24) 
ylabel('x','FontSize',24,'Rotation',0,'HorizontalAlignment','right')

figure
plot(Ftvec(tstart:tend),xadet(tstart:tend),'g');
hold on
plot(Ftvec(tstart:tend),xasto(tstart:tend),'k');
xlabel('tvec','FontSize',24) 
ylabel('f','FontSize',24,'Rotation',0,'HorizontalAlignment','right')
end

end
