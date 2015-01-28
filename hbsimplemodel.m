function [Xdet, Xsto, Fext2,Po,cdet] = hbsimplemodel(Fc,k,S,noiselevel,Fextmax,fr,tvec)
%
% This function simulates the hair-buyndle model from PNAS 2012.
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

%Stochasic HB model integration

%EM Euler-Maruyama method
%Ito integral

a = 3.5;
%b > 1 has unbounded solutions
b = 0.5;
tau = 10;
ghb=1;
xzero = 1; 
fzero = 0;
d=0.17;
x0=0.3;
Pc=0.15;
tauc=5e4;
delta=0.05;
gm=280;
gp=0.75;
taum=1;
taup=375;
kp=0.6;
cmax=1;


%Decrease tvec step size by factor of Dtfac to ensure convergence
Dtfac = 10^1;
Dt = (tvec(2)-tvec(1))/Dtfac;

N = round(tvec(end)/Dt);

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments
fdW = sqrt(Dt)*randn(1,N); % White noise increments

xdet = zeros(1,N);
fdet = zeros(1,N);
xsto = zeros(1,N);
fsto = zeros(1,N);
cdet = zeros(1,N);
csto = zeros(1,N);

xdet(1) = xzero;
xsto(1) = xzero;

fdet(1) = fzero;
fsto(1) = fzero;

cdet(1) = 0.3;
csto(1) = 0.3;

%Not using FD theorem
xNoiseSTD = noiselevel; fNoiseSTD = noiselevel; % equal noise levels

Ftvec = linspace(tvec(1),tvec(end),N);
Fext = Fextmax*cos(2*pi*fr*Ftvec);
Fcin=Fc;
for j = 2:N
%Deterministic integral
xdet(j) = xdet(j-1) + Dt*(-k/ghb*xdet(j-1) + a/ghb*(xdet(j-1) - fdet(j-1)) - 1/ghb*(xdet(j-1)-fdet(j-1))^3 + Fc/ghb + Fext(j-1));
fdet(j) = fdet(j-1) + Dt*(b*xdet(j-1) - fdet(j-1)*(S*cdet(j-1)/cmax+1))/tau;
Po(j) = 0.5+0.5*xdet(j)/sqrt(1+xdet(j)^2);
cdet(j) = cdet(j-1) + Dt*(taum/tauc*gm*Po(j)-taup/tauc*gp*cdet(j-1)^2/(cdet(j-1)^2+kp^2));
%{
if j > N1*N/100 && j < N2*N/100
    Fc =Fc2;
else
    Fc=Fcin;
end
%}

%{
if j > 5*N/15 && j < 6*N/15
    cdet(j) = 1;
end
%}

%cdet(j) = cdet(j-1) + Dt*(-cdet(j-1)/tauc*(xdet(j-1)^2/delta^2-1));

%Stochastic integral
xsto(j) = xsto(j-1) + Dt*(-k/ghb*xsto(j-1) + a/ghb*(xsto(j-1) - fsto(j-1)) - 1/ghb*(xsto(j-1)-fsto(j-1))^3 + Fc/ghb)+ xNoiseSTD*xdW(j);
fsto(j) = fsto(j-1) + Dt*(b*xsto(j-1) - fsto(j-1)*(1+S*csto(j-1)/cmax))/tau + fNoiseSTD*fdW(j);
Po(j) = 0.5+0.5*xsto(j)/sqrt(1+xsto(j)^2);
csto(j) = csto(j-1) + Dt*(taum/tauc*gm*Po(j)-taup/tauc*gp*csto(j-1)^2/(csto(j-1)^2+kp^2));
%{
if j > 4*N/12 && j < 5*N/12
    csto(j) = 0;
end
if j > 7*N/12 && j < 8*N/12
    csto(j) = 1;
end
%}
end

Xdet = zeros(2,length(tvec)-1);
Xsto = zeros(2,length(tvec)-1);

%Return vectors at tvecs specified by tvec.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = fdet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Xsto(2,:) = fsto(1:Dtfac:N);
Fext2(1,:) = Fext(1:Dtfac:N);
Cdet=cdet(1:Dtfac:N);

plotyn = 0;
if plotyn == 1
figure
%plot(Ftvec(1:end),xsto,'r');
hold on
subplot(2,1,1);
%plot(Ftvec(1:end),xsto,'r');hold on;
plot(Ftvec(1:end),xdet,'k');hold on;
%plot(Ftvec(1:end),abs(hilbert(xdet)),'r');
xlabel('tvec','FontSize',24) 
axis([Ftvec(1000) Ftvec(end-1000) 1.1*min(xdet) 1.1*max(xdet)])
ylabel('x','FontSize',24,'Rotation',0,'HorizontalAlignment','right')
subplot(2,1,2);
%plot(Ftvec(1:end),csto,'r');hold on;
plot(Ftvec(1:end),cdet,'g');
xlabel('tvec','FontSize',24) 
ylabel('C','FontSize',24,'Rotation',0,'HorizontalAlignment','right')
%axis([Ftvec(1) Ftvec(end) -0.1 1.1])
%{
figure
plot(Ftvec(1:end),fsto,'g');
hold on
plot(Ftvec(1:end),fdet,'k');
xlabel('tvec','FontSize',24) 
ylabel('f','FontSize',24,'Rotation',0,'HorizontalAlignment','right')
%}
end

end
