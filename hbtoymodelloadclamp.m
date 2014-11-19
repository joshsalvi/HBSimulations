function [Xdet, Xsto, Deltadet, Deltasto, Xcdet,Xcsto, Fext2,xdetdd] = hbtoymodelloadclamp(Fc,k,noiselevel,Fextmax,fr,ksf,Fe,kv,ev,mv,G,tvec)
%
% This function simulates the hair-buyndle model from PNAS 2012.
%
% [Xdet, Xsto, Delta, Fext2] = hbtoymodelloadclamp(Fc,k,noiselevel,Fextmax,fr,ksf,Fe,kv,gv,mv,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
% Deltasto/det : motion of the base of the fiber
% Fext2 : external force independent of fiber
%
%  
% tvec : tvec vector
% Fc,k : control parameters
% noiselevels : standard deviation of stochastic noise in x and y
% fr : frequency of oscillation on the unstable side of the bifurcation
% Fextmax : amplitude in force of sinusoidal stimulation.
% ksf : fiber stiffness
% kv,ev,mv : virtual impedances from fiber
% Fe : external force from fiber
% G : proportional gain of clamp
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
xzero = 1; 
fzero = 0;
Deltadet(1:3) = 0;Xcdet(1:3) = 0;
Deltasto(1:3) = 0;Xcsto(1:3) = 0;
Vc(1:2) = 0;
alpha = 10;
beta = .1;
edx = 1e-3;
exx = 1e-3;


%Decrease tvec step size by factor of Dtfac to ensure convergence
Dtfac = 1;
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

xdet(1) = xzero;
xsto(1) = xzero;

fdet(1) = fzero;
fsto(1) = fzero;

%Not using FD theorem
xNoiseSTD = noiselevel; fNoiseSTD = noiselevel; % equal noise levels

Ftvec = linspace(tvec(1),tvec(end),N);
Fext = Fextmax*cos(2*pi*fr*Ftvec);

for j = 2:N
    
if j > 3
%Load Clamp - Deterministic
xdetd(j) = diff([xdet(j-1) xdet(j-2)]);
xdetdd(j) = diff([xdet(j-1) xdet(j-2) xdet(j-3)],2);
xstod(j) = diff([xsto(j-1) xsto(j-2)]);
xstodd(j) = diff([xsto(j-1) xsto(j-2) xsto(j-3)],2);
Xcdet(j) = (edx*Xcdet(j-1))/((edx-ksf*Dt)*alpha) + Dt*((mv*xdetdd(j)+(ev-exx+alpha*beta*G*edx)*xdetd(j)+(kv-(1+alpha*beta*G)*ksf)*xdet(j-1)-Fe)/(alpha*beta*G*(edx-ksf*Dt)));
Xcsto(j) = (edx*Xcsto(j-1))/((edx-ksf*Dt)*alpha) + Dt*((mv*xstodd(j)+(ev-exx+alpha*beta*G*edx)*xstod(j)+(kv-(1+alpha*beta*G)*ksf)*xsto(j-1)-Fe)/(alpha*beta*G*(edx-ksf*Dt)));
Deltadet(j) = G*beta*alpha*(Xcdet(j)-xdet(j-1));
Deltasto(j) = G*beta*alpha*(Xcsto(j)-xsto(j-1));
%Vcdet(j) = (edx*Vcdet(j-1))/(edx-ksf*Dt) + Dt*((mv*xdetdd(j)+(ev-exx+G*edx)*xdetd(j)+(kv-(1+G)*ksf)*xdet(j)-Fe)/(G*(edx-ksf*Dt)));
%Vcsto(j) = (edx*Vcsto(j-1))/(edx-ksf*Dt) + Dt*((mv*xstodd(j)+(ev-exx+G*edx)*xstod(j)+(kv-(1+G)*ksf)*xsto(j)-Fe)/(G*(edx-ksf*Dt)));
%Deltadet(j) = (Vcdet(j)) - xdet(j-1);
%Deltasto(j) = (Vcsto(j)) - xdet(j-1);
if G == 0
    Deltadet(j) = 0;
    Deltasto(j) = 0;
end
end

%Deterministic integral
fsfdet(j) = ksf*(Deltadet(j)-xdet(j-1));
xsfdet(j) = fsfdet(j)/k;
xdet(j) = xdet(j-1) + Dt*(-k*xdet(j-1) + a*(xdet(j-1)-fdet(j-1)) - (xdet(j-1)-fdet(j-1))^3 + Fc + Fext(j) + fsfdet(j));
fdet(j) = fdet(j-1) + Dt*(b*xdet(j-1) - fdet(j-1))/tau;

%Stochastic integral
fsfsto(j) = ksf*(Deltasto(j)-xsto(j-1));
xsfsto(j) = fsfsto(j)/k;
xsto(j) = xsto(j-1) + Dt*(-k*xsto(j-1) + a*(xsto(j-1)-fsto(j-1)) - (xsto(j-1)-fsto(j-1))^3 + Fc + Fext(j) + fsfsto(j)) + xNoiseSTD*xdW(j);
fsto(j) = fsto(j-1) + Dt*(b*xsto(j-1) - fsto(j-1))/tau + fNoiseSTD*fdW(j)/tau;

end


Xdet = zeros(2,length(tvec)-1);
Xsto = zeros(2,length(tvec)-1);

%Return vectors at tvecs specified by tvec.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = fdet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Xsto(2,:) = fsto(1:Dtfac:N);
Fext2(1,:) = Fext(1:Dtfac:N);

plotyn = 1;
if plotyn == 1
figure
plot(Ftvec(1:end),xsto,'r');
hold on
plot(Ftvec(1:end),xdet,'k');
xlabel('tvec','FontSize',24) 
ylabel('x','FontSize',24,'Rotation',0,'HorizontalAlignment','right')
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
