function [Xdet, Xsto, Deltadet, Deltasto,Fext,fsfdet] = hbtoymodelloadclampvariablesampling(Fc,k,noiselevel,Fextmax,fr,ksf,Fe,kv,ev,mv,G,Fsample,tvec)
%
% This function simulates the hair-buyndle model from PNAS 2012.
%
% function [Xdet, Xsto, Deltadet, Deltasto,Fext] = hbtoymodelloadclampvariablesampling(Fc,k,noiselevel,Fextmax,fr,ksf,Fe,kv,ev,mv,G,Fsample,tvec)
%
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
% Fsample : sampling rate of clamp
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
mhb = 0.1;
ghb = 1;
if Fsample < 2
    Fsample = 2;
end
Dtsample = 1/Fsample;

%Decrease tvec step size by factor of Dtfac to ensure convergence
Dtfac = 1e2;
Dt = (tvec(2)-tvec(1))/Dtfac;
Dtsample = Dtsample/Dtfac;

N = round(tvec(end)/Dt);

% Add external forcing if desired
sinusoidalstim = 1;

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



%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments
fdW = sqrt(Dt)*randn(1,N); % White noise increments

xdet = zeros(1,N);
fdet = zeros(1,N);
xsto = zeros(1,N);
fsto = zeros(1,N);

xdet(1) = xzero;xdets(1)=xzero;
xsto(1) = xzero;xstos(1)=xzero;
vdet(1) = -xzero;
vsto(1) = -xzero;

fdet(1) = fzero;
fsto(1) = fzero;
fsfdet(1:5) = 1e-10;
fsfsto(1:5) = 1e-10;
xdetd(1:4) = 1e-10;xdetdd(1:4) = -1e-10;
xstod(1:4) = 1e-10;xstodd(1:4) = -1e-10;
xsfdet(1:4) = 1e-10;xsfsto(1:4) = 1e-10;

%Not using FD theorem
xNoiseSTD = noiselevel; fNoiseSTD = noiselevel; % equal noise levels

ns=1;
ks = 3;

for j = 2:N
    
if ks >= 3
if ns == (1/Fsample)*Dtfac
%Load Clamp - Deterministic
xdetd(j) = diff([xdets(ks-1) xdets(ks-1-1)]);
xdetdd(j) = diff([xdets(ks-1) xdets(ks-1-1) xdets(ks-1-2)],2);
xstod(j) = diff([xstos(ks-1) xstos(ks-1-1)]);
xstodd(j) = diff([xstos(ks-1) xstos(ks-1-1) xstos(ks-1-2)],2);
Xcdet(j) = (edx*Xcdet(ks-1))/((edx-ksf*Dt)*alpha) + Dt*((mv*xdetdd(j)+(ev-exx+alpha*beta*G*edx)*xdetd(j)+(kv-(1+alpha*beta*G)*ksf)*xdets(ks-1)-Fe)/(alpha*beta*G*(edx-ksf*Dt)) + real(Fext(j)));
Xcsto(j) = (edx*Xcsto(ks))/((edx-ksf*Dt)*alpha) + Dt*((mv*xstodd(j)+(ev-exx+alpha*beta*G*edx)*xstod(j)+(kv-(1+alpha*beta*G)*ksf)*xstos(ks-1)-Fe)/(alpha*beta*G*(edx-ksf*Dt)) + imag(Fext(j)));
Deltadet(j) = G*beta*alpha*(Xcdet(j)-xdet(j-1));
Deltasto(j) = G*beta*alpha*(Xcsto(j)-xsto(j-1));
%Vcdet(j) = (edx*Vcdet(j-1))/(edx-ksf*Dt) + Dt*((mv*xdetdd(j)+(ev-exx+G*edx)*xdetd(j)+(kv-(1+G)*ksf)*xdet(j)-Fe)/(G*(edx-ksf*Dt)));
%Vcsto(j) = (edx*Vcsto(j-1))/(edx-ksf*Dt) + Dt*((mv*xstodd(j)+(ev-exx+G*edx)*xstod(j)+(kv-(1+G)*ksf)*xsto(j)-Fe)/(G*(edx-ksf*Dt)));
%Deltadet(j) = (Vcdet(j)) - xdet(j-1);
%Deltasto(j) = (Vcsto(j)) - xdet(j-1);
fsfdet(j) = ksf*(Deltadet(j)-xdet(j-1));
xsfdet(j) = fsfdet(j)/k;
fsfsto(j) = ksf*(Deltasto(j)-xsto(j-1));
xsfsto(j) = fsfsto(j)/k;
ns=1;
if G == 0
    Deltadet(j) = 0;
    Deltasto(j) = 0;
end
else
    ns = ns+1;
    xdetd(j) = xdetd(j-1);
    xdetdd(j) = xdetdd(j-1);
    xstod(j) = xstod(j-1);
    xstodd(j) = xstodd(j-1);
    Xcdet(j) = Xcdet(j-1);
    Xcsto(j) = Xcsto(j-1);
    Deltadet(j) = Deltadet(j-1);
    Deltasto(j) = Deltasto(j-1);
    if G == 0
    Deltadet(j) = 0;
    Deltasto(j) = 0;
    end
    fsfdet(j) = fsfdet(j-1);
    fsfsto(j) = fsfsto(j-1);
    xsfdet(j) = xsfdet(j-1);
    xsfsto(j) = xsfsto(j-1);
    Deltadet(j) = 0;
    Deltasto(j) = 0;
    
end
end

%Deterministic integral
xdet(j) = xdet(j-1) + Dt*(vdet(j-1));
vdet(j) = vdet(j-1) + Dt*(-ghb/mhb*vdet(j-1)-k/mhb*xdet(j-1)+a/mhb*(xdet(j-1)-fdet(j-1))-1/mhb*(xdet(j-1)-fdet(j-1))^3 + Fc/mhb + Fext(j)/mhb + fsfdet(j)/mhb);
%vdet(j) = vdet(j-1) + Dt*(-ghb*vdet(j-1)-k*xdet(j-1)+a*(xdet(j-1)-fdet(j-1))-(xdet(j-1)-fdet(j-1))^3 + Fc + Fext(j) + fsfdet(j));
%xdet(j) = xdet(j-1) + Dt*(-k*xdet(j-1) + a*(xdet(j-1)-fdet(j-1)) - (xdet(j-1)-fdet(j-1))^3 + Fc + Fext(j) + fsfdet(j));
fdet(j) = fdet(j-1) + Dt*(b*xdet(j-1) - fdet(j-1))/tau;

%Stochastic integral
xsto(j) = xsto(j-1) + Dt*(vsto(j-1)) + xNoiseSTD*xdW(j);
vsto(j) = vsto(j-1) + Dt*(-ghb/mhb*vsto(j-1)-k/mhb*xsto(j-1)+a/mhb*(xsto(j-1)-fsto(j-1))-1/mhb*(xsto(j-1)-fsto(j-1))^3 + Fc/mhb + Fext(j)/mhb + fsfsto(j)/mhb);
%vsto(j) = vsto(j-1) + Dt*(-ghb*vsto(j-1)-k*xsto(j-1)+a*(xsto(j-1)-fsto(j-1))-(xsto(j-1)-fsto(j-1))^3 + Fc + Fext(j) + fsfsto(j));
%xsto(j) = xsto(j-1) + Dt*(-k*xsto(j-1) + a*(xsto(j-1)-fsto(j-1)) - (xsto(j-1)-fsto(j-1))^3 + Fc + Fext(j) + fsfsto(j)) + xNoiseSTD*xdW(j);
fsto(j) = fsto(j-1) + Dt*(b*xsto(j-1) - fsto(j-1))/tau + fNoiseSTD*fdW(j)/tau;

if j == 1
    xdets(2) = xdet(j);xstos(2)=xsto(j);
end
if ns == (1/Fsample)*Dtfac
    xdets(ks) = xdet(j);xstos(ks)=xsto(j);
    ks=ks+1;
end

end

Xdet = zeros(2,length(tvec)-1);
Xsto = zeros(2,length(tvec)-1);

%Return vectors at tvecs specified by tvec.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = fdet(1:Dtfac:N);
Xsto(1,:) = xsto(1:Dtfac:N);
Xsto(2,:) = fsto(1:Dtfac:N);

plotyn = 1;
if plotyn == 1
figure
plot(Ftime(1:end),xsto,'r');
hold on
plot(Ftime(1:end),xdet,'k');
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
