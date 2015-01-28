function [Xdet] = vilfanduke(F,k,tvec)
%
% This function simulates the Ca2+ dynamics model from Vilfan/Duke (2003).
%
% where mu is the control parameter. For mu>0, the system will oscillate a
%
% Here we simulate both the deterministic and stochastic cases for the
% supercritical Hopf bifurcation 
%
% [Xdet] = vilfanduke(F,k,tvec)
%
% jsalvi@rockefeller.edu
%

% Initial conditions and time vector
xzero = 1e-9;p0zero=0.2;p1zero=0.15;p2zero=0.5;czero=1e-6;zzero=1;
ghb = 1;kgs=816e-6;
d1=0;d2=0;d0=98e-9;
C0=3.6e-6;Cb=0.1e-6;Cm=7e-6;
G0=0;G1=0;G2=70e-21;
T=295;kb=4.1e-21;
lambda = 6.4e-3;D=28e-9;
N=50;xm=9.5e-9;

% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^1;
Dt = (tvec(2)-tvec(1))/Dtfac;
tvec2 = interp(tvec,round(1/Dt));
N = tvec(end)/Dt;

xdet = zeros(1,N); xdet(1) = xzero;
p0det = zeros(1,N); p0det(1) = p0zero;
p1det = zeros(1,N); p1det(1) = p1zero;
p2det = zeros(1,N); p2det(1) = p2zero;
cdet = zeros(1,N); cdet(1) = czero;
zdet = zeros(1,N); zdet(1) = zzero;

% oscillate calcium
Ftime = linspace(tvec(1),tvec(end),N);
cdet = C0*cos(2*pi*1*Ftime);

% Euler-Murayama Method with Ito Integration
for j = 2:N

    %xdet(j) = xdet(j-1) + Dt*( -k/ghb*(xdet(j-1) - xm) - kgs/ghb*(xdet(j-1) - d0*p0det(j-1) - d1*p1det(j-1) - d2*p2det(j-1)) + F/ghb );
    xdet(j) = xdet(j-1) + Dt*( -k/ghb*(xdet(j-1) - xm) - kgs/ghb*(xdet(j-1) - d0*p0det(j-1)) + F/ghb );
    %cdet(j) = cdet(j-1) + Dt*( exp(sqrt(-1)*cdet(j-1)*tvec2(j-1)) );
    %cdet(j) = cdet(j-1) + Dt*( -lambda*(cdet(j-1) - Cb - Cm*p0det(j-1)) );
    p0det(j) =  1/zdet(j-1)*exp( (-kgs*(xdet(j-1)-d0-D*log(cdet(j-1)/C0))^2 - 2*G0) / (2*N*kb*T) ) ;
    p1det(j) =  1/zdet(j-1)*exp( (-kgs*(xdet(j-1)-d1-D*log(cdet(j-1)/C0))^2 - 2*G1) / (2*N*kb*T) ) ;
    p2det(j) =  1/zdet(j-1)*exp( (-kgs*(xdet(j-1)-d2-D*log(cdet(j-1)/C0))^2 - 2*G2) / (2*N*kb*T) ) ;
    zdet(j) =  exp( (-kgs*(xdet(j-1)-d0-D*log(cdet(j-1)/C0))^2 - 2*G0) / (2*N*kb*T) ) + exp( (-kgs*(xdet(j-1)-d1-D*log(cdet(j-1)/C0))^2 - 2*G1) / (2*N*kb*T) ) + exp( (-kgs*(xdet(j-1)-d2-D*log(cdet(j-1)/C0))^2 - 2*G2) / (2*N*kb*T) ) ;

end



%Return vectors at times specified by Time.
Xdet(1,:) = xdet(1:Dtfac:end);
Xdet(2,:) = cdet(1:Dtfac:end);
Xdet(3,:) = p0det(1:Dtfac:end);
Xdet(4,:) = p0det(1:Dtfac:end);
Xdet(5,:) = p0det(1:Dtfac:end);

% Make a plot of the data?
plotyn=1;
if plotyn==1
    figure;
    subplot(1,2,1);plot(Xdet(1,:),'k');title('Bundle Motion');
    subplot(1,2,2);plot(Xdet(2,:),'k');title('Calcium Dynamics');
end


end
