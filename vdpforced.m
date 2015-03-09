function [Xout, Xin3] = vdpforced(mu,fosc,Xin,Fs)
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
% [Xout, Xin] = hopfstoforced(mu,fosc,Xin,Fs)
%
% Xout : signal filtered by a Hopf oscillator
% fosc : frequency of Hopf oscillator
% Xin : original signal
% Fs : sample rate (Hz)
%
%
% jsalvi@rockefeller.edu
%

% Initial conditions and time vector
xzero = 0;yzero=0;
tvec = linspace(0,length(Xin)/Fs-1/Fs,length(Xin));

% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 30^2;
Dt = (tvec(2)-tvec(1))/Dtfac;

Xin2 = interp(Xin,round(1/Dt));
tvec = interp(tvec,round(1/Dt));
N = length(Xin2);
%N = tvec(end)/Dt;

xdet = zeros(1,N); xdet(1) = xzero;
ydet = zeros(1,N); ydet(1) = yzero;

% Scale the signal
xmax = max(abs(Xin2));
if xmax ~=0
    Xin2 = 1e-8.*(Xin2./xmax);
end

% Euler-Murayama Method with Ito Integration
for j = 2:N
    
    xdet(j) = xdet(j-1) + Dt*(mu*xdet(j-1) - mu*1/3*xdet(j-1)^3 - mu*ydet(j-1) - (Xin2(j)));
    ydet(j) = ydet(j-1) + Dt*(1/mu*xdet(j-1) - (Xin2(j)));
    
end



%Return vectors at times specified by Time.
Xout(1,:) = xdet(1:Dtfac:end);
Xout(2,:) = ydet(1:Dtfac:end);
Xin3 = Xin2(1:Dtfac:end);

% Make a plot of the data?
plotyn=1;
if plotyn==1
    figure;
    subplot(1,2,1);plot(Xin3,'k');title('Input Signal');
    subplot(1,2,2);plot(Xout(1,:),'k');title('Output Signal');
end


end
