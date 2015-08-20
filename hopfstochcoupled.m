function [X1det, X1sto, X2det, X2sto, Fext] = hopfstochcoupled(mu1,mu2,fosc1,fosc2,K,noiselevel,Fextmax,fr,tvec)
%
% This function simulates the normal form of the supercritical Hopf
% bifurcation, given by two planar equations:
%
% x_dot = mu*x - omega*y - x*(x^2 + y^2)
% y_dot = omega*x + mu*y - y*(x^2 + y^2)
%
% where mu is the control parameter. For mu>0, the system will oscillate at
% an amplitude that grows with sqrt(mu).  Alternatively, one may express the
% above equations in polar coordinates, making the amplitude relationship
% with respect to mu more apparent:
%
% rho_dot = rho*(mu + i*omega - rho^2)
%
% Here we simulate both the deterministic and stochastic cases for the
% supercritical Hopf bifurcation for two coupled oscillators.
%
% [X1det, X1sto, X2det, X2sto, Fext] = hopfstochcoupled(mu1,mu2,fosc1,fosc2,K,noiselevel,tvec)
%
% Xdet : deterministic result
% Xsto : stochastic result
%
%  
% tvec : time vector
% mu : control parameter
% noiselevel : standard deviation of stochastic noise
% K : coupling constant
% fosc : frequency of oscillation on the unstable side of the bifurcation
%
% By modifying the code, you can also add a step function or external
% forcing. 
%
% jsalvi@rockefeller.edu
%

% All noise levels equal
x1NoiseSTD=noiselevel;
y1NoiseSTD=noiselevel;
x2NoiseSTD=noiselevel;
y2NoiseSTD=noiselevel;
% Initial condition
x1zero = 0.1;x2zero = -0.1;
y1zero = -0.1;y2zero = 0.1;

% Add external forcing if desired
sinusoidalstim = 1; pulsestim = 0;  % pulse or sinusoid?
%Fextmax = 1;        % amplitude of sinusoidal stim OR pulse
%fr = 5;             % frequency of stimulation
pulsestart = 1;     % start of pulse
pulseend = 2;       % end of pulse

% Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^2;
Dt = (tvec(2)-tvec(1))/Dtfac;

N = tvec(end)/Dt;

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
x1dW = sqrt(Dt)*randn(1,N); % White noise increments
y1dW = sqrt(Dt)*randn(1,N); % White noise increments
x2dW = sqrt(Dt)*randn(1,N); % White noise increments
y2dW = sqrt(Dt)*randn(1,N); % White noise increments


x1det = zeros(1,N); x1det(1) = x1zero;
x1sto = zeros(1,N); x1sto(1) = x1zero;
y1det = zeros(1,N); y1det(1) = y1zero;
y1sto = zeros(1,N); y1sto(1) = y1zero;
x2det = zeros(1,N); x2det(1) = x2zero;
x2sto = zeros(1,N); x2sto(1) = x2zero;
y2det = zeros(1,N); y2det(1) = y2zero;
y2sto = zeros(1,N); y2sto(1) = y2zero;

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
x1det(j) = x1det(j-1) + Dt*(mu1*x1det(j-1) - 2*pi*fosc1*y1det(j-1) - x1det(j-1)*(x1det(j-1)^2 + y1det(j-1)^2) + K*(x2det(j-1)-x1det(j-1)) + real(Fext(j)));
y1det(j) = y1det(j-1) + Dt*(2*pi*fosc1*x1det(j-1) + mu1*y1det(j-1) - y1det(j-1)*(x1det(j-1)^2 + y1det(j-1)^2) + K*(y2det(j-1)-y1det(j-1)) + imag(Fext(j)));
x2det(j) = x2det(j-1) + Dt*(mu2*x2det(j-1) - 2*pi*fosc2*y2det(j-1) - x2det(j-1)*(x2det(j-1)^2 + y2det(j-1)^2) + K*(x1det(j-1)-x2det(j-1)) + real(Fext(j)));
y2det(j) = y2det(j-1) + Dt*(2*pi*fosc2*x2det(j-1) + mu2*y2det(j-1) - y2det(j-1)*(x2det(j-1)^2 + y2det(j-1)^2) + K*(y1det(j-1)-y2det(j-1)) + imag(Fext(j)));

%Stochastic integral
x1sto(j) = x1sto(j-1) + Dt*(mu1*x1sto(j-1) - 2*pi*fosc1*y1sto(j-1) - x1sto(j-1)*(x1sto(j-1)^2 + y1sto(j-1)^2) + K*(x2det(j-1)-x1sto(j-1)) + real(Fext(j))) + x1NoiseSTD*x1dW(j);
y1sto(j) = y1sto(j-1) + Dt*(2*pi*fosc1*x1sto(j-1) + mu1*y1sto(j-1) - y1sto(j-1)*(x1sto(j-1)^2 + y1sto(j-1)^2) + K*(y2sto(j-1)-y1sto(j-1)) + imag(Fext(j))) + y1NoiseSTD*y1dW(j);
x2sto(j) = x2sto(j-1) + Dt*(mu2*x2sto(j-1) - 2*pi*fosc2*y2sto(j-1) - x2sto(j-1)*(x2sto(j-1)^2 + y2sto(j-1)^2) + K*(x1sto(j-1)-x2sto(j-1)) + real(Fext(j))) + x2NoiseSTD*x2dW(j);
y2sto(j) = y2sto(j-1) + Dt*(2*pi*fosc2*x2sto(j-1) + mu2*y2sto(j-1) - y2sto(j-1)*(x2sto(j-1)^2 + y2sto(j-1)^2) + K*(y1sto(j-1)-y2sto(j-1)) + imag(Fext(j))) + y2NoiseSTD*y2dW(j);

end


X1det = zeros(2,length(tvec)-1);X2det = zeros(2,length(tvec)-1);
X1sto = zeros(2,length(tvec)-1);X2sto = zeros(2,length(tvec)-1);

%Return vectors at times specified by Time.
X1det(1,:) = x1det(1:Dtfac:N);
X1det(2,:) = y1det(1:Dtfac:N);
X1sto(1,:) = x1sto(1:Dtfac:N);
X1sto(2,:) = y1sto(1:Dtfac:N);
X2det(1,:) = x2det(1:Dtfac:N);
X2det(2,:) = y2det(1:Dtfac:N);
X2sto(1,:) = x2sto(1:Dtfac:N);
X2sto(2,:) = y2sto(1:Dtfac:N);

% Make a plot of the data?
plotyn=0;

if plotyn==1
    figure;
    subplot(1,2,1);hold on;plot(tvec(2:end),X1sto(1,:),'r');plot(tvec(2:end),X1det(1,:),'k');title('Black=deterministic; Red=stochastic; real part only');
    subplot(1,2,2);hold on;plot(X1sto(1,:),X1sto(2,:),'r');plot(X1det(1,:),X1det(2,:),'k');title('Black=deterministic; Red=stochastic');
end


end
