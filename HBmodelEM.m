function [Xdet Xem] = HBmodelEM(Time,Fextmax,fr)
%Stochasic HB model integration

%EM Euler-Maruyama method
%Ito integral

a = 3.5;
%b > 1 has unbounded solutions
b = 0.5;
tau = 10;
Fc = 0;
k = 3;

xzero = 1; 
fzero = 0;

%Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10^2;
Dt = (Time(2)-Time(1))/Dtfac;

%N = numel(Time)*Dtfac;
N = 1+Time(end)/Dt;

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments
fdW = sqrt(Dt)*randn(1,N); % White noise increments

xdet = zeros(1,N);
fdet = zeros(1,N);
xem = zeros(1,N);
fem = zeros(1,N);

xdet(1) = xzero;
xem(1) = xzero;

fdet(1) = fzero;
fem(1) = fzero;

%Not using FD theorem
xNoiseSTD = 0;
fNoiseSTD = 0;

Ftime = linspace(Time(1),Time(end),N);
Fext = Fextmax*cos(2*pi*fr*Ftime);

for j = 2:N
%Deterministic integral
xdet(j) = xdet(j-1) + Dt*(-k*xdet(j-1) + a*(xdet(j-1)-fdet(j-1)) - (xdet(j-1)-fdet(j-1))^3 + Fc + Fext(j));
fdet(j) = fdet(j-1) + Dt*(b*xdet(j-1) - fdet(j-1))/tau;

%Stochastic integral
xem(j) = xem(j-1) + Dt*(-k*xem(j-1) + a*(xem(j-1)-fem(j-1)) - (xem(j-1)-fem(j-1))^3 + Fc + Fext(j)) + xNoiseSTD*xdW(j);
fem(j) = fem(j-1) + Dt*(b*xem(j-1) - fem(j-1))/tau + fNoiseSTD*fdW(j)/tau;
end

Xdet = zeros(2,length(Time));
Xem = zeros(2,length(Time));

%Return vectors at times specified by Time.
Xdet(1,:) = xdet(1:Dtfac:N);
Xdet(2,:) = fdet(1:Dtfac:N);
Xem(1,:) = xem(1:Dtfac:N);
Xem(2,:) = fem(1:Dtfac:N);

plotcheck = 1;
if plotcheck == 1
figure
plot(Ftime,xem,'r');
hold on
plot(Ftime,xdet,'k');
xlabel('Time','FontSize',24) 
ylabel('x','FontSize',24,'Rotation',0,'HorizontalAlignment','right')

figure
plot(Ftime,fem,'g');
hold on
plot(Ftime,fdet,'k');
xlabel('Time','FontSize',24) 
ylabel('f','FontSize',24,'Rotation',0,'HorizontalAlignment','right')
end

end
