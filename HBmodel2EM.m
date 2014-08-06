function [Xdet Xem] = HBmodel2EM( Time )
%Stochasic HB model integration

%EM Euler-Maruyama method
%Ito integral

a = 3.5;
%b > 1 has unbounded solutions
b = 0.5;
tau = 60;
Fc = 0;
gamma = 0.5;
k = 2.5;
xzero = 1; 
fzero = 1;

%Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 100;
Dt = (Time(2)-Time(1))/Dtfac;

N = numel(Time)*Dtfac;

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
xNoiseSTD = 0.1;
fNoiseSTD = 0.1;

for j = 2:N
%Deterministic integral
xdet(j) = xdet(j-1) + Dt*((-k*xdet(j-1) + a*(xdet(j-1)-fdet(j-1)) - (xdet(j-1)-fdet(j-1))^3 + Fc)/gamma);
fdet(j) = fdet(j-1) + Dt*((b*xdet(j-1) - fdet(j-1))/tau);

%Stochastic integral
xem(j) = xem(j-1) + Dt*((-k*xem(j-1) + a*(xem(j-1)-fem(j-1)) - (xem(j-1)-fem(j-1))^3 + Fc)/gamma) + xNoiseSTD*xdW(j);
fem(j) = fem(j-1) + Dt*((b*xem(j-1) - fem(j-1))/tau) + fNoiseSTD*fdW(j);
end

Xdet = zeros(2,N/Dtfac);
Xem = zeros(2,N/Dtfac);

%Return vectors at times specified by Time.
Xdet(1,:) = xdet(find(mod(1:N,Dtfac) == 0));
Xdet(2,:) = fdet(find(mod(1:N,Dtfac) == 0));
Xem(1,:) = xem(find(mod(1:N,Dtfac) == 0));
Xem(2,:) = fem(find(mod(1:N,Dtfac) == 0));

plotcheck = 0;
if plotcheck == 1
figure
plot(Time,Xem(1,:),'r');
hold on
plot(Time,Xdet(1,:),'k');
xlabel('Time','FontSize',24) 
ylabel('x','FontSize',24,'Rotation',0,'HorizontalAlignment','right')

figure
plot(Time,Xem(2,:),'g');
hold on
plot(Time,Xdet(2,:),'k');
xlabel('Time','FontSize',24) 
ylabel('f','FontSize',24,'Rotation',0,'HorizontalAlignment','right')
end

end
