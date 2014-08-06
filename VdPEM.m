function [Xdet Xem] = VdPEM( Time )
%Stochasic HB model integration

%EM Euler-Maruyama method
%Ito integral

%Period is greater than 2*pi/sqrt(k/m)
m = 2;
%Increasing gamma increases nonlinearity
gamma = -5;
k = 0.1;
Fc = 0;
xzero = 1; 
vzero = 1;

%Decrease time step size by factor of Dtfac to ensure convergence
Dtfac = 10000;
Dt = (Time(2)-Time(1))/Dtfac;

N = numel(Time)*Dtfac;

%Set the default random number stream
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1))
xdW = sqrt(Dt)*randn(1,N); % White noise increments

xdet = zeros(1,N);
vdet = zeros(1,N);
xem = zeros(1,N);
vem = zeros(1,N);

xdet(1) = xzero;
xem(1) = xzero;

vdet(1) = vzero;
vem(1) = vzero;

%Not using FD theorem
xNoiseSTD = 0.1;

for j = 2:N
%Deterministic integral
xdet(j) = xdet(j-1) + Dt*vdet(j-1);
vdet(j) = vdet(j-1) + Dt*((-gamma*(1-xdet(j-1)^2)*vdet(j-1) - k*xdet(j-1) + Fc)/m);

%Stochastic integral
xem(j) = xem(j-1) + Dt*vem(j-1) + xNoiseSTD*xdW(j);
vem(j) = vem(j-1) + Dt*((-gamma*(1-xem(j-1)^2)*vem(j-1) - k*xem(j-1) + Fc)/m);
end

Xdet = zeros(2,N/Dtfac);
Xem = zeros(2,N/Dtfac);

%Return vectors at times specified by Time.
Xdet(1,:) = xdet(find(mod(1:N,Dtfac) == 0));
Xdet(2,:) = vdet(find(mod(1:N,Dtfac) == 0));
Xem(1,:) = xem(find(mod(1:N,Dtfac) == 0));
Xem(2,:) = vem(find(mod(1:N,Dtfac) == 0));

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
