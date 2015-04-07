function [Sig DetSig Time Fext DetSigf Sigf] = ONHmodel(Fextmax,fr)

global FEXTMAX;
global FR;
global TIME;

FEXTMAX = Fextmax;
FR = fr;

t_range = 1000;
%100 sample points per period
Time = 0:1:t_range;

TIME = Time;

lin = 1;
if lin == 1
%Sinusoidal function
TestT = 10;%period length for sinusoidal function
w = 2*pi/TestT;
DetSig = cos(w.*Time);
Sig = DetSig + 0.1.*randn(length(Time),1)';
else

%%%%%%%Simulations%%%%%%%%%

%TimeS = 0:1:2*t_range;
%Time = TimeS(find(TimeS == t_range):find(TimeS == 2*t_range));

%Stochastic Ito integration
%[Xdet Xem] = HBmodelEM(TimeS,Fextmax,fr);
[Xdet Xem] = HBmodelEM(Time,Fextmax,fr);

%%%%%Use Matlab's Integrator%%%%%%
%options = odeset('MaxStep', 10^-1);
%sol = ode45(@system,[0 t_range],[1 0],options);
%[Y,YP] = deval(sol,Time);`     `` \
 

xlabel('Time','FontSize',24); ylabel('x','FontSize',24);
%}
%savefile = '/Users/dmelody/Work/Ear/Hair Bundle Expts/Data/April13/Sigs';
%save([savefile,'.mat'], 'Sig', 'DetSig', 'Time');




end
function dx = system(t,x)


global FEXTMAX;
global FR;
global TIME;

 dx = zeros(2,1);
 Fext = FEXTMAX*cos(2*pi*FR*TIME);
 
 a = 500;
%b > 1 has unbounded solutions
b = 0.96;
tau = 100;
Fshift = 2000;
Fc = 0 + Fshift;
ke = 100; 
gamma = 5*10^2;
ksp = 00;
k = ksp + ke;
 
 Fext = interp1(TIME,Fext,t);
 dx(1) = (-k*x(1) + a*(x(1)-x(2)) - (x(1)-x(2))^3 + Fc + Fext)/gamma;
 dx(2) = (b*x(1) - x(2))/tau;
 
end
end
