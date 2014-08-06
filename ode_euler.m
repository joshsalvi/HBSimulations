% ode_euler(x,f_o)
% There are two input arguments, x and f_o
% x defines a VECTOR with time points for the function f
% f_o is the initial condition
% Function returns VECTOR f representing approximate solution to
% differential equation, df/dx=2x with f(0)=f_o


function f = ode_euler(x,f_o)
% There are two input arguments, x and f_o
% x defines a VECTOR with time points for the function f
% f_o is the initial condition
% Function returns VECTOR f representing approximate solution to
% differential equation, df/dx=2x with f(0)=f_o

% Delta_x is difference between successive x values
delta_x = x(2) - x(1);

% Number of points necessary to approximate x
l_x = length(x);

% Initialize f 
f= zeros(1,l_x);
f(1) = f_o;

for i = 1:(l_x)
    f(i+1) = f(i) + delta_x*2*x*(i);
end
