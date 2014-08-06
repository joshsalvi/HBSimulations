function [t,z,z_dot]=hopf2(mu,mu_c,w0,z0,t_length,a)

% [t,z,z_dot]=hopf2(mu,mu_c,w0,z0,t_length,a)

t = 0:0.1:t_length;
z(1)=z0;
z_dot(1) = ((mu-mu_c)+i*w0)*z(1)+a*z(1).^3;

for j = 2:length(t)
    z(j) = z(j-1)+z_dot(j-1)*0.1;
    z_dot(j) = ((mu-mu_c)+i*w0)*z(j)+a*z(j)*abs(z(j)).^2;
end
