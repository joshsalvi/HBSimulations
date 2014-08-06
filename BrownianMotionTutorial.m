% Brownian Motion Tutorial.

% 

N = 1000; displacement = randn(1,N); plot(displacement);

figure(1);
hist(displacement, 25);


figure(2);
x = cumsum(displacement);
plot(x);
ylabel('position');
xlabel('time step');
title('Position of 1D Particle versus Time');

figure(3);
particle = struct();
particle.x = cumsum( randn(N, 1) );
particle.y = cumsum( randn(N, 1) );
plot(particle.x, particle.y);
ylabel('Y Position');
xlabel('X Position');
title('position versus time in 2D');

figure(4);
dsquared = particle.x .^ 2 + particle.y .^ 2;
plot(dsquared);
title('displacement squared');


%% Theoretical value of diffusion coefficient
close all;

d    = 1.0e-6;              % radius in meters
eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
kB   = 1.38e-23;            % Boltzmann constant
T    = 293;                 % Temperature in degrees Kelvin

D    = kB * T / (3 * pi * eta * d)


%% More realistic particle
close all;

dimensions = 2;         % two dimensional simulation
tau = .1;               % time interval in seconds
time = tau * 1:N;       % create a time vector for plotting

k = sqrt(D * dimensions * tau);
dx = k * randn(N,1);
dy = k * randn(N,1);

x = cumsum(dx);
y = cumsum(dy);

dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
 squaredDisplacement = ( x .^ 2) + ( y .^ 2);

 figure;
plot(x,y);
title('Particle Track of a Single Simulated Particle');

figure;

hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3);      % plot theoretical line

plot(time, squaredDisplacement);
hold off;
xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time for 1 Particle in 2 Dimensions');

simulatedD = mean( dSquaredDisplacement ) / ( 2 * dimensions * tau )
standardError = std( dSquaredDisplacement ) / ( 2 * dimensions * tau * sqrt(N) )
actualError = D - simulatedD


figure;
dx = dx + 0.2 * k;
dy = dy + 0.05 * k;

x = cumsum(dx);
y = cumsum(dy);

dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
 squaredDisplacement = ( x .^ 2) + ( y .^ 2);

simulatedD    = mean( dSquaredDisplacement ) / ( 2 * dimensions * tau )
standardError = std(  dSquaredDisplacement ) / ( 2 * dimensions * tau * sqrt(N) )
actualError = D - simulatedD

plot(x,y);
title('Particle Track of a Single Simulated Particle with Bulk Flow');

%%

clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3);      % plot theoretical line
plot(time, squaredDisplacement);
hold off;

xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time with Bulk Flow');

particleCount = 10;
N = 50;
tau = .1;
time = 0:tau:(N-1) * tau;
particle = { };             % create an empty cell array to hold the results

for i = 1:particleCount
    particle{i} = struct();
    particle{i}.dx = k * randn(1,N);
    particle{i}.x = cumsum(particle{i}.dx);
    particle{i}.dy = k * randn(1,N);
    particle{i}.y = cumsum(particle{i}.dy);
    particle{i}.drsquared = particle{i}.dx .^2 + particle{i}.dy .^ 2;
    particle{i}.rsquared = particle{i}.x .^ 2 + particle{i}.y .^ 2;
    particle{i}.D = mean( particle{i}.drsquared ) / ( 2 * dimensions * tau );
    particle{i}.standardError = std( particle{i}.drsquared ) / ( 2 * dimensions * tau * sqrt(N) );
end

figure;
clf;
hold on;
for i = 1:particleCount
    plot(particle{i}.x, particle{i}.y, 'color', rand(1,3));
end

xlabel('X position (m)');
ylabel('Y position (m)');
title('Combined Particle Tracks');
hold off;

%%
% compute the ensemble average
rsquaredSum = zeros(1,N);

for i = 1:particleCount
    rsquaredSum = rsquaredSum + particle{i}.rsquared;
end

ensembleAverage = rsquaredSum / particleCount;

% create the plot
clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'b', 'LineWidth', 3);      % plot theoretical line

plot(time, ensembleAverage , 'k', 'LineWidth', 3);          % plot ensemble average
legend('Theoretical','Average','location','NorthWest');

for i = 1:particleCount
    plot(time, particle{i}.rsquared, 'color', rand(1,3));   % plot each particle track
end

xlabel('Time (seconds)');
ylabel('Displacement Squared (m^2)');
title('Displacement Squared vs Time');
hold off;

%%
clear D e dx;

% extract the D value from each simulation and place them all into a single
% matrix called 'D'
for i = 1:particleCount
    D(i) = particle{i}.D;
    dx(i,:) = particle{i}.dx;
    e(i) = particle{i}.standardError;
end

% compute the estimate of D and the uncertainty
averageD = mean(D)
uncertainty = std(D)/sqrt(particleCount)

% plot everything
clf;
hold on;

plot(averageD * ones(1,particleCount), 'b', 'linewidth', 3);                    % plot estimated D
plot((averageD + uncertainty) * ones(1,particleCount), 'g-', 'linewidth', 1);   % plot upper error bar
plot((averageD - uncertainty) * ones(1,particleCount), 'g-', 'linewidth', 1);   % plot lower error bar
errorbar(D,e,'ro');                                                             % plot D values with error bars

xlabel('Simulation Number');
ylabel('Estimated Diffusion Coefficient');
title('Estimated Diffusion Coefficient with Error Bars')
legend('Average Value of D', 'location', 'NorthWest');

hold off;
