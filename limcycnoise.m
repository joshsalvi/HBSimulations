% Limit Cycle Oscillations; all are in arbitrary units


%standard deviation of noise
sigma=0.5;

%Range for oscillation amplitudes
% Easiest to visualize if low=high
low=1;
high=1;
N=high/low;
for i=1:N
    amp(i)=low*i;
end

% Create wavelet for relaxation oscillator, initial amplitude is ±2
for i=1:N
    [t,y]=call_osc();
    t_tot(:,:,i)=t;
    y2=amp(i).*y;
    y_amp(:,:,i)=y2;
%    figure(1); hold on;
%    plot(t_tot(:,:,i),y_amp(:,1,i)); hold off;
    figure(3); subplot(2,1,1); hold on;
    plot(fliplr(t_tot(:,:,i)'),fliplr(y_amp(:,1,i)')); xlabel('Time (a.u.)');ylabel('Amplitude (a.u.)');
end
title('Noise-Free Simulation [TOP], Simulation Plus Noise [BOTTOM]');
hold off;

for i=1:length(y_amp(:,1,1))
    for j=1:length(y_amp(1,1,:))
        y_amp(i,1,j)=y_amp(i,1,j)+sigma*(3.92*random('Normal',0,0.2)); % 95-percent confidence interval
    end
end

for i=1:N
%    figure(2); hold on;
%    plot(t_tot(:,:,i),y_amp(:,1,i)); hold off;
    figure(3); subplot(2,1,2); hold on;
    plot(fliplr(t_tot(:,:,i)'),fliplr(y_amp(:,1,i)')); hold off;
    xlabel('Time (a.u.)');ylabel('Amplitude (a.u.)');
end
title({sprintf('Standard Deviation of Noise = %0.5g',sigma)});
hold off;
