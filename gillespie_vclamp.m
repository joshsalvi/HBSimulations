% gillespie_vclamp.m
% runs the gillespie simultation for a two-state channel under voltage
% clamp
% 
function [output] = gillespie_vclamp

    NV = [1]; % number of channels

    ts = 10;  % time to switch V, ms
    tend = ts;    % time to end simulation, ms

    figure
    hold
    
    for i = 1:length(NV)
        [t, g] = gillespie_function(NV(i));
        output{i,1} = [t' g'];
        stairs(t,g/NV(i)+(length(NV)+1-i)*1.2) % plot normalized, offset conductance
        % axis([0 tend*1.05 -0.1 1.1]) % normalize plot axes
    end
    
    % plot deterministic solution
    tdet = [0:0.05:20];
    vi = -80;   % initial voltage, mV
    vf = -20;   % final voltage, mV
    [a0 b0 m0 ta0] = mparams(vi); % initial values of alpha, beta, minf and taum
    [a1 b1 m1 ta1] = mparams(vf);  % final values of alpha, beta, minf, and taum

    gdet = m0 + (tdet>ts).*(m1-m0).*(1-exp(-(tdet-ts)/ta1));
    output{length(NV)+1,1} = [tdet' gdet'];
    plot(tdet,gdet,'k')
    
    tscale = [15 20];
    gscale = [-0.1 -0.1];
    plot(tscale,gscale,'k')
    
    gdivide = [-1 6];
    plot([ts ts],gdivide,':')

end % main function


function [tv,N_open] = gillespie_function(N)
    gamma = 10; % open-channel conductance, pS
    %N  = 100;  % number of channels
    vi = -80;   % initial voltage, mV
    vf = -20;   % final voltage, mV
    ts = 10;  % time to switch V, ms
    tend = 222*ts;    % time to end simulation, ms
    [a0 b0 m0 ta0] = mparams(vi); % initial values of alpha, beta, minf and taum
    [a1 b1 m1 ta1] = mparams(vf);  % final values of alpha, beta, minf, and taum

    % seed random number generator
    rand('twister',sum(1000*clock));

    No = round(m0*N);   % initial number of open channels
    Nc = N-No;          % initial number of closed channels
	t = 0;  % initial time
	i = 1;  % counter for vector of transition times
	% simulate first portion
	tv(1) = 0;
	N_open(1)=No;
	while (t<tend)
        i = i + 1;
        if (t<ts)
            lm(i) = Nc*a0 + No*b0;
            beta = b0;
        else
            lm(i) = Nc*a1 + No*b1;
            beta = b1;
        end
        r = rand(2,1);
        dt(i) = 1/lm(i)*log(1/r(1));
        if (No*beta > r(2)*lm(i))
            No = No - 1;
            Nc = Nc + 1;
        else
            No = No + 1;
            Nc = Nc - 1;
        end
        N_open(i)=No;
        N_closed(i)=Nc;
        t = t+dt(i);
        tv(i) = t;
    end % while loop
    
end % gillespie_function

    
% mparams(v)
% returns taum and minf as functions of v
% HH (6 deg C) parameters

function [a,b,minf,taum] = mparams(v)
    a = alm(v);
    b = bem(v);
    taum = 1/(alm(v)+bem(v));
    minf = alm(v)*taum;
end

function alm = alm(v)
    x = -(v+40);
    y = 10;
    if abs(x/y) < 1.e-6
        alm = 0.1*y*(1 - x/y/2);    % this "trap" catches divide by zero errors
    else
        alm = 0.1*x/(exp(x/y) - 1);
    end % if/then function
end % function call

function bem = bem(v)
bem = 4*exp(-(v+65)/18);
end