% Simple routine to plot r as a function of x_star

% Range of x_star (here x_star = x)
     x = -4:0.1:4;

% Auxiliary variable for plotting y-axis
     y = -4:0.1:4;

% Working out what r is
     r = x + exp(-x);

% Auxiliary variable for plotting x- and y-axis
     ay = 0.0.*x
     ax = 0.0.*y

% Plotting bifurcation curve
     plot(x,r,'r')

% Determining range for plot
     axis([-4 4 -4 4])
     hold on

% Plotting x- and y- axes
     plot(x,ay,'--')
     plot(ax,y,'--')

% Labelling
     xlabel('x_*','Fontsize',20)
     ylabel('r','FontSize',20)
     text(-3,-2,'r = x_* + e^{-x_*}','Color','r','FontSize',20)
