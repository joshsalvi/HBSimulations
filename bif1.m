%% Animation for a Saddle Node bifurcation

% Range of x
     xmin = -4;
     xmax = 4;
     x = xmin:0.1:xmax
% Some auxiliary variable to plot x-axis
     ax = (0.0).*x;

% Range of r
     rmin = -4
     rmax = 4

% number of frames in the movie
     imax = 21

% Rage of dx/dt (here y=dx/dt)
     ymin = rmin ;
     ymax = rmax + xmax.^2;

% Some auxiliary variables to plot the y-axis
     ky = ymin:0.1:ymax
     ay = (0.0).*ky

%% Start of the loop for acquisition of the movie frames

for i = 1:imax

%calculate r value
     r = rmin + (i-1)*(rmax-rmin)/(imax-1)

%create string for labelling r
     strr = num2str(r);
     strt = ['r=' strr]; 

%calculate dxdt
     y = r + x.^2

%plot x-axis
     plot(x,ax,'--');

% determine range of plot
     axis([xmin xmax ymin ymax])
     hold on

%plot y-axis
     plot(ay,ky,'--')

%plot f(x)
     plot(x,y,'r');

%plot fixed points
if (r < 0) 
     plot(-(-r).^(0.5),0,'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
     plot((-r).^(0.5),0,'o','MarkerSize',10,'MarkerEdgeColor','k');
end
if(r == 0.0) 
     plot((-r).^(0.5),0,'*','MarkerSize',10,'MarkerEdgeColor','k');
end

%plot flow direction
if(r < 0)
     plot(-3,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
     plot(0,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(3,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end

if(r == 0) 
     plot(-3,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');     
plot(3,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end

if(r > 0) 
     plot(-3,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');     plot(0,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(3,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end

%plot labels
     xlabel('x','FontSize',20)
     ylabel('dx/dt','FontSize',20)
     text(-3,15,'dx/dt = r + x^2','Color','r','FontSize',20)
     text(1,15,strt,'FontSize',20)
     title('Saddle Node Bifurcation','FontSize',20)

%get movie frame
     h = gcf;
     M(i) = getframe(h,[5 5 480 380]);
     hold off

% end of loop
end

%play and save movie
     movie(M);
     movie2avi(M,'SaddleNode','fps',1);
