%% Animation for a Subcritical Pitchfork bifurcation

     xmin = -3;
     xmax  = 3;
     x = xmin:0.1:xmax
     ax = (0.0).*x;

     rmin = -4
     rmax = 4
     imax = 21

     ymin = -20;
     ymax = 20; 
     ky = ymin:0.1:ymax
     ay = (0.0).*ky

for i = 1:imax
%calculate r value
     r = rmin + (i-1)*(rmax-rmin)/(imax-1)
%create title string
     strr = num2str(r);
     strt = ['r=' strr]; 

%calculate dxdt
     y = r.*x + x.^3

%plot x-axis
     plot(x,ax,'--');
     axis([xmin xmax ymin ymax])
     hold on
%plot y-axis
     plot(ay,ky,'--')

%plot f(x)
     plot(x,y,'r');

%plot fixed points
if (r < 0) 
     plot(0,0,'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
     plot((-r).^(0.5),0,'o','MarkerSize',10,'MarkerEdgeColor','k');
     plot(-(-r).^(0.5),0,'o','MarkerSize',10,'MarkerEdgeColor','k');
end
if(r == 0.0) 
     plot(0,0,'*','MarkerSize',10,'MarkerEdgeColor','k');
end
if (r > 0) 
     plot(0,0,'o','MarkerSize',10,'MarkerEdgeColor','k');
end

%plot flow direction
if(r < 0)
     plot(-2.5,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
     plot(-(-r).^(0.5)/2,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
     plot((-r).^(0.5)/2,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
   plot(2.5,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end

if(r >= 0) 
     plot(-2.5,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');     
     plot(2.5,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end

%plot labels
     xlabel('x','FontSize',20)
     ylabel('dx/dt','FontSize',20)
     text(-2.5,15,'dx/dt = rx + x^3','Color','r','FontSize',20)
     text(1,15,strt,'FontSize',20)
     title('Subcritical Pitchfork Bifurcation','FontSize',20)

%getframe
     h = gcf;
     M(i) = getframe(h,[5 5 480 380]);
     hold off
end

%play and save movie
     movie(M);
     movie2avi(M,'SubPitchfork','fps',1);
