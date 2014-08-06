%% Animation for a Transcritical Bifurcation (same structure as Saddle-node, limited annotations for this program)

     xmin = -6;
     xmax = 6;
     x = xmin:0.1:xmax
     ax = (0.0).*x;

     rmin = -4
     rmax = 4
     imax = 21

     ymin = -6 ;
     ymax = 6;
     ky = ymin:0.1:ymax
     ay = (0.0).*ky

for i = 1:imax
%calculate r value
     r = rmin + (i-1)*(rmax-rmin)/(imax-1)
%create title string
     strr = num2str(r);
     strt = ['r=' strr]; 

%calculate dxdt
     y = r.*x - x.^2

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
     plot(r,0,'o','MarkerSize',10,'MarkerEdgeColor','k');
     plot(0,0,'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end
if(r == 0.0) 
     plot((-r).^(0.5),0,'*','MarkerSize',10,'MarkerEdgeColor','k');
end
if (r > 0) 
     plot(r,0,'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
     plot(0,0,'o','MarkerSize',10,'MarkerEdgeColor','k');
end

%plot flow direction
if(r < 0)
     plot(-4.5,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
     plot(0.5.*r,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(4.5,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end

if(r == 0) 
     plot(-4.5,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');     
plot(4.5,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end

if(r > 0) 
     plot(-4.5,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');     plot(0.5.*r,0,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(4.5,0,'<','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
end

%plot labels
     xlabel('x','FontSize',20)
     ylabel('dx/dt','FontSize',20)
     text(-5,5,'dx/dt = rx - x^2','Color','r','FontSize',20)
     text(2,5,strt,'FontSize',20)
     title('Transcritical Bifurcation','FontSize',20)

%getframe
     h = gcf;
     M(i) = getframe(h,[5 5 480 380]);
     hold off
end

%play and save movie
     movie(M);
     movie2avi(M,'Transcritical','fps',1);
