%% Animation for a different Saddle-Node

     xmin = -3;
     xmax = 3;
     x = xmin:0.1:xmax
     ax = (0.0).*x;

     rmin = -4
     rmax = 4
     imax = 41

     ymin = -10;
     ymax = 10;
     ky = ymin:0.1:ymax
     ay = (0.0).*ky

for i = 1:imax
%calculate r value
     r = rmin + (i-1)*(rmax-rmin)/(imax-1)
%create title string
     strr = num2str(r);
     strt = ['r=' strr]; 

%calculate dxdt
     y1 = r - x
     y2 = exp(-x)

%plot x-axis
     plot(x,ax,'--');
     axis([xmin xmax ymin ymax])
     hold on
%plot y-axis
     plot(ay,ky,'--')

%plot f(x)
     plot(x,y1,'r');
     plot(x,y2,'g');

%plot labels
     xlabel('x','FontSize',20)
     ylabel('dx/dt','FontSize',20)
     text(-2,-8,'dx/dt = r - x','Color','r','FontSize',20)
     text(-2,-6,'dx/dt = e^{-x}','Color','g','FontSize',20)
     text(2,5,strt,'FontSize',20)


%getframe
     h = gcf;
     M(i) = getframe(h,[5 5 480 380]);
     hold off
end

%play and save movie
     movie(M);
     movie2avi(M,'Graphbif','fps',2);
