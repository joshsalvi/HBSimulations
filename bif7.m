%% Animation for a zoomed in saddlenode bifurcation

     xmin = -0.5;
     xmax = 0.5;
     x = xmin:0.001:xmax
     ax = (0.0).*x;

     rmin = 1-0.2
     rmax = 1+0.2
     imax = 21

     ymin = -0.5;
     ymax = 0.5;
     ky = ymin:0.1:ymax
     ay = (0.0).*ky

for i = 1:imax
%calculate r value
     r = rmin + (i-1)*(rmax-rmin)/(imax-1)
%create title string
     strr = num2str(r);
     strt = ['r=' strr]; 

%calculate dxdt
     y = r - x - exp(-x);

%plot x-axis
     plot(x,ax,'--');
     axis([xmin xmax ymin ymax])
     hold on
%plot y-axis
     plot(ay,ky,'--')

%plot f(x)
     plot(x,y,'r');


%plot labels
     xlabel('x','FontSize',20)
     ylabel('dx/dt','FontSize',20)
     text(-0.3,0.4,'dx/dt = r - x - e^{-x}','Color','r','FontSize',20)
     text(0.3,0.4,strt,'FontSize',20)


%getframe
     h = gcf;
     M(i) = getframe(h,[5 5 480 380]);
     hold off
end

%play and save movie
     movie(M);
     movie2avi(M,'Graphbif3','fps',2);
