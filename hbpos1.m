function ddx = hbpos1(time,kappags,xa,po,xbase,x)
% x'' =
% (-1/(masshb+masssf))*((xihb+xisf)*x'+ngs*gamma*kappags(t)*(gamma*x-xa(t)+xcinit-po(t)*d)+kappasp*(x-xs)-kappasf*(xbase(t)-x);
ddx = zeros(2,1);
ddx(1) = ddx(2);
ddx(2) = -1/(masshb+masssf))*((xihb+xisf)*x(2)+ngs*gamma*kappags(1)*(gamma*x-xa(1)+xcinit-po(1)*d)+kappasp*(x(1)-xs)-kappasf*(xbase(1)-x);
end
