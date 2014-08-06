
% Range of x
x = -3:0.1:3

% Input value of r
r = -1.0
strr = num2str(r);
strl = ['r=' strr];

% Calculate h
h = x.^3 - r.*x;

% Stable branch (all stable): h starts very low and increases until hc. 

% Plot the stable branch
plot(h,x)
axis([-3 3 -2 2])
hold on

% Plot the x- and h- axes.
t = -5:0.1:5
line(0,t)
line(t,0)

% Label axes and add comments
xlabel('h','FontSize',20)
     ylabel('x_*','FontSize',20)
     text(-2,1,strl,'Fontsize',20)

