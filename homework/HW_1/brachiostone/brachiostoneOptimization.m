clear
close all

nx = 64;
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;

x = linspace(xMin,xMax,nx);
dx = x(2) - x(1);
y = interp1([0 1],[1 0],x);
y(1) = yMax;
y(end) = yMin;

mu = 0.0;
g = 9.81;
m = 1;

% figure
hold on
plot(x,y,'--o')

% Solve via Bernoulli Method
options = optimset;
options.MaxFunEvals = 1e6;
% options.PlotFcns = ["optimplotx","optimplotfval","optimplotfunccount"];

y = fminunc(@(y)belgundu(y,dx,m,g,mu),y,options);

plot(x,y,'-o')

% Solve via Ning method
H = yMax - yMin;
options = optimset;
options.MaxFunEvals = 1e6;
% options.PlotFcns = ["optimplotx","optimplotfval","optimplotfunccount"];
y = fminunc(@(y)ning(y,x,H,mu),y,options);
plot(x,y,'-o')
