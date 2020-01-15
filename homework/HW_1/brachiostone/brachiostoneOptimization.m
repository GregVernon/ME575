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

%% Bezier Geometry
clear
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
domain = [xMin yMax xMax yMin];

nDOFS = 2;
x0 = linspace(xMin,xMax,nDOFS+2);
x0(1) = [];
x0(end) = [];
y0 = interp1([0 1],[1 0],x0);

dofNodes = [x0' y0'];
dofNodes = reshape(dofNodes',1,numel(dofNodes));

nPts = 64;
mu = 0.;

options = optimset;
options.MaxFunEvals = 1e6;
dofNodes = fminunc(@(dofNodes)ningBezier(dofNodes,domain,mu,nPts),dofNodes,options);

bezNodes = [domain(1) domain(2) dofNodes domain(3) domain(4)];
bezNodes = transpose(reshape(bezNodes,2,length(bezNodes)/2));
[f,B] = bezier([bezNodes(:,1) bezNodes(:,2)],nPts);
x = f(:,1);
y = f(:,2);
plot(x,y,'-o')

plot(x,y,'-o')