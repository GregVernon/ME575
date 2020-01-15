clear

nx = 16;
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
options.MaxIter = 1e4;
% options.PlotFcns = ["optimplotx","optimplotfval","optimplotfunccount"];

y = fminunc(@(y)belgundu(y,dx,m,g,mu),y,options);


%% Ning Method
clear

nx = 256;
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

% Solve via Ning method
H = yMax - yMin;
options = optimset;
options.MaxFunEvals = 1e6;
options.MaxIter = 1e4;
% options.PlotFcns = ["optimplotx","optimplotfval","optimplotfunccount"];
y = fminunc(@(y)ning(y,x,H,mu),y,options);

% Exact Solution
[xExact, yExact] = exactSolution(mu,1e4);

figure 
hold on
plot(xExact,yExact)
plot(x,y)

%% Bezier Geometry
clear
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
domain = [xMin yMax xMax yMin];

bezDegree = 3;
bezDomain = [0 1];
nPts = 64;
bezVariate = linspace(bezDomain(1), bezDomain(2),nPts);
B = bernsteinPolynomial(bezDegree,bezDomain,bezVariate);
nDOFS = (bezDegree+1) - 2;

x0 = linspace(xMin,xMax,nDOFS+2);
x0(1) = [];
x0(end) = [];
y0 = interp1([0 1],[1 0],x0);

dofNodes = [x0' y0'];
dofNodes = reshape(dofNodes',1,numel(dofNodes));

mu = 0.;

options = optimset;
options.MaxFunEvals = 1e6;
dofNodes = fminunc(@(dofNodes)ningBezier(dofNodes,domain,mu,B),dofNodes,options);

bezNodes = [domain(1) domain(2) dofNodes domain(3) domain(4)];
bezNodes = transpose(reshape(bezNodes,2,length(bezNodes)/2));
f = bezier([bezNodes(:,1) bezNodes(:,2)],B);
x = f(:,1);
y = f(:,2);

% Exact Solution
[xExact, yExact] = exactSolution(mu,1e4);

figure;
hold on
plot(exactX,exactY);
plot(x,y,'-o')
plot(bezNodes(:,1),bezNodes(:,2),'LineStyle','-','Color','k','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')

%% Scaling of Ning Method
clear
maxNumCompThreads(1) % Run Single-Thread

nx = 2.^[2:8];
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;

mu = 0;

for iter = 1:length(nx)
    x = linspace(xMin,xMax,nx(iter));
    dx = x(2) - x(1);
    y = interp1([0 1],[1 0],x);
    y(1) = yMax;
    y(end) = yMin;

    H = yMax - yMin;
    options = optimset;
    options.MaxFunEvals = 1e6;
    options.MaxIter = 1e4;
    
    [~,fval,~,output] = fminunc(@(y)ning(y,x,H,mu),y,options);
    solverData(iter) = output;
    
    benchFun = @() fminunc(@(y)ning(y,x,H,mu),y,options);
    t(iter) = timeit(benchFun);
end
plot(nx,t,"Color",'k',"LineStyle",'-',"Marker",'o',"MarkerFaceColor",'k',"MarkerEdgeColor",'k')
xlabel("Number of Points (nDOF + 2)")
ylabel("Solve Time (sec)")

%% Scaling of Bezier-enhanced Ning
clear
maxNumCompThreads(1) % Run Single-Thread

xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
domain = [xMin yMax xMax yMin];

bezDegree = 4;
bezDomain = [0 1];
nPts = 2.^[2:12];
for iter = 1:length(nPts)
    bezVariate = linspace(bezDomain(1), bezDomain(2),nPts(iter));
    B = bernsteinPolynomial(bezDegree,bezDomain,bezVariate);
    nDOFS = (bezDegree+1) - 2;
    
    x0 = linspace(xMin,xMax,nDOFS+2);
    x0(1) = [];
    x0(end) = [];
    y0 = interp1([0 1],[1 0],x0);

    dofNodes = [x0' y0'];
    dofNodes = reshape(dofNodes',1,numel(dofNodes));

    mu = 0.;

    options = optimset;
    options.MaxFunEvals = 1e6;
    [~,fval,~,output] = fminunc(@(dofNodes)ningBezier(dofNodes,domain,mu,B),dofNodes,options);
    solverData(iter) = output;
    
    benchFun = @() fminunc(@(dofNodes)ningBezier(dofNodes,domain,mu,B),dofNodes,options);
    t(iter) = timeit(benchFun);
end
figure
plot(nPts,t,"Color",'k',"LineStyle",'-',"Marker",'o',"MarkerFaceColor",'k',"MarkerEdgeColor",'k')
xlabel("Number of Points")
ylabel("Solve Time (sec)")

%% "Coarse-Grid Preconditioning"
clear

figure 
hold on

xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;

nx = 2.^[2:8];
for iter = 1:length(nx)
    if iter == 1
        x = linspace(xMin,xMax,nx(iter));
        dx = x(2) - x(1);
        y = interp1([0 1],[1 0],x);
    else
        xOld = x;
        yOld = y;
        x = linspace(xMin,xMax,nx(iter));
        y = interp1(xOld,yOld,x);
    end
    y(1) = yMax;
    y(end) = yMin;
    
    mu = 0.0;
    
%     % figure
    plot(x,y,'--o')
    
    % Solve via Bernoulli Method
    options = optimset;
    options.MaxFunEvals = 1e6;
    options.MaxIter = 1e4;
    
    H = yMax - yMin;
    options = optimset;
    options.MaxFunEvals = 1e6;
    % options.PlotFcns = ["optimplotx","optimplotfval","optimplotfunccount"];
    y = fminunc(@(y)ning(y,x,H,mu),y,options);
%     plot(x,y,'-')
end