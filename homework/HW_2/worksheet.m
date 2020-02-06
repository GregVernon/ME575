%% Quadratic Conjugate Gradient
clear
n = 100;
fun = @(x) x.^2;
x0 = rand(1,n)*100;

% Analytic Gradient
gfun = @(x) 2.*x;
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Analytic Gradient -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end))) 
% Centered Difference
gfun = "Centered-Difference";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Centered-Difference -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))
% Complex Step
gfun = "Complex-Step";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Complex-Step -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))
%% Matyas Conjugate Gradient
clear
fun = @(x) 0.26 * (x(1)^2 + x(2)^2) - 0.48 * prod(x); 
x0 = [20 -20];

% Analytic Gradient
gfun = @(x) [(13*x(1))/25 - (12*x(2))/25; (13*x(2))/25 - (12*x(1))/25];
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Analytic Gradient -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end))) 
% Centered Difference
gfun = "Centered-Difference";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Centered-Difference -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))
% Complex Step
gfun = "Complex-Step";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Complex-Step -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))

%% Rosenbrock Conjugate Gradient
clear
fun = @(x) (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
x = sym('x',[1,2],'real');
f(x) = (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
g(x) = gradient(f);
% g = @(x) [2*x(1) - 400*x(1)*(- x(1)^2 + x(2)) - 2; -200*x(1)^2 + 200*x(2)]
x0 = [20,20];

% Analytic Gradient
gfun = g;
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Analytic Gradient -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))
% Centered Difference
gfun = "Centered-Difference";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Centered-Difference -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))
% Complex Step
gfun = "Complex-Step";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Complex-Step -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))

%% Brachistochrone Conjugate Gradient
fun = @(x) brachistochrone(x);
x0 = sort(linspace(0,1,128),"descend");

% Integrated Gradient
gfun = "Integrated";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Integrated -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))
% Centered-Difference
gfun = "Centered-Difference";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Centered-Difference -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))
% Complex Step
gfun = "Complex-Step";
[x,fopt,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Complex-Step -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))

%% Simple BFGS
clear
n = 100;
fun = @(x) x.^2;
x0 = rand(1,n)*100;

% Analytic Gradient
gfun = @(x) 2.*x;
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Analytic Gradient -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end))) 
% Centered Difference
gfun = "Centered-Difference";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Centered-Difference -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))
% Complex Step
gfun = "Complex-Step";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Complex-Step -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))

%% Matyas BFGS
clear
fun = @(x) 0.26 * (x(1)^2 + x(2)^2) - 0.48 * prod(x); 
x0 = [20 -20];

% Analytic Gradient
gfun = @(x) [(13*x(1))/25 - (12*x(2))/25; (13*x(2))/25 - (12*x(1))/25];
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Analytic Gradient -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end))) 
% Centered Difference
gfun = "Centered-Difference";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Centered-Difference -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))
% Complex Step
gfun = "Complex-Step";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-5,"max_iter",1e4);
disp("Complex-Step -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))

%% Rosenbrock BFGS
clear
fun = @(x) (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
x = sym('x',[1,2],'real');
f(x) = (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
g(x) = gradient(f);
% g = @(x) [2*x(1) - 400*x(1)*(- x(1)^2 + x(2)) - 2; -200*x(1)^2 + 200*x(2)]
x0 = [20,20];

% Analytic Gradient
gfun = g;
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Analytic Gradient -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))
% Centered Difference
gfun = "Centered-Difference";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Centered-Difference -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))
% Complex Step
gfun = "Complex-Step";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Complex-Step -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))

%% Brachistochrone BFGS
fun = @(x) brachistochrone(x);
x0 = sort(linspace(0,1,128),"descend");

% Integrated Gradient
gfun = "Integrated";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Integrated -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))
% Centered-Difference
gfun = "Centered-Difference";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Centered-Difference -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))
% Complex Step
gfun = "Complex-Step";
[x,fopt,res,iter,funEvals] = BFGS(fun,x0,"gradFun",gfun,"nl_tol",1e-4,"max_iter",1e4);
disp("Complex-Step -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)) + " Residual = " + num2str(res(end)))

%% Quadratic FMINUNC
fun = @(x) sum(x.^2);
n = 100;
x0 = rand(1,n)*100;

options = optimoptions("fminunc");
options.MaxIterations = 400;
options.OptimalityTolerance = 1e-4;
options.FiniteDifferenceType = "central";

[xopt,fopt,~,output] = fminunc(fun,x0,options);
res = output.firstorderopt;
iter = output.iterations;
funEvals = output.funcCount;
disp("fminunc -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))

%% Matyas FMINUNC
clear
fun = @(x) 0.26 * (x(1)^2 + x(2)^2) - 0.48 * prod(x); 
x0 = [20 -20];

options = optimoptions("fminunc");
options.MaxIterations = 400;
options.OptimalityTolerance = 1e-4;
options.FiniteDifferenceType = "central";

[xopt,fopt,~,output] = fminunc(fun,x0,options);
res = output.firstorderopt;
iter = output.iterations;
funEvals = output.funcCount;
disp("fminunc -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))

%% Rosenbrock FMINUNC
clear
fun = @(x) (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
x = sym('x',[1,2],'real');
f(x) = (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
g(x) = gradient(f);
x0 = [20,20];

options = optimoptions("fminunc");
options.MaxIterations = 400;
options.OptimalityTolerance = 1e-4;
options.FiniteDifferenceType = "central";

[xopt,fopt,~,output] = fminunc(fun,x0,options);
res = output.firstorderopt;
iter = output.iterations;
funEvals = output.funcCount;
disp("fminunc -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))

%% Brachistochrone FMINUNC
fun = @(x) brachistochrone(x);
x0 = sort(linspace(0,1,128)',"descend");

% Integrated
options = optimoptions("fminunc");
options.MaxIterations = 400;
options.OptimalityTolerance = 1e-4;
options.SpecifyObjectiveGradient = true;

[xopt,fopt,~,output] = fminunc(fun,x0,options);
res = output.firstorderopt;
iter = output.iterations;
funEvals = output.funcCount;
disp("fminunc -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))

% Centered Diff
options = optimoptions("fminunc");
options.MaxIterations = 400;
options.OptimalityTolerance = 1e-4;
options.FiniteDifferenceType = "central";
options.MaxFunctionEvaluations = 1e6;
[xopt,fopt,~,output] = fminunc(fun,x0,options);
res = output.firstorderopt;
iter = output.iterations;
funEvals = output.funcCount;
disp("fminunc -- " + "# Iterations = " + num2str(iter) + " # funEvals = " + num2str(funEvals) + " Residual = " + num2str(res(end)))
