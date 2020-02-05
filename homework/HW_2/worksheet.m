%% Test Simple Quadratic
clear
n = 100;
fun = @(x) x.^2;
x0 = rand(1,n)*100;
fun(x0);
[x,res,funEvals] = conjugateGradient(fun,x0,"nl_tol",1e-5,"max_iter",1e4)


%% Test Matyas Function
clear
fun = @(x) 0.26 * (x(1)^2 + x(2)^2) - 0.48 * prod(x); 
x0 = [20 -20];
fun(x0);
[x,res,funEvals] = conjugateGradient(fun,x0)

%% Test Rosenbrock Function
clear
clc
fun = @(x) (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
x = sym('x',[1,2],'real');
f(x) = (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
g(x) = gradient(f);
% g = @(x) [2*x(1) - 400*x(1)*(- x(1)^2 + x(2)) - 2; -200*x(1)^2 + 200*x(2)]
x0 = [10,10];
fun(x0);
[x,res,iter,funEvals] = conjugateGradient(fun,x0,"gradFun",g,"searchDirection","Fletcher-Reeves","nl_tol",1e-9,"max_iter",1e4)

%% Test Simple Quadratic
clear
n = 100;
fun = @(x) x.^2;
x0 = rand(1,n)*100;
fun(x0);
[x,res,funEvals] = BFGS(fun,x0,"nl_tol",1e-5,"max_iter",1e4)

%% Test Matyas Function
clear
fun = @(x) 0.26 * (x(1)^2 + x(2)^2) - 0.48 * prod(x); 
x0 = [20 -20];
fun(x0);
[x,res,funEvals] = BFGS(fun,x0)

%% Test Rosenbrock Function
clear
clc
fun = @(x) (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
x = sym('x',[1,2],'real');
f(x) = (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
g(x) = gradient(f);
% g = @(x) [2*x(1) - 400*x(1)*(- x(1)^2 + x(2)) - 2; -200*x(1)^2 + 200*x(2)]
x0 = [4,4];
fun(x0);
[x,res,iter,funEvals] = BFGS(fun,x0,"gradFun",g,"nl_tol",1e-9,"max_iter",1e4)

