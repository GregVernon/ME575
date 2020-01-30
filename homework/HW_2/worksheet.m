%% Test Simple Quadratic
clear
n = 100;
fun = @(x) x.^2;
x0 = rand(n,1)*100;
fun(x0);
[x,res] = conjugateGradient(fun,x0,1e-5,1e4)


%% Test Matyas Function
clear
fun = @(x) 0.26 * (x(1)^2 + x(2)^2) - 0.48 * prod(x); 
x0 = [20 -20];
fun(x0);
[x,res] = conjugateGradient(fun,x0,1e-5,1e4)

%% Test Rosenbrock Function
clear
fun = @(x) (1-x(1))^2 + 100 * (x(2) - x(1)^2)^2;
x0 = [10,10];
fun(x0);
[x,res] = conjugateGradient(fun,x0,1e-9,1e5)