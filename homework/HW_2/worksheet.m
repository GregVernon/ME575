%% Test Simple Quadratic
clear
n = 100;
fun = @(x) x.^2;
x0 = rand(n,1)*100;
fun(x0);
[x,res] = conjugateGradient(fun,x0,1e-5,1e4)


%%