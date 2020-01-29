%% Test Simple Quadratic
n = 100;
fun = @(x) x.^2;
x0 = rand(n,1)*100;
fun(x0);
[x,res] = nl_unconstrainedOptimize(fun,x0)


%%