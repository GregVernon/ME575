function [xopt, fopt, outputs] = uncon(func, x0, epsilon_g, options)
%{
An algorithm for unconstrained optimization.

Parameters
----------
func : function handle
    function handle to a function of the form: [f, g] = func(x)
    where f is the function value and g is a *column* vector containing
    the gradient.  x are design variables only.
x0 : *column* vector
    starting point
epsilon_g : float
    convergence tolerance.  you should terminate when
    norm(g, inf) <= epsilon_g.  (the infinity norm of the gradient)
options : struct
    a struct containing options.  You can use this to try out different
    algorithm choices.  I will not pass anything in, so you should 
    check if nargin=3 and if so setup some defaults.

Outputs
-------
xopt : column vector
    the optimal solution
fopt : float
    the corresponding function value
outputs : struct
    other miscelaneous outputs that you might want, for example an array
    containing a convergence metric at each iteration.
%}

if nargin == 3  % means no options were passed in
    % set defaults here for how you want me to run it. 
    method = "BFGS"; % CG % fminunc
    gradFun = "Integrated"; % Complex-Step % Centered-Difference
    max_iter = 1e4;
end

if strcmpi(method,"CG")
    [xopt,fopt,res,iter,funEvals] = conjugateGradient(func,x0,"gradFun",gradFun,"nl_tol",epsilon_g,"max_iter",max_iter);
elseif strcmpi(method,"BFGS")
    [xopt,fopt,res,iter,funEvals] = BFGS(func,x0,"gradFun",gradFun,"nl_tol",epsilon_g,"max_iter",max_iter);
elseif strcmpi(method,"fminunc")
    options = optimoptions("fminunc");
    options.MaxIterations = max_iter;
    options.OptimalityTolerance = epsilon_g;
    [xopt,fopt,~,output] = fminunc(func,x0,options);
    res = output.firstorderopt;
    iter = output.iterations;
    funEvals = output.funcCount;
end

outputs.res = res;
outputs.iter = iter;
outputs.funEvals = funEvals;
end