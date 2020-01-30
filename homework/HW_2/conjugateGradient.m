function [x,nl_res,iter,funEvals] = conjugateGradient(fun,x0,NameValueArgs)
arguments
    fun
    x0
    NameValueArgs.gradFun = "Centered Finite Difference";
    NameValueArgs.nl_tol = 1e-8;
    NameValueArgs.max_iter = 1e3;
end
x0 = reshape(x0,length(x0),1);
gradFun = NameValueArgs.gradFun;
nl_tol = NameValueArgs.nl_tol;
max_iter = NameValueArgs.max_iter;

if contains(class(gradFun),"sym")
    gradFun = matlabFunction(gradFun);
end

funEvals = 0;
fval = fun(x0); funEvals = funEvals + 1;
x = x0;
iter = 0;
nl_res = inf;
while nl_res > nl_tol && iter <= max_iter
    iter = iter + 1;
    
    % Compute Search Direction
    if strcmpi(gradFun, "Centered Finite Difference") == true
        [Gk,funEvals] = computeGradient_FDM(fun,x,funEvals);
    elseif strcmpi(class(gradFun), "function_handle") == true
        if nargin(gradFun) == 1
            Gk = gradFun(x);
        elseif nargin(gradFun) > 1
            X = num2cell(x);
            Gk = gradFun(X{:});
        end
    end
    
    if iter == 1
        % Do Steepest Descent
        fval_last = fval;
        Pk = -Gk;
    else
        % Do Conjugate Gradient
        Bk = (transpose(Gk) * Gk) / (transpose(Gk_last) * Gk_last);
        Pk = -Gk + Bk * Pk_last;
    end
    
    % Compute Initial Step Size
    stepSize = 1; %stepSize * (norm(Gk_last,2)/norm(Gk,2)).^2;
    
    % Sufficient decrease
    isSufficient = false;
    mu = 1e-4;
    while isSufficient == false
        x_test = x + stepSize * Pk;
        fval_test = fun(x_test); funEvals = funEvals + 1;
        condition = fval_test <= fval_last + mu * stepSize * Pk;
        if condition == true
            isSufficient = true;
        else
            stepSize = stepSize / 2;
        end
    end
    x = x_test;
    fval = fun(x); funEvals = funEvals + 1;
    
    
    % Set values for next iteration
    Gk_last = Gk;
    Pk_last = Pk;
    fval_last = fval;
    
    nl_res = norm(Gk,2);
end



end

function [g,funEvals] = computeGradient_FDM(fun,x0, funEvals)
% Test if function returns 1 or multiple results
if length(fun(x0)) == 1
    singleOutput = true;
else
    singleOutput = false;
end
g = zeros(length(x0),1);
dx = 1e-6;
for dim = 1:length(x0)
    DX = zeros(length(x0),1);
    DX(dim) = dx;
    % Central Difference Approximation
    if singleOutput == true
        g(dim) = (fun(x0+DX) - fun(x0-DX))./(2*dx); funEvals = funEvals + 2;
    elseif singleOutput == false
        g = g + (fun(x0+DX) - fun(x0-DX))./(2*dx); funEvals = funEvals + 2;
    end
end
end