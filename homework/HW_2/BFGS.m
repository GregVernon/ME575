function [x,fval,nl_res,iter,funEvals] = BFGS(fun,x0,NameValueArgs)
arguments
    fun
    x0
    NameValueArgs.gradFun = "Centered Finite Difference";
    NameValueArgs.nl_tol = 1e-8;
    NameValueArgs.max_iter = 1e3;
end
x0 = reshape(x0,length(x0),1); % Make sure we're using column vectors
gradFun = NameValueArgs.gradFun;
nl_tol = NameValueArgs.nl_tol;
max_iter = NameValueArgs.max_iter;
ndof = length(x0);
I = eye(ndof);

if contains(class(gradFun),"sym")
    gradFun = matlabFunction(gradFun);
end

funEvals = 0;
x = x0;
iter = 0;
nl_res = inf;
isStagnated = false;
while nl_res > nl_tol && iter <= max_iter && isStagnated == false
    iter = iter + 1;
    if iter == 1
        if strcmpi(gradFun, "Centered Finite Difference") == true
            fval = fun(x0);
            fval_last = fval;
            [Gk,funEvals] = computeGradient_FDM(fun,x,funEvals);
        elseif strcmpi(class(gradFun), "function_handle") == true
            fval = fun(x0);
            fval_last = fval;
            if nargin(gradFun) == 1
                Gk = gradFun(x);
            elseif nargin(gradFun) > 1
                X = num2cell(x);
                Gk = gradFun(X{:});
            end
        elseif strcmpi(gradFun,"integrated") == true
            [fval_last,Gk] = fun(x0);
        end
        Vk = I;
    end
    searchDirection = -Vk * Gk;
    stepSize = 1;
    % Sufficient decrease
    isSufficient = false;
    mu = 1e-4;
    lineSearch_iter = 0;
    while isSufficient == false
        lineSearch_iter = lineSearch_iter + 1;
        x_test = x + stepSize * searchDirection;
        fval_test = fun(x_test); funEvals = funEvals + 1;
        condition = fval_test <= fval_last + mu * stepSize * searchDirection;
        if all(condition == true)
            isSufficient = true;
        else
            stepSize = stepSize / 2;
        end
    end
    sk = stepSize * searchDirection;
    x_next = sk + x;
    if strcmpi(gradFun, "Centered Finite Difference") == true
        fval = fun(x_next);
        [Gk_next,funEvals] = computeGradient_FDM(fun,x_next,funEvals);
    elseif strcmpi(class(gradFun), "function_handle") == true
        fval = fun(x_next);
        if nargin(gradFun) == 1
            Gk_next = gradFun(x_next);
        elseif nargin(gradFun) > 1
            X = num2cell(x_next);
            Gk_next = gradFun(X{:});
        end
    elseif strcmpi(gradFun,"integrated") == true
        [fval,Gk] = fun(x0);
    end
    yk = Gk_next - Gk;
    
    Vk_next = [I - (sk*transpose(yk))/(transpose(sk)*yk)] * Vk * [I - (yk * transpose(sk))/(transpose(sk)*yk)] + (sk*transpose(sk)) / (transpose(sk) * yk);
    
    
    fval_last = fval;
    x = x_next;
    Gk = Gk_next;
    Vk = Vk_next;
    nl_res = norm(Gk,2);
    isStagnated = any(isnan(Vk(:)));
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
dx = 1e-14;
for dim = 1:length(x0)
    DX = zeros(length(x0),1);
    DX(dim) = 1i * (x0(dim) * dx + eps);
    % Central Difference Approximation
    if singleOutput == true
        g(dim) = (fun(x0+DX) - fun(x0-DX))./(2*DX(dim)); funEvals = funEvals + 2;
    elseif singleOutput == false
        g = g + (fun(x0+DX) - fun(x0-DX))./(2*DX(dim)); funEvals = funEvals + 2;
    end
end
end