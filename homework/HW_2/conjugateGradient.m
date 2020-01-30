function [x,nl_res] = conjugateGradient(fun,x0,nl_tol,max_iter)
x0 = reshape(x0,length(x0),1);
fval = fun(x0);
x = x0;
iter = 0;
nl_res = inf;
while nl_res > nl_tol && iter <= max_iter
    iter = iter + 1;
    
    % Compute Search Direction
    Gk = computeGradient_FDM(fun,x);
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
        fval_test = fun(x_test);
        condition = fval_test <= fval_last + mu * stepSize * Pk;
        if condition == true
            isSufficient = true;
        else
            stepSize = stepSize / 2;
        end
    end
    x = x_test;
    fval = fun(x);
    
    
    % Set values for next iteration
    Gk_last = Gk;
    Pk_last = Pk;
    fval_last = fval;
    
    nl_res = norm(Gk,2);
end



end

function g = computeGradient_FDM(fun,x0)
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
        g(dim) = (fun(x0+DX) - fun(x0-DX))./(2*dx);
    elseif singleOutput == false
        g = g + (fun(x0+DX) - fun(x0-DX))./(2*dx);
    end
end
end