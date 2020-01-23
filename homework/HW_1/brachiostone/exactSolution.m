function [x,y] = exactSolution(x,fCoeff,xType)
if xType == "x-locations"
    a = 0.572917;
    % Solve for theta
    options = optimoptions('fminunc');
    options.MaxFunEvals = 1e6;
    options.MaxIter = 1e4;
    options.TolFun = 1e-14;
    options.TolX = 1e-14;
    options.OptimalityTolerance = 1e-14;
    options.StepTolerance = 1e-14;
    theta = reshape(linspace(0,2.412,length(x)),size(x));
    theta = fminunc(@(theta)thetaObjFun(theta,a,fCoeff,x),theta,options);
    y = 1-a*(1-cos(theta) + fCoeff*(theta + sin(theta)));
elseif xType == "#Pts"
    a = 0.572917;
    theta = linspace(0,2.412,nPts);
    x = a*(theta - sin(theta) + fCoeff*(1-cos(theta)));
    y = 1-a*(1-cos(theta) + fCoeff*(theta + sin(theta)));
end
end

function objVal = thetaObjFun(theta,a,fCoeff,x)
    objVal = sum((a*(theta - sin(theta) + fCoeff*(1-cos(theta))) - x).^2);
end