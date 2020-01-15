function [x,y] = exactSolution(fCoeff,nPts)
a = 0.572917;
theta = linspace(0,2.412,nPts);
x = a*(theta - sin(theta) + fCoeff*(1-cos(theta)));
y = 1-a*(1-cos(theta) + fCoeff*(theta + sin(theta)));
end