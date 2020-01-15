function [x,y] = exactSolution(nPts)
a = 0.572917;
theta = linspace(0,2.412,nPts);
x = a*(theta - sin(theta));
y = 1-a*(1-cos(theta));
end