function [J,t] = ning(y,x,H,mu)
y(end) = 0;
y(1) = y(end) + H;

J = zeros(length(x)-1,1);
for ii = 1:length(x)-1
    dx = x(ii+1) - x(ii);
    dy = y(ii+1) - y(ii);
    J(ii) = sqrt(dx^2 + dy^2) / (sqrt(H - y(ii+1) - mu*x(ii+1)) + sqrt(H - y(ii) - mu*x(ii)));
end
J = sum(J);
t = computeElapsedTime(J,9.81);
end

function t = computeElapsedTime(J,g)
    t = sqrt(2/g) * J;
end