function [J,t] = ningBezier(dofNodes,domain,mu,B)
% dofNodes = [x2,y2,x3,y3,x(n-1),y(n-1)]
% domain = [x0,y0,xn,yn]
bezNodes = [domain(1) domain(2) dofNodes domain(3) domain(4)];
% nodes = [x1,y1,x2,y2,x3,y3,...xn,yn]
bezNodes = transpose(reshape(bezNodes,2,length(bezNodes)/2));
bez = bezier([bezNodes(:,1) bezNodes(:,2)],B);
x = abs(bez(:,1));
y = bez(:,2);

H = domain(2) - domain(4);
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