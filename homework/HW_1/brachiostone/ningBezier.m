function J = ningBezier(dofNodes,domain,mu,nPts)
% dofNodes = [x2,y2,x3,y3,x(n-1),y(n-1)]
% domain = [x0,y0,xn,yn]
bezNodes = [domain(1) domain(2) dofNodes domain(3) domain(4)];
% nodes = [x1,y1,x2,y2,x3,y3,...xn,yn]
bezNodes = transpose(reshape(bezNodes,2,length(bezNodes)/2));
[f,B] = bezier([bezNodes(:,1) bezNodes(:,2)],nPts);
x = f(:,1);
y = f(:,2);
% y(end) = 0;
% y(1) = y(end) + H;

H = domain(2) - domain(4);
J = zeros(length(x)-1,1);
for ii = 1:length(x)-1
    dx = x(ii+1) - x(ii);
    dy = y(ii+1) - y(ii);
    J(ii) = sqrt(dx^2 + dy^2) / (sqrt(H - y(ii+1) - mu*x(ii+1)) + sqrt(H - y(ii) - mu*x(ii)));
end
J = sum(J);
end