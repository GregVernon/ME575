function [conINEQ] = bezierConstraint(bezNodes)
% bezNodes_x = [domain(1) dofNodes(1:2:end) domain(3)];
bezNodes_x = bezNodes(1:2:end);
conINEQ = bezNodes_x(2:end) - bezNodes_x(1:end-1);
conINEQ = conINEQ(conINEQ<0);
conINEQ = length(bezNodes_x) * sum(abs(conINEQ));
% conEQ = [];
end