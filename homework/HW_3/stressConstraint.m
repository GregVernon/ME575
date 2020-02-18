function [conINEQ,conEQ] = stressConstraint(area)
[~,stress] = truss(area);
yieldStress = 25e3 * ones(10,1);
yieldStress(9) = 75e3;

conINEQ = [stress - yieldStress; -yieldStress - stress];
conEQ = [];
end
