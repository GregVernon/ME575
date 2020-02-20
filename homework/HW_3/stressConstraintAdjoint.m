function [conINEQ,conEQ,gradConINEQ,gradConEQ] = stressConstraintAdjoint(area)
[~,~,stress,DG] = trussAdjoint(area);
yieldStress = 25e3 * ones(10,1);
yieldStress(9) = 75e3;

conINEQ = [stress - yieldStress; -yieldStress - stress];
conEQ = [];
gradConINEQ = [DG -DG];
gradConEQ = [];
end
