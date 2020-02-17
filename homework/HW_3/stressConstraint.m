function [conINEQ,conEQ,gradConINEQ,gradConEQ] = stressConstraint(area, gradMethod)
[~,~,stress] = truss(area);
yieldStress = 25e3 * ones(10,1);
yieldStress(9) = 75e3;

conINEQ = stress.^2 - yieldStress.^2;
conEQ = [];

if exist('gradMethod','var')
    if     strcmpi(gradMethod, "Forward-Difference")
        dx = 1e-6 * norm(area,inf);
        [conINEQ_1,conEQ_1] = stressConstraint(area+dx);
        gradConINEQ = (conINEQ_1 - conINEQ) / dx;
        gradConEQ = (conEQ_1 - conEQ) / dx;
    elseif strcmpi(gradMethod, "Centered-Difference")
        dx = 1e-8 * norm(area,inf);
        [conINEQ_1,conEQ_1] = stressConstraint(area+dx);
        [conINEQ_2,conEQ_2] = stressConstraint(area-dx);
        gradConINEQ = (conINEQ_1 - conINEQ_2) / (2*dx);
        gradConEQ = (conEQ_1 - conEQ_2) / dx;
    elseif strcmpi(gradMethod, "Complex-Step")
        dx = 1e-16 * 1i;
        [conINEQ_1,conEQ_1] = stressConstraint(area+dx);
        gradConINEQ = imag(conINEQ_1) / imag(dx);
        gradConEQ = imag(conEQ_1) / imag(dx);
    end
else
    gradConINEQ = [];
    gradConEQ = [];
end
end
