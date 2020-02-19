function [conINEQ,conEQ,gradConINEQ,gradConEQ] = stressConstraintWithDerivatives(area,gradMethod)

% Compute constraints at current configuration
[conINEQ,conEQ] = stressConstraint(area);

% Compute constraint derivatives at current configuration
if exist('gradMethod','var')
    if strcmpi(gradMethod,"Auto-Diff")
        x = amatinit(area);
        [conINEQ_tmp,conEQ_tmp] = stressConstraintAutoDiff(x);
        gradConINEQ = transpose(ajac(conINEQ_tmp,0));
        gradConEQ = [];
        
    else
        gradConINEQ = zeros(length(area),length(conINEQ));
        for ii = 1:length(area)
            dx = zeros(length(area),1);
            if     strcmpi(gradMethod, "Forward-Difference")
                dx(ii) = 1e-6 * area(ii);
                [conINEQ_1,conEQ_1] = stressConstraint(area+dx);
                gradConINEQ(ii,:) = (conINEQ_1 - conINEQ) / dx(ii);
                gradConEQ = []; %(conEQ_1 - conEQ) / dx(ii);
            elseif strcmpi(gradMethod, "Centered-Difference")
                dx(ii) = 1e-8 * area(ii);
                [conINEQ_1,conEQ_1] = stressConstraint(area+dx);
                [conINEQ_2,conEQ_2] = stressConstraint(area-dx);
                gradConINEQ(ii,:) = (conINEQ_1 - conINEQ_2) / (2*dx(ii));
                gradConEQ = [];
            elseif strcmpi(gradMethod, "Complex-Step")
                dx(ii) = 1e-16 * 1i;
                [conINEQ_1,conEQ_1] = stressConstraint(area+dx);
                gradConINEQ(ii,:) = imag(conINEQ_1) / imag(dx(ii));
                gradConEQ = [];
            end
        end
    end
else
    gradConINEQ = [];
    gradConEQ = [];
end
end