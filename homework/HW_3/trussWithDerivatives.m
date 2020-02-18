function [mass,gradMass,stress] = trussWithDerivatives(area,gradMethod)

% Compute values at current configuration
[mass, stress] = truss(area);

% Compute derivatives at current configuration
if exist("gradMethod","var")
    gradMass = zeros(length(area),1);
    for ii = 1:length(area)
        dx = zeros(size(area));
        if     strcmpi(gradMethod, "Forward-Difference")
            dx(ii) = 1e-6 * norm(area,2);
            mass_0 = truss(area);
            mass_1 = truss(area+dx);
            gradMass(ii) = (mass_1 - mass_0) / dx(ii);
        elseif strcmpi(gradMethod, "Centered-Difference")
            dx(ii) = 1e-8 * norm(area,2);
            mass_1 = truss(area+dx);
            mass_2 = truss(area-dx);
            gradMass(ii) = (mass_1 - mass_2) / (2*dx(ii));
        elseif strcmpi(gradMethod, "Complex-Step")
            dx(ii) = 1e-12 * norm(area,2) * 1i;
            mass_1 = truss(area+dx);
            gradMass(ii) = imag(mass_1) / imag(dx(ii));
        elseif strcmpi(gradMethod, "Auto-Diff")
            x = amatinit(area);
            [m,s] = trussAutoDiff(x);
            gradMass = ajac(m,0);
%             dg = ajac(s,0);
%             area = ainit(area,1);
%             gradMass = adiff(trussAutoDiff(area));
%             gradMass = gradMass.c(:,2);
        end
    end
else
    gradMass = [];
end


end

