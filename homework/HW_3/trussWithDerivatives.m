function [mass,gradMass,stress] = trussWithDerivatives(area,gradMethod)

% Compute values at current configuration
[mass, stress] = truss(area);

% Compute derivatives at current configuration
if exist("gradMethod","var")
    if strcmpi(gradMethod, "Auto-Diff")
        x = amatinit(area);
        [m,~] = trussAutoDiff(x);
        gradMass = ajac(transpose(m),0);
    elseif strcmpi(gradMethod,"Adjoint")
        [~,gradMass,~] = trussAdjoint(area);
    else
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
            end
        end
    end
else
    gradMass = [];
end


end

