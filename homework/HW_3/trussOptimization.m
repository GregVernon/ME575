x0 = ones(10,1);
minArea = 0.1*ones(size(x0));
[mass,stress] = truss(x0);

options = optimoptions('fmincon');
% options.Display = "iter-detailed";
options.MaxFunctionEvaluations = 1e6;
options.MaxIterations = 1e6;
options.OptimalityTolerance = 1e-12;
options.ConstraintTolerance = 1e-12;
% options.PlotFcns = {"optimplotx","optimplotfval","optimplotfunccount"};

% Internal Calculations
options.SpecifyObjectiveGradient = false;
options.SpecifyConstraintGradient = false;
[A,fval] = fmincon(@(A)truss(A),x0,[],[],[],[],minArea,[],@(A)stressConstraint(A),options);

% Forward Difference
options.SpecifyObjectiveGradient = true;
options.SpecifyConstraintGradient = true;
gradMethod = "Forward-Difference";
[A,fval] = fmincon(@(x0)truss(x0,gradMethod),x0,[],[],[],[],minArea,[],@(x0)stressConstraint(x0,gradMethod),options);

% Centered Difference
options.SpecifyObjectiveGradient = true;
options.SpecifyConstraintGradient = true;
gradMethod = "Centered-Difference";
[A,fval] = fmincon(@(x0)truss(x0,gradMethod),x0,[],[],[],[],minArea,[],@(x0)stressConstraint(x0,gradMethod),options);

% Complex Step
options.SpecifyObjectiveGradient = true;
options.SpecifyConstraintGradient = true;
gradMethod = "Complex-Step";
[A,fval] = fmincon(@(x0)truss(x0,gradMethod),x0,[],[],[],[],minArea,[],@(x0)stressConstraint(x0,gradMethod),options);
