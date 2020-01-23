A = ones(10,1);
minArea = 0.1*ones(size(A));
[mass,stress] = truss(A);

options = optimoptions('fmincon');
options.Display = "iter-detailed";
options.MaxFunctionEvaluations = 1e6;
options.MaxIterations = 1e6;
options.OptimalityTolerance = 1e-12;
options.ConstraintTolerance = 1e-12;
% options.PlotFcns = {"optimplotx","optimplotfval","optimplotfunccount"};

[A,fval] = fmincon(@truss,A,[],[],[],[],minArea,[],@stressConstraint,options);

