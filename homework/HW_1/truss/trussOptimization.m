A = ones(10,1);
minArea = 0.1*ones(size(A));
[mass,stress] = truss(A);

options = optimset;
% options.Display = "iter-detailed";
options.TolFun = 1e-9;
% options.PlotFcns = ["optimplotx","optimplotfval","optimplotfunccount"];

[A,fval] = fmincon(@truss,A,[],[],[],[],minArea,[],@stressConstraint,options);

