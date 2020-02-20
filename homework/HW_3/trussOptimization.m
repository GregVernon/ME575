function solverData = trussOptimization(gradMethod)
x0 = ones(10,1);
minArea = 0.1*ones(size(x0));
[mass,stress] = truss(x0);

options = optimoptions('fmincon');
% options.Display = "iter-detailed";
options.MaxFunctionEvaluations = 1e4;
options.MaxIterations = 1e4;
options.OptimalityTolerance = 1e-12;
options.ConstraintTolerance = 1e-12;
options.OutputFcn = @outFun;
% options.PlotFcns = {"optimplotx","optimplotfval","optimplotfunccount"};

solverData.FirstOrderOpt = [];
solverData.FunCount = [];

if strcmpi(gradMethod,"FMINCON")
    % Internal Calculations
    disp("FMINCON")
    options.SpecifyObjectiveGradient = false;
    options.SpecifyConstraintGradient = false;
    [AREA,fval,exitflag,output,lambda] = fmincon(@(A)truss(A),x0,[],[],[],[],minArea,[],@(A)stressConstraint(A),options);
else
    options.SpecifyObjectiveGradient = true;
    options.SpecifyConstraintGradient = true;
    if strcmpi(gradMethod,"FORWARD-DIFFERENCE")
        % Forward Difference
        disp("FORWARD-DIFFERENCE")
        options.CheckGradients = false;
        [AREA,fval,exitflag,output,lambda] = fmincon(@(A)trussWithDerivatives(A,gradMethod),x0,[],[],[],[],minArea,[],@(A)stressConstraintWithDerivatives(A,gradMethod),options);
        
    elseif strcmpi(gradMethod,"CENTERED-DIFFERENCE")
        % Centered Difference
        disp("CENTERED-DIFFERENCE")
        options.CheckGradients = false;
        [AREA,fval,exitflag,output,lambda] = fmincon(@(A)trussWithDerivatives(A,gradMethod),x0,[],[],[],[],minArea,[],@(A)stressConstraintWithDerivatives(A,gradMethod),options);
        
    elseif strcmpi(gradMethod,"COMPLEX-STEP")
        % Complex Step
        disp("COMPLEX-STEP")
        options.CheckGradients = false;
        [AREA,fval,exitflag,output,lambda] = fmincon(@(A)trussWithDerivatives(A,gradMethod),x0,[],[],[],[],minArea,[],@(A)stressConstraintWithDerivatives(A,gradMethod),options);
        
    elseif strcmpi(gradMethod,"AUTO-DIFF")
        % Auto Differentiation
        disp("AUTO-DIFF")
        options.CheckGradients = false;
        [AREA,fval,exitflag,output,lambda] = fmincon(@(A)trussWithDerivatives(A,gradMethod),x0,[],[],[],[],minArea,[],@(A)stressConstraintWithDerivatives(A,gradMethod),options);
        
    elseif strcmpi(gradMethod,"ADJOINT")
        % Adjoint
        disp("ADJOINT")
        options.CheckGradients = false;
        [AREA,fval,exitflag,output,lambda] = fmincon(@(A)trussWithDerivatives(A,gradMethod),x0,[],[],[],[],minArea,[],@(A)stressConstraintWithDerivatives(A,gradMethod),options);
    end
end
solverData.AREA = AREA;
solverData.fval = fval;
solverData.output = output;
solverData.lambda = lambda;

save(gradMethod,"solverData");


    function stop = outFun(x,optimValues,state)
        stop = false;
        
        solverData.FirstOrderOpt = [solverData.FirstOrderOpt  optimValues.firstorderopt];
        solverData.FunCount = [solverData.FunCount optimValues.funccount];
    end
end