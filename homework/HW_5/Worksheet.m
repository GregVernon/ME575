clear
close all

n = 2:20;
time = nan(7,length(n));
err = nan(7,length(n));
fCount = nan(7,length(n));

for ii = 1:length(n)
    N = n(ii)
    x0 = 20*(rand(N,1)-0.5);
    LB = -10 * ones(size(x0));
    UB =  10 * ones(size(x0));
    
    objFun = @(x) rosenbrock(x);
    [f,df] = objFun(rand(N)); % Warmup
    exactSolution = ones(N,1);

%% Gradient-Based
%%Exact Gradient
    options = optimoptions('fminunc');
    options.MaxFunctionEvaluations = 1e6;
    options.MaxIterations = 1e6;
    options.SpecifyObjectiveGradient = true;
    options.Display = "off";
    tic;
    [sol,~,~,output] = fminunc(objFun,x0,options);
    t_elapsed = toc;
    err(1,ii) = norm(exactSolution - sol);
    time(1,ii) = t_elapsed;
    fCount(1,ii) = output.funcCount;
%%Forward Difference
    options = optimoptions('fminunc');
    options.MaxFunctionEvaluations = 1e6;
    options.MaxIterations = 1e6;
    options.FiniteDifferenceType = "forward";
    options.OptimalityTolerance = 1e-12;
    options.Display = "off";
    tic;
    [sol,~,~,output] = fminunc(objFun,x0,options);
    t_elapsed = toc;
    err(2,ii) = norm(exactSolution - sol);
    time(2,ii) = t_elapsed;
    fCount(2,ii) = output.funcCount;
%%Centered Difference
    options = optimoptions('fminunc');
    options.MaxFunctionEvaluations = 1e6;
    options.MaxIterations = 1e6;
    options.FiniteDifferenceType = "central";
    options.OptimalityTolerance = 1e-12;
    options.Display = "off";
    tic;
    [sol,~,~,output] = fminunc(objFun,x0,options);
    t_elapsed = toc;
    err(3,ii) = norm(exactSolution - sol);
    time(3,ii) = t_elapsed;
    fCount(3,ii) = output.funcCount;

%% Gradient Free
%%Simulated Annealing
    options = optimoptions("simulannealbnd");
    options.MaxFunctionEvaluations = 1e30;
    options.MaxIterations = 1e30;
    options.Display = "off";
    tic;
    [sol,~,~,output] = simulannealbnd(objFun, x0, LB, UB, options);
    t_elapsed = toc;
    err(4,ii) = norm(exactSolution - sol);
    time(4,ii) = t_elapsed;
    fCount(4,ii) = output.funccount;
%%Particle Swarm
    options = optimoptions("particleswarm");
    options.MaxIterations = 1e30;
    options.Display = "off";
    tic;
    [sol,~,~,output] = particleswarm(objFun,N,LB,UB,options);
    t_elapsed = toc;
    err(5,ii) = norm(exactSolution - sol);
    time(5,ii) = t_elapsed;
    fCount(5,ii) = output.funccount;
%%Genetic Algorithm
    options = optimoptions("ga");
    options.MaxGenerations = 1e30;
    options.Display = "off";
    tic;
    [sol,~,~,output] = ga(objFun,N,[],[],[],[],LB,UB,[],[],options);
    t_elapsed = toc;
    err(6,ii) = norm(exactSolution - sol);
    time(6,ii) = t_elapsed;
    fCount(6,ii) = output.funccount;
%%Pattern Search
    options = optimoptions("patternsearch");
    options.MaxFunctionEvaluations = 1e30;
    options.MaxIterations = 1e30;
    options.Display = "off";
    tic;
    [sol,~,~,output] = patternsearch(objFun,x0,[],[],[],[],[],[],[],options);
    t_elapsed = toc;
    err(7,ii) = norm(exactSolution - sol);
    time(7,ii) = t_elapsed;
    fCount(7,ii) = output.funccount;

end

legendNames = ["Exact Gradient", "Forward Difference", "Centered Difference", "Simulated Annealing", "Particle Swarm", "Genetic Algorithm", "Pattern Search"];

figure();  % Plot Time per method vs nVar
plot(n, time)
legend(legendNames)
xlabel("# Variables")
ylabel("Elapsed Time (s)")
ax = gca;
ax.YScale = "log";

figure();  % Plot Err per method vs nVar
plot(n, err)
legend(legendNames)
xlabel("# Variables")
ylabel("Error Norm")
ax = gca;
ax.YScale = "log";

figure();  % Plot nFuncEval per method vs nVar
plot(n,fCount)
legend(legendNames)
xlabel("# Variables")
ylabel("# Function Evaluations")
ax = gca;
ax.YScale = "log";



%% Functions
function [f, df] = rosenbrock(x)
% Compute function value
f = 0;
for ii = 1:length(x)-1
    f = f + 100*(x(ii+1) - x(ii)^2)^2 + (1-x(ii))^2;
end

if nargout > 1 % gradient required
    N = length(x);
    df = zeros(N,1);
    for ii = 1:N
        if ii == 1
            df(ii) = 2*x(ii) - 400*x(ii)*(x(ii+1) - x(ii)^2) - 2;
        elseif ii == N
            df(ii) = 200*x(ii) - 200*x(ii-1)^2;
        else
            df(ii) = -200*x(ii-1)^2 + 202*x(ii) - 400*x(ii)*(x(ii+1)-x(ii)^2) - 2;
        end
    end
end
end
