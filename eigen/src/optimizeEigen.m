Opt.nDOF = 4;
Opt.nChildren = 4;
Opt.diffWeight = 0.8;
Opt.crossProb = 0.95;

main(Opt);

%% Helper Functions
function main(Opt)
Parent = initialize(Opt);

for g = 1:100
    disp("Generation: " + num2str(g))
    Child = createNextGeneration(Opt,Parent);
    
    f = zeros(Opt.nChildren,1);
    for c = 1:Opt.nChildren
        disp("Evaluating function: " + num2str(c))
        f(c) = functionEvaluation(Child(c).x);
        Child(c).objVal = f(c);
    end
    
    Parent = tournament(Parent,Child);
end
end

%%%%%%%%%%%

function Parent = tournament(Parent,Child)
for ii = 1:length(Parent)
    if Child(ii).objVal < Parent(ii).objVal
        Parent(ii) = Child(ii);
    end
end
end

%%%%%%%%%%%

function Child = createNextGeneration(Opt,Parent)
nChildren = Opt.nChildren;
nDOF = Opt.nDOF;

parentIDS = 1:nChildren;
for ii = 1:nChildren
    % Pick three distinct parents
    availParents = true(nChildren,1);
    availParents(ii) = false;
    availParents = parentIDS(availParents);
    pID = availParents(randperm(nChildren-1,3));
    
    % Evaluate differential evolution
    A = Parent(pID(1)).x;
    B = Parent(pID(2)).x;
    C = Parent(pID(3)).x;
    for n = 1:nDOF
        Child(ii).x(n,:) = A(n,:) + Opt.diffWeight * (B(n,:) - C(n,:));
    end
    
    % Evaluate crossover
    for n = 1:nDOF
        r = rand();
        if r >= Opt.crossProb
            Child(ii).x(n,:) = Parent(ii).x(n,:);
        end
    end
end
end

%%%%%%%%%%%

function Parent = initialize(Opt)
nChildren = Opt.nChildren;
nDOF = Opt.nDOF;

maxRadius = 1;
rTol = 0.02;

x = cell(nChildren,1);
y = cell(nChildren,1);
for ii = 1:nChildren
    x{ii} = zeros(nDOF,2);
    for n = 1:nDOF
        theta = rand * 2*pi;
        radius = rand * maxRadius;
        x{ii}(n,:) = [radius*cos(theta) radius*sin(theta)];
    end
end

f = zeros(nChildren,1);
for ii = 1:nChildren
    disp("Evaluating function: " + num2str(ii))
    f(ii) = functionEvaluation(x{ii});
end

for ii = 1:nChildren
    Parent(ii).x = x{ii};
    Parent(ii).objVal = f(ii);
end
end

%%%%%%%%%%%

function f = functionEvaluation(x)
%% Write input data file
% fData = fopen("C:\Users\gregj\Documents\Abaqus\Temp\inputData.csv","w+");
fName = "C:\Users\gregj\Documents\Abaqus\Temp\inputData.csv";
csvwrite(fName,x);
%% Create Mesh in Cubit
command = '"C:\Program Files\Cubit 15.4\bin\claro.exe" -nobanner -nographics -nojournal -noecho -information off -batch C:\Users\gregj\Documents\GitHub\ME575\eigen\src\runCubit.py';
system(command);

%% Run Abaqus Simulation
% command = string(['"C:/Program Files\SIMULIA\Commands\abaqus.bat" cae noGUI=C:\Users\gregj\Documents\GitHub\ME575\eigen\src\runAbaqus.py -- ', num2str(x), " ", num2str(y)]);
command = string(['"C:/Program Files\SIMULIA\Commands\abaqus.bat" cae noGUI=C:\Users\gregj\Documents\GitHub\ME575\eigen\src\runAbaqus.py']);
command = strjoin(command);
system(command);

%% Post-process Abaqus Simulation
command = '"C:\Program Files\SIMULIA\Commands\abaqus.bat" cae noGUI=C:\Users\gregj\Documents\GitHub\ME575\eigen\src\postAbaqus.py';
system(command);

%% Read in objective function
f = fileread("C:\Users\gregj\Documents\Abaqus\Temp\objectiveFunction.csv");
f = str2double(f);
f = -f;
end
