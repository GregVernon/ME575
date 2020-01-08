Opt.nDOF = 4;
Opt.nChildren = 20;
% Opt.diffWeight = 0.8;
Opt.crossProb = 0.95;
Opt.InitFile = "";

VID = VideoWriter("drumEigen.mp4",'MPEG-4');
VID.FrameRate = 4;
VID.Quality = 100;
open(VID);

Parent = main(Opt,VID);

%% Helper Functions
for g = 1:1000
function Parent = main(Opt,VID)
    disp("Generation: " + num2str(g))
    Opt.generation = g;
    if g == 1
        if Opt.InitFile == ""
            Parent = initialize(Opt);
            Child = Parent;
        else
            m = matfile(Opt.InitFile);
            Parent = m.Parent;
            Child = Parent;
        end
    else
        Child = createNextGeneration(Opt,Parent);
        parfor c = 1:Opt.nChildren
            disp("Evaluating function: " + num2str(c))
            f = functionEvaluation(Child(c));
            Child(c).objVal = f;
        end
        
        Parent = tournament(Parent,Child);
    end
    writeOutput(Opt,Parent,Child);
    plotParents(Opt,Parent,VID);
end
close(VID)
end

%%%%%%%%%%%

function writeOutput(Opt,Parent, Child)
if Opt.generation == 1
    save("drumEigen_results.mat",'-v7.3',"Child","Opt","Parent");
else
    m = matfile("drumEigen_results.mat",'Writable',true);
    m.Opt = Opt;
    m.Parent = Parent;
    m.Child = Child;
end

for ii = 1:length(Parent)
    if Parent(ii).generation == Opt.generation
        copyfile("drumEigen_"+num2str(ii)+".cae","drumEigen_Parent_"+num2str(ii)+".cae");
        copyfile("drumEigen_"+num2str(ii)+".odb","drumEigen_Parent_"+num2str(ii)+".odb");
    end
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
    Child(ii).generation = Opt.generation;
    Child(ii).individual = ii;
    isValid = false;
    while isValid == false
        F = rand + 0.5;
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
            Child(ii).x(n,:) = A(n,:) + F * (B(n,:) - C(n,:));
        end
        
        % Evaluate crossover
        for n = 1:nDOF
            r = rand();
            if r >= Opt.crossProb
                Child(ii).x(n,:) = Parent(ii).x(n,:);
            end
        end
        
        radius = sqrt(Child(ii).x(:,1).^2 + Child(ii).x(:,2).^2);
        if all(radius <= 0.975)
            isValid = true;
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
    Parent(ii).x = x{ii};
    Parent(ii).objVal = NaN;
    Parent(ii).generation = Opt.generation;
    Parent(ii).individual = ii;
end

parfor ii = 1:nChildren
    disp("Evaluating function: " + num2str(ii))
    f = functionEvaluation(Parent(ii));
    Parent(ii).objVal = f
end
end

%%%%%%%%%%%

function plotParents(Opt,Parent,VID)
X = {Parent.x};
close all
figure
hold on;

t = linspace(0,2*pi,10000);
x = cos(t);
y = sin(t);
plot(x,y);

for ii = 1:length(X)
    scatter(X{ii}(:,1),X{ii}(:,2),60,'.')
end
axis equal
title("Parents @ Generation: " + num2str(Opt.generation))
drawnow

f = getframe(gcf);
writeVideo(VID,f);
end

%%%%%%%%%%%

function f = functionEvaluation(Child)
x = Child.x;
simID = Child.individual;
abqPath = 'C:\Program Files\SIMULIA\Commands\';
cubPath = 'C:\Program Files\Cubit 15.4\bin\';
%% Write input data file
fName = "inputData_" + num2str(simID) + ".csv";
csvwrite(fName,x);
%% Create Mesh in Cubit
command = string(['"' cubPath 'claro.exe" -nobanner -nographics -nojournal -noecho -information off -batch simID=' num2str(simID) ' runCubit.py']);
command = strjoin(command);
[status,~] = system(command);

%% Run Abaqus Simulation
% command = string(['"C:/Program Files\SIMULIA\Commands\abaqus.bat" cae noGUI=C:\Users\gregj\Documents\GitHub\ME575\eigen\src\runAbaqus.py -- ', num2str(x), " ", num2str(y)]);
workerID = labindex;
pause(workerID/10); % Helps avoid .rec file name clashing
command = string(['"' abqPath 'abaqus.bat" cae noGUI=runAbaqus.py -- ', num2str(simID)]);
command = strjoin(command);
[status,~] = system(command);
% pause(1); % Give some time for things to catch up?
%% Post-process Abaqus Simulation
command = string(['"C:\Program Files\SIMULIA\Commands\abaqus.bat" cae noGUI=postAbaqus.py -- ', num2str(simID)]);
command = strjoin(command);
[status,~] = system(command);
%% Read in objective function
success = isfile("objectiveFunction_" + num2str(simID) + ".csv");
if success == true
    f = fileread("objectiveFunction_" + num2str(simID) + ".csv");
    f = str2double(f);
    f = -f;
    delete("objectiveFunction_" + num2str(simID) + ".csv");
else
    f = inf;
end
end
