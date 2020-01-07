clc
clear

x = 0.25;
y = 0.25;

%% Run Abaqus Simulation
command = string(['"C:/Program Files\SIMULIA\Commands\abaqus.bat" cae noGUI=C:\Users\gregj\Documents\GitHub\ME575\eigen\src\runAbaqus.py -- ', num2str(x), " ", num2str(y)]);
command = strjoin(command);
system(command);

%% Post-process Abaqus Simulation
command = '"C:\Program Files\SIMULIA\Commands\abaqus.bat" cae noGUI=C:\Users\gregj\Documents\GitHub\ME575\eigen\src\postAbaqus.py';
system(command);

%% Read in objective function
f = fileread("C:\Users\gregj\Documents\Abaqus\Temp\objectiveFunction.csv");
f = str2double(f);