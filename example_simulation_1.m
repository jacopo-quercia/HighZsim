%-------- Testbench for Tracking simulator class methods 
close all; 
clear all; 
clc;

% add path for simulator repository
addpath("\\10.79.34.83\Laboratorio\4 - Quercia\SDIClab\matlab\Detector Simulator\OOPv2")

myGeometry = geometryClass(); 

pitch = 500e-6; 
arrayDims = [1, 3];
interpixelGap = 25e-6;
guardDepth = 5e-4;
thickness = 1.85e-3; 
worldOffset = 1e-3;

myGeometry = buildPixelDetector(myGeometry, pitch, arrayDims, interpixelGap, guardDepth, thickness, worldOffset); 


plotVolumes(myGeometry)

mySolver1 = FEMComsolSolverClass();

mySolver1 = createModel(mySolver1,"model1","det","geom1");


mySolver1 = importParams(mySolver1,myGeometry);


mySolver1 = importGeometry(mySolver1);

mySolver1 = ComsolFEMInitialization(mySolver1,'CZT',10.4,'Air',1);

mySolver1 = ComsolComputeWeigthingField(mySolver1,"electrode_pixel12");


figure("Name","Check mesh")
mphmesh(mySolver1.model,"mesh1", "facealpha", 0.5)


figure("Name","Weighting potential of electrode 1")
pg1 = mySolver1.model.result.create('pg1', 'PlotGroup3D');
pg1.feature.create('slice1', 'Slice');
mphplot(mySolver1.model,'pg1','rangenum',1)

mu_e_array = 0.1;
tau_e_array = 2e-5; 

Vbias =  700;


energy = 122.06e3;
counts = 1;

for j = 1 : 3

nBin = 4000;
eMax = 130e3; 
energyThreshold = 1e3; 

ENC = 26.6;
myInstrument = instrumentSimulatorClass(nBin, eMax, ENC); 
ctr = 1; 

accTime = 0; 

figure("Name","energy spectra")

for i = 1 : 25
tic 
% Tracking simulation object 
myTracker = trackingSimulatorClass(); 

% Specify the fraction by weight of each element
myTracker = importMaterialProperties(myTracker, "0.4 Cd 0.1 Zn 0.5 Te", 5.86, 1e3); 

myTracker = importDetectorVolume(myTracker, myGeometry); 

spectrum.energy = energy;
spectrum.counts = counts;

myTracker = defineInitialDistribution(myTracker,"uniform rectangular",[guardDepth/2 (3/2)*guardDepth 0],[1e-3 1e-3 0], spectrum,[0 0]); 

nPhoton = 10000; 

myTracker.energyThreshold = 1e-3; 

myTracker = initializeEEM(myTracker,nPhoton); 
myTracker = particleTracking(myTracker,nPhoton);

myTracker = particleTracking(myTracker,nPhoton);

% figure("Name","Tracking simulation")
% plotVolumes(myGeometry)
% hold on 
% scatter3(myTracker.energyDepositionMap{1,2}.x, myTracker.energyDepositionMap{1,2}.y, myTracker.energyDepositionMap{1,2}.z)
% hold on
% scatter3(myTracker.energyDepositionMap{1,3}.x, myTracker.energyDepositionMap{1,3}.y, myTracker.energyDepositionMap{1,3}.z)
% 
% figure("Name","energy deposition 1")
% hist(myTracker.energyDepositionMap{1,2}.energy,200)

%%%%%%%%%%%%%%%%%%%%%%%%% transport simulation parameter definition %%%%%%%%%%%%%%%%%%%%%%%%%

ehE = 4.46;       % e-h pair mean generation energy, expressed in eV           
mu_e = mu_e_array;      % mobility of electrons, expressed in m^2/(Vs)
mu_h  = 0.0115;     % mobility of electrons, expressed in m^2/(Vs)

tau_e  = tau_e_array;    % mean lifetime of electrons, expressed in s  
tau_h  = 1.5e-6;    % mean lifetime of electrons, expressed in s

T = 298;          % Absolute temperature expressed in K 

epsr = 10.2; 

Fano = 0.11;
generationEnergy = 4.46; 

A = 0.95;             % 10^-6 m/keV
B = 0.98;             % dimensionless
C = 0.003;            % keV^-1


CloudRadiusParameters = [A B C];            % Initial cloud radius parameter (gaussian distribution hypothesis)
timeStep = 2e-9;
eField = [0 0 Vbias/(1.85e-3)]; 
% Construt a transporter obejct 
myTransporter = transportSimulatorClass(ehE, mu_e, tau_e, mu_h, tau_h, T,myTracker,epsr,CloudRadiusParameters, Fano, generationEnergy);
myTransporter = initSimulationParameters(myTransporter,timeStep);
nPseudoCarrier = 100; 
myTransporter = constantElectricFieldStatic(myTransporter, myTracker, mySolver1, nPseudoCarrier, eField); 


myInstrument = recordEnergySpectra(myInstrument, myTransporter,energyThreshold); 

semilogy(myInstrument.energySpectrum.energy,myInstrument.energySpectrum.counts)  
drawnow;  % Refresh plot 

disp("Simulation "+num2str(ctr)+" completed in "+num2str(toc)+" seconds.")
accTime = accTime + toc; 
ctr = ctr +1; 
end 

s = seconds(accTime);
s.Format = 'hh:mm:ss'; 
disp("Simulation ended, total elapsed time is "+char(s)+" min")

% Save simulated spectra in specified folder 
% Specify the directory where you want to create the new folder
% For example, creating the folder in the current working directory:
baseDir = "\\10.79.34.83\Laboratorio\4 - Quercia\Work\11 - Publications\quercia2024montecarlo\simulator validation\Simulation results"; % Change this to any directory you prefer

% Get the current date and time in a specific format
% For example, 'yyyy-mm-dd_HH-MM-SS' format
dateTimeNow = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

% Create the folder name by appending the date and time
folderName = ['SimResults_' dateTimeNow];

% Full path for the new folder
fullFolderPath = fullfile(baseDir, folderName);

% Create the folder
mkdir(fullFolderPath);

% Save the current directory to return later
previousDir = pwd;

% Change the current directory to the new folder
cd(fullFolderPath);

SimulationReport.Tracker = myTracker;
SimulationReport.Transporter = myTransporter;
SimulationReport.Geometry = myGeometry;
SimulationReport.FEMSolver = mySolver1; 
SimulationReport.Instrument = myInstrument; 
save("SimulationReport.mat", "SimulationReport")

energyspectrum(:,1) = myInstrument.energySpectrum.energy./1000; 
energyspectrum(:,2) = myInstrument.energySpectrum.counts;
save("energyspectrum.mat", "energyspectrum")
% Return to the previous directory
cd(previousDir);

end 




