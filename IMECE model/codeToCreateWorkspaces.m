clc
clear all
close all
addpath(genpath('fluid properties'))

conditionSelection = 'Zanganeh2015'; % options are Zanganeh2015, Beasley1989-low, Beasley1989-high, ...
includePCM = true;

switch conditionSelection
    case 'Zanganeh2015'
        load('inputsZanganeh2015.mat')
        chargeTimeHours = 3;                % 3 is standard\
        chargeTimePartialHour = 1/6;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 5;             % 5 standard
    case 'Beasley2015-low'
        load('inputsBeasley1989.mat')
        chargeTimeHours = 3;                % 3 is standard\
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 0;             % 5 standard
    case 'Beasley2015-high'
        load('inputsBeasley1989.mat')
        chargeTimeHours = 2;                % 3 is standard\
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 0;             % 5 standard
    otherwise
        disp('Invalid condition selection')
        return
end

% inletTemperatureProfile = cell2mat(struct2cell(load('inletTempProfile40Hz3hrcharge75Cmin.mat')));
% useRockTempProfileForPCM = true;

includeConduction = true;
includeRadiation = true;
includeHeatLoss = true;
includeWallLoss = true;
includeCapRadiation = true;

fluidName = 'DryAirAtmo';           % options are CO2HP, CO2LP, DryAir, DryAirAtmo, MoistAir
referenceFluid = 'DryAirAtmo';      % When using CO2HP for fluid at controlled vdot, ref fluid should be DryAir
                                    % for controlled mdot, ref fluid should be same as fluidName
bypassFraction = 0 ;%.15;           % 0.15 used in Zanganeh 2015 but unclear whether it was applied before or after figure plot. Seems to be applied before.

increaseFluidTemps = false;         % Always true (via override later) for CO2, 
increaseBaseTemp = false;           % can be true or false for air

runTransient = true;
hvPCMModifier = 1;

if strcmp(conditionSelection,'Zanganeh2015')
    if includePCM
        plateTempIncreaseCharge = 38;
        plateTempIncreaseDischarge = -7.7;
        wallTempIncreaseCharge = 11.6;
        wallTempIncreaseDischarge = 2.08;
    else
        plateTempIncreaseCharge = 24.6;
        plateTempIncreaseDischarge = 19;
        wallTempIncreaseCharge = 12;
        wallTempIncreaseDischarge = 2;
    end
end

%% Plot selector
plotAverage = false;
plotPCM = false;
transVis= false;

%% ============== Time setup

if strcmp(fluidName,'DryAirAtmo')
    timeStepInverse = 40;                % 41        % samples per second
else
    timeStepInverse = 40;                % 41        % samples per second
end

timeStep = 1/timeStepInverse;                       % seconds per sample
chargeTime = floor((chargeTimeHours+chargeTimePartialHour)*3600*timeStepInverse);        % 3+1/6 hours in seconds
dischargeTime = floor(dischargeTimeHours*3600*timeStepInverse);        % Normally 5h, shortened for quicker runtime while testing charge mode
numCycles = 1;                                      % #*(charge+discharge)
totalTime = floor(numCycles*(chargeTime+dischargeTime)) + 1;
onePercentTime = round(0.01*totalTime);
progressBarPercentTime = onePercentTime;
t = 1;
useMode = 1;

logTime = chargeTime-1 ;%60*60*timeStepInverse; %chargeTime-1; %130*60*timeStepInverse;

%% ============== Tank dimensions

numRows = 40+4;                         % 2 extra rows for inlet/outlet and 2 for endcaps
% tankHeight = 0.2779776;
if includePCM
    heightPCM = tankHeight; %0.09*1; %0.18; % 0.09
    rockHeight = tankHeight - heightPCM;
    numRowsPCM = numRows; %round(4*heightPCM/0.09);
    rowCutoffPCM = numRows; %2+numRowsPCM;
    rowStartRock = rowCutoffPCM+1; 
    
    numRowsRock = 0; %numRows-numRowsPCM-4;
else
    rowStartRock = 1; 
    rockHeight = tankHeight;      
    PCMHeight = 0;
    numRowsPCM = 0;
    rowCutoffPCM = 3;
    numRowsRock = numRows-4;
end

tankRadius = 0.2081784/2;                             % m
tankArea = pi*(tankRadius)^2;               % m2

%% ============== Discretization
rowLengthRock = rockHeight/numRowsRock;
rowVolumeRock = rowLengthRock*pi*(tankRadius)^2;
rowWallSurfaceAreaRock = pi*2*tankRadius*rowLengthRock;

if includePCM
    rowLengthPCM = heightPCM/numRowsPCM;
    rowVolumePCM = rowLengthPCM*pi*(tankRadius)^2;
    rowWallSurfaceAreaPCM = pi*2*tankRadius*rowLengthPCM;
end

%% ============== Insulation parameters
%{
kMicrotherm = (0.026+0.038)/2;                  % W/mK
kFelt = (0.046+0.078)/2;                        % W/mK
kRockwool = 0.038;

tTankWall = 0.003;                              % m
kTankWall = 38;                                 % kplate = 41-35 W/mK from Incropera 2007, 
tMicrotherm(1:3) = [0.02 0 0];
tFelt(1:3) = [0.04 0.06 0.05];
tRockwool(1:3) = [0.1 0.1 0];
            % tank assumed to be made of same material as plate, average value of range used

lossInsulationTermTop = (tankRadius)*(1/kTankWall*log((tankRadius+tTankWall)/(tankRadius))+ ...
    1/kMicrotherm * log((tankRadius+tTankWall+tMicrotherm(1))/(tankRadius+tTankWall)) + ...
    1/kFelt * log((tankRadius+tTankWall+tMicrotherm(1)+tFelt(1))/(tankRadius+tTankWall+tMicrotherm(1))) + ...
    1/kRockwool * log((tankRadius+tTankWall+tMicrotherm(1)+tFelt(1)+tRockwool(1))/(tankRadius+tTankWall+tMicrotherm(1)+tFelt(1))));

lossInsulationTermMiddle = (tankRadius)*(1/kTankWall*log((tankRadius+tTankWall)/(tankRadius))+ ...
    1/kMicrotherm * log((tankRadius+tTankWall+tMicrotherm(2))/(tankRadius+tTankWall)) + ...
    1/kFelt * log((tankRadius+tTankWall+tMicrotherm(2)+tFelt(2))/(tankRadius+tTankWall+tMicrotherm(2))) + ...
    1/kRockwool * log((tankRadius+tTankWall+tMicrotherm(2)+tFelt(2)+tRockwool(2))/(tankRadius+tTankWall+tMicrotherm(2)+tFelt(2))));

lossInsulationTermBottom = (tankRadius)*(1/kTankWall*log((tankRadius+tTankWall)/(tankRadius))+ ...
    1/kMicrotherm * log((tankRadius+tTankWall+tMicrotherm(3))/(tankRadius+tTankWall)) + ...
    1/kFelt * log((tankRadius+tTankWall+tMicrotherm(3)+tFelt(3))/(tankRadius+tTankWall+tMicrotherm(3))) + ...
    1/kRockwool * log((tankRadius+tTankWall+tMicrotherm(3)+tFelt(3)+tRockwool(3))/(tankRadius+tTankWall+tMicrotherm(3)+tFelt(3))));
lossInsulationTerms = [lossInsulationTermTop lossInsulationTermMiddle lossInsulationTermBottom];


% lossInsideTerm = 0.1;
%}

%% ============== Rock properties
%{
densityRocks = 1*2635;                                        % kg/m3

voidFractionRocks = 0.4;                                    % .4
voidFractionRocksInverse = 1-voidFractionRocks;
rockAverageDiameter = 0.032;                                % m
rockAverageRadius = rockAverageDiameter/2;
rockAverageSurfaceArea = 4*pi*(rockAverageRadius)^2;    % m2
rockAverageVolume = 4/3*pi*(rockAverageRadius)^3;       % m3
rocksPerRow = rowVolumeRock*voidFractionRocksInverse / rockAverageVolume;
rowSurfaceAreaRocks = rocksPerRow*rockAverageSurfaceArea;   % m2/row
emissivityRocks = 0.83;
hvRocksCoeff = (1/rockAverageDiameter)^0.92 * 824;
% rockMass = densityRocks * voidFractionRocksInverse*rowVolumeRock*numRowsRock;

% getkRocks = @(x) -0.000000006329*x.^3+1.0794e-5*x.^2-0.007761*x+3.758;
getkRocks = @(x) 5+(1-5)*(x-25)/(700-25);
getRockEnergy = @(temp) 747.0995*temp+0.2838*temp.^2;
calculateRockTemp = @(energy)(-747.0995+sqrt(747.0995^2+4*0.2838*energy))/2/0.2838;
% calculateFluidTemp = @calculateFluidTempSingleRow;

%}

%% ============== Metal Plate properties
%{
voidFractionPlate = 0.142;
voidFractionPlateInverse = 1- voidFractionPlate;
thicknessPlate = 0.02;          % m
kPlate = 38;                    % W/mK 41-35
plateArea = tankArea;
emissivityPlate = 0.7;
densityPlate = 7930;
surfaceAreaPlate = tankArea*(1-voidFractionPlate)*2+221*pi*0.01*0.02;
%}

%% Fluid selection
% fluidPressure = 101.3e3; %12.38e6;                  % 12.38 MPa (critical point is at 7.38 MPa)
switch fluidName
    case 'CO2LP'        
            fluidPressure=12.38e6;
%             disp('Pressure too low to avoid CO2 critical point. Stopping simulation.')
            disp('CO2 properties set to Low Pressure Mode (P=12.38 MPa = 12.38e6 Pa).')
            getEnthalpy = @getEnthalpyCO2;
            getSpecificHeat = @getSpecificHeatCO2;
            getConductivity = @getConductivityCO2;
            getViscosity = @getViscosityCO2;
            getDensity = @getDensityCO2;
            getEnergyDerivative = @getEnergyDerivativeCO2;
        useCO2 = true; 
        increaseTemps = true;
        fluidTempIncreaseAmount = 50;
    case 'CO2HP'
        fluidPressure=25e6;
%             disp('Pressure too low to avoid CO2 critical point. Stopping simulation.')
        disp('CO2 properties set to High Pressure Mode (P=25 MPa = 25e6 Pa).')
        getEnthalpy = @getEnthalpyCO2HP;
        getSpecificHeat = @getSpecificHeatCO2HP;
        getConductivity = @getConductivityCO2HP;
        getViscosity = @getViscosityCO2HP;
        getDensity = @getDensityCO2HP;
        getEnergyDerivative = @getSpecificHeatCO2HP;
        useCO2 = true; 
        increaseTemps = true;
        fluidTempIncreaseAmount = 50;
    case 'DryAir'
        fluidPressure = 25e6;
        getEnthalpy = @getEnthalpyDryAir;
        getSpecificHeat = @getSpecificHeatDryAir;
        getConductivity = @getConductivityDryAir;
        getViscosity = @getViscosityDryAir;
        getDensity = @getDensityDryAir;
        getEnergyDerivative = @getEnergyDerivativeDryAir;
        useCO2 = false;
    case 'DryAirAtmo'
        fluidPressure = 101.3e3;
        getEnthalpy = @getEnthalpyDryAirAtmo;
        getSpecificHeat = @getSpecificHeatDryAirAtmo;
        getConductivity = @getConductivityDryAirAtmo;
        getViscosity = @getViscosityDryAirAtmo;
        getDensity = @getDensityDryAirAtmo;
        getEnergyDerivative = @getSpecificHeatDryAirAtmo;
        useCO2 = false; 
    otherwise
        disp('Invalid fluid. Terminating operation.')
        return
end

switch referenceFluid
    case 'DryAir'
        getReferenceDensity = @getDensityDryAir;
    case 'DryAirAtmo'
        getReferenceDensity = @getDensityDryAirAtmo;
    case 'CO2HP'
        getReferenceDensity = @getDensityCO2HP;        
    otherwise
        disp('Invalid reference fluid')
        return
end
%% PCM/encapsulation Properties

switch conditionSelection
    case 'Zanganeh2015'
        if includePCM
            calculateTemps = @calculateTempsWithPCM;
            tempProfile = @tempProfileNoPCM;
            getmdot = @getmdotWithPCM;
            getStorageTopConduction = @getPCMConduction;
        else
            pcmLowEnergyCutoff = 0;
            pcmHighEnergyCutoff = 0;
            TmeltStart = 0;
            TmeltEnd = 0;
            ceffPCMSolid = 0;
            ceffPCMPhaseChange = 0;
            ceffPCMLiquid = 0;

            calculateTemps = @calculateTempsNoPCM;
            tempProfile = @tempProfileNoPCM;
            getmdot = @getmdotNoPCM;

            getStorageTopConduction = @getRockConduction;
        end
    otherwise
            calculateTemps = @calculateTempsNoRocks;
            tempProfile = @tempProfileBeasley;
            getmdot = @getmdotBeasley;
            getStorageTopConduction = @getPCMConduction;
end

% if includePCM
%     thicknessEncapsulation = 0; %0.001; % m
%     ST = 0.0253;
%     SD = 0.0253;
%     emissivityEncapsulation = 0.7;
%     surfaceAreaPCM = 103.45;                % m2/m3
%     heightPCM = 0.09;
%     heightPCMRow = heightPCM/4;
%     voidFractionPCM = 0.369; %0.549;
%     voidFractionPCMInverse = 1 - voidFractionPCM;

%     dEncap = 0.02; % 0.018;
%     rEncap = dEncap/2;
%     pcmAverageSurfaceArea = 4*pi*(rEncap)^2;    % m2
%     pcmAverageVolume = 4/3*pi*(rEncap)^3;       % m3
%     pcmPerRow = rowVolumePCM*voidFractionPCMInverse / pcmAverageVolume;
%     rowSurfaceAreaPCM = pcmPerRow*pcmAverageSurfaceArea;   % m2/row

%     Crow = 0.95;
%     densityPCM = 811.81572; % 2650;                      % kg/m3
%     densityEncap = 7930;                    % kg/m3
%     enthalpyFusion = 197.71*1000; %466*1000;                   % J/kg
%     TmeltStart = 47.13; %573 ; %+ fluidTempIncreaseAmount;
%     dTMelt = 5; %4;                             % K
%     TmeltEnd = TmeltStart + dTMelt;
    
%     massPCMRow=2.3975;
%     massEncapsRow = 0; %3.2925;
%     volumePCMRow = 1.23825e-3;
%     cPCMSolid = 8.3736; %1070;                   % J/kgK
%     cPCMPhaseChange = enthalpyFusion/dTMelt;
%     cPCMLiquid = 2.0934;%1170;
%     kPCMSolid = 160;                    % W/mK
    
%     cEncapsSolid = 535;                 % J/kgK
%     cEncapsPhaseChange = 581;
%     cEncapsLiquid = 590;
%     kEncapsSolid = 0.24; %17.7;                % W/mK
%     kEncapsPhaseChange = 0.24; %22;
%     kEncapsLiquid = 0.24; %23.1;
    
%     fcontPCM = 0.05;
%     fcontWall = 0.02;
%     rowSurfaceAreaPCM = rowVolumePCM/surfaceAreaPCM;
    
%     ceffPCMSolid = cPCMSolid; %(cPCMSolid*massPCMRow + cEncapsSolid*massEncapsRow)/(massPCMRow+massEncapsRow); % J/kgK
%     ceffPCMPhaseChange = cPCMPhaseChange; %(cPCMPhaseChange*massPCMRow + cEncapsPhaseChange*massEncapsRow)/(massPCMRow+massEncapsRow); % J/kgK
%     ceffPCMLiquid = cPCMLiquid; %(cPCMLiquid*massPCMRow + cEncapsLiquid*massEncapsRow)/(massPCMRow+massEncapsRow); % J/kgK
    
%     pcmLowEnergyCutoff = TmeltStart*ceffPCMSolid;
%     pcmHighEnergyCutoff = pcmLowEnergyCutoff + ceffPCMPhaseChange*dTMelt;
    
%     densityPCMEff = densityPCM; % (massPCMRow + massEncapsRow)/volumePCMRow;  
%     calculateTemps = @calculateTempsNoRocks; %@calculateTempsWithPCM;
%     tempProfile = @tempProfileWithPCM;
%     getmdot = @getmdotWithPCM;
%     tempProfile = @tempProfileBeasley; %@tempProfileNoPCM;
%     getmdot = @getmdotBeasley; %@getmdotNoPCM;
    
%     getStorageTopConduction = @getPCMConduction;
% else
%     pcmLowEnergyCutoff = 0;
%     pcmHighEnergyCutoff = 0;
%     TmeltStart = 0;
%     TmeltEnd = 0;
%     ceffPCMSolid = 0;
%     ceffPCMPhaseChange = 0;
%     ceffPCMLiquid = 0;
%     
%     calculateTemps = @calculateTempsNoPCM;
%     tempProfile = @tempProfileNoPCM;
%     getmdot = @getmdotNoPCM;
%     
%     getStorageTopConduction = @getRockConduction;
% end


%% ============= Temperature conditions
if increaseFluidTemps || useCO2               % Set tempIncreaseAmount
    fluidTempIncreaseAmount = 50;        % How much (degC) base and air temps are increased 
else
    fluidTempIncreaseAmount = 0;
end

if increaseBaseTemp
    baseTempIncreaseAmount = 50;
else
    baseTempIncreaseAmount = 0;
end

switch conditionSelection
    case 'Zanganeh2015'
        dischargeTemp = 25+baseTempIncreaseAmount;                                 % C - temperature of HTF at bottom inlet during discharge
        baseTemp = 25+baseTempIncreaseAmount;                                      % C - starting temperature of rock bed
        ambientTemp = 25;
    otherwise
        baseTemp = 26.5;
        ambientTemp = 25;
end

%% ============== Temp/energy/step Matrix setup

tempMatrixStorage = ones(numRows,totalTime);    % C
tempMatrixStorage(:,1) = baseTemp;    % C
intEnergyMatrixStorage = zeros(numRows,2);

% rock internal energy, J/kg
% intEnergyRock = 747.0995*tempRock+0.2838*tempRock.^2;       % J
if includePCM
    intEnergyMatrixStorage(1:rowCutoffPCM,1) = ceffPCMSolid * tempMatrixStorage(1:rowCutoffPCM,1);
    intEnergyMatrixStorage(rowStartRock:numRows,1) = getRockEnergy(tempMatrixStorage(rowStartRock:numRows,1));
    totalEnergyRocksStart = sum(intEnergyMatrixStorage(rowStartRock:numRows-2,1));
    totalEnergyPCMStart = sum(intEnergyMatrixStorage(3:rowCutoffPCM,1));
else
    intEnergyMatrixStorage(1:numRows,1) = getRockEnergy(tempMatrixStorage(1:numRows,1));
    totalEnergyRocksStart = sum(intEnergyMatrixStorage(3:numRows-2,1));
end

tempMatrixFluid = ones(numRows,totalTime);    % C
tempMatrixFluid(:,1) = baseTemp;    % C
intEnergyMatrixFluid = zeros(numRows,2);

intEnergyMatrixFluid(1:numRows,1) = getEnthalpy(tempMatrixFluid(1:numRows,1));

intEnergyMatrixStorage(:,2) = intEnergyMatrixStorage(:,1);
intEnergyMatrixFluid(:,2) = intEnergyMatrixFluid(:,1);

%% ============== Fluid Base properties and flow rates
% densityFluidDischarge = getDensity(dischargeTemp);
    densityFluid(numRows) = 0;
    viscosityFluid(numRows) = 0;
    kFluid(numRows) = 0;
    specificHeatFluid(numRows) = 0;

for i = 1:numRows
    densityFluid(i) = getDensity(tempMatrixFluid(i,t));
    viscosityFluid(i) = getViscosity(tempMatrixFluid(i,t));
    kFluid(i) = getConductivity(tempMatrixFluid(i,t));
    specificHeatFluid(i) = getSpecificHeat(tempMatrixFluid(i,t));
end

mdot = zeros(1,totalTime);
mdot2 = zeros(1,totalTime);
vdot = zeros(1,totalTime);
bypassFractionInverse = 1-bypassFraction;               % bypass fraction 
mdot(1,1)=0.008*bypassFractionInverse;

G = mdot(1,1)/tankArea;

%% ============== Energy and heat transfer variables

fluidStepConvection = 0;
fluidStepAdvection = 0;
fluidStepHeatLoss = 0;

storageStepConvection = 0;
storageStepConduction = 0;
storageStepRadiation = 0;

% getFluidConvection = @(hveff,fluiddensity,ts,tf) timeStep*hveff/voidFractionRocks/fluiddensity*(ts-tf);
% getFluidAdvection = @(mdot,fluiddensity,hin,hout) timeStep*mdot/fluiddensity/rowVolumeRock/voidFractionRocks*(hin-hout);
% getFluidHeatLoss = @(uwall,fluiddensity,tf,tinf) timeStep*uwall*rowWallSurfaceAreaRock*(tf-tinf)/voidFractionRocks/fluiddensity/rowVolumeRock;
% 
% getRockConvection = @(hveff,tf,ts) timeStep*hveff/fluidVolumeFraction/densityRocks*(tf-ts);
% getRockConduction = @(keff,tsprev,ts,tsnext) timeStep*rowSurfaceAreaRocks*keff/fluidVolumeFraction/densityRocks/rowLengthRock/rowVolumeRock*(tsprev+tsnext-2*ts);
% getRockRadiation = @(ts) 0*ts;

energyInput = 0;                           
energyRadiatedPlateCharging = 0;
energyRadiatedPlateDischarging = 0;
energyWallLossPCMCharging = 0;
energyWallLossPCMDischarging = 0;
energyWallLossRocksCharging = 0;           
energyWallLossRocksDischarging = 0;
energyStoredRocks = 0;                     
energyStoredPCM = 0;


%% Loggers
pressureDrop(totalTime) = 0;
logEnergyFlow(numRows,13) = 0;
logEnergyFlowJ(numRows,13) = 0;
energyIn = 0;
logConvectionStorage(totalTime) = 0;
logConduction(totalTime) = 0;
logConvectionFluid(totalTime) = 0;
logRadiationPlate(totalTime) = 0;
logRockLoss(totalTime) = 0;
logPCMLoss(totalTime) = 0; 
logEnergyIn(totalTime) = 0;
logEnergyOut(totalTime) = 0;
logEnergyStored(totalTime) = 0;
logQRadPlate(totalTime) = 0;
logQWallPCM(totalTime) = 0;

% save('inputsBeasley1989fsafasf.mat')

if ~runTransient
    disp('Setup complete. Transient analysis not performed.')
    return
end