clc
clear all
close all
addpath(genpath('fluid properties'))


conditionSelection = 'Zanganeh2015'; % options are Zanganeh2015, Beasley1989-low, Beasley1989-high, ...

includePCM = false;

includeConduction = true;
includeRadiation = true;
includeHeatLoss = true;
includeWallLoss = true;
includeCapRadiation = true;
useInternalEnergy = false;

fluidName = 'DryAirAtmo';           % options are CO2HP, CO2LP, DryAir, DryAirAtmo, MoistAir
referenceFluid = 'DryAirAtmo';      % When using CO2HP for fluid at controlled vdot, ref fluid should be DryAir
                                    % for controlled mdot, ref fluid should be same as fluidName
bypassFraction = 0 ;%.15;           % 0.15 used in Zanganeh 2015 but unclear whether it was applied before or after figure plot. Seems to be applied before.

increaseFluidTemps = false;         % Always true (via override later) for CO2, 
increaseBaseTemp = false;           % can be true or false for air

runTransient = true;
hvPCMModifier = 1;

if strcmp(fluidName,'CO2HP') || strcmp(fluidName,'CO2LP')
    load('inletTempProfile40Hz3hrcharge75Cmin.mat')
end
useRockTempProfileForPCM = false;

%% Plot selector
plotAverage = true;
plotPCM = false;
transVis= false;

%% Condition set
switch conditionSelection
    case 'Zanganeh2015'
        load('inputsZanganeh2015.mat')
        chargeTimeHours = 3;                % 3 is standard\
        chargeTimePartialHour = 1/6;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 5;             % 5 standard
        
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
    case 'Beasley2015-low'
        load('inputsBeasley1989.mat')
        chargeTimeHours = 3;                % 3 is standard\
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 0;             % 5 standard
        includePCM = true;
    case 'Beasley2015-high'
        load('inputsBeasley1989.mat')
        chargeTimeHours = 2;                % 3 is standard\
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 0;             % 5 standard
        includePCM = true;
    otherwise
        disp('Invalid condition selection')
        return
end

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

logTime = 78706; %60*60*timeStepInverse; %chargeTime-1; %130*60*timeStepInverse;
breakAfterLogTime = false;

%% ============== Tank dimensions
% numRows = 40+4;                         % 2 extra rows for inlet/outlet and 2 for endcaps

% tankHeight = 0.2779776;
if includePCM
    heightPCM = 0.09*1; %0.18; % 0.09
    rockHeight = tankHeight - heightPCM;
    numRowsPCM = round(4*heightPCM/0.09);
    rowCutoffPCM = 2+numRowsPCM;
    rowStartRock = rowCutoffPCM+1; 
    
    numRowsRock = numRows-numRowsPCM-4;
else
    rowStartRock = 1; 
    rockHeight = tankHeight;      
    PCMHeight = 0;
    numRowsPCM = 0;
    rowCutoffPCM = 3;
    numRowsRock = numRows-4;
end

% tankRadius = 0.2081784/2;                             % m
% tankArea = pi*(tankRadius)^2;               % m2

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
% {
% densityRocks = 1*2635;                                        % kg/m3

% voidFractionRocks = 0.4;                                    % .4
% voidFractionRocksInverse = 1-voidFractionRocks;
% rockAverageDiameter = 0.032;                                % m
% rockAverageRadius = rockAverageDiameter/2;
% rockAverageSurfaceArea = 4*pi*(rockAverageRadius)^2;    % m2
% rockAverageVolume = 4/3*pi*(rockAverageRadius)^3;       % m3
rocksPerRow = rowVolumeRock*voidFractionRocksInverse / rockAverageVolume;
rowSurfaceAreaRocks = rocksPerRow*rockAverageSurfaceArea;   % m2/row
% emissivityRocks = 0.83;
% hvRocksCoeff = (1/rockAverageDiameter)^0.92 * 824;
% rockMass = densityRocks * voidFractionRocksInverse*rowVolumeRock*numRowsRock;

% getkRocks = @(x) -0.000000006329*x.^3+1.0794e-5*x.^2-0.007761*x+3.758;
% getkRocks = @(x) 5+(1-5)*(x-25)/(700-25);
% getRockEnergy = @(temp) 747.0995*temp+0.2838*temp.^2;
% calculateRockTemp = @(energy)(-747.0995+sqrt(747.0995^2+4*0.2838*energy))/2/0.2838;
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

%% Transient
tic
disp('Starting transient analysis...')
startTime = clock;
fprintf('Time: %g:%g:%g\n',startTime(4),startTime(5),floor(startTime(6)))

while t < totalTime
    
    intEnergyMatrixStorage(:,1) = intEnergyMatrixStorage(:,2);
    intEnergyMatrixFluid(:,1) = intEnergyMatrixFluid(:,2);
    
    % Fluid properties
    for i = 1:numRows
        densityFluid(i) = getDensity(tempMatrixFluid(i,t));
        viscosityFluid(i) = getViscosity(tempMatrixFluid(i,t));
        kFluid(i) = getConductivity(tempMatrixFluid(i,t));
        specificHeatFluid(i) = getSpecificHeat(tempMatrixFluid(i,t));
    end
    
    % Storage conductivity
    if includePCM
        % [rowCutoffPCM rowStartRock]
        kRocks (rowStartRock:numRows) = getkRocks(tempMatrixStorage(rowStartRock:numRows,t));
        for i=1:1:rowCutoffPCM
            if tempMatrixStorage(i,t) < TmeltStart
                kRocks(i) = 17.7;
            elseif tempMatrixStorage(i,t) < TmeltEnd
                kRocks(i) = 22;
            else
                kRocks(i) = 23.1;
            end
        end
    else
        kRocks = getkRocks(tempMatrixStorage(1:numRows,t));
    end
    
    % Dispersion and loss coefficients
   [keffRocks,heffWall] = calculatekeff(kRocks,kFluid,tempMatrixFluid(:,t),tempMatrixStorage(:,t),voidFractionRocks,rockAverageDiameter,numRows,emissivityRocks);
    
   % Switch to discharge mode, log energy before proceeding
    if mod(t,(chargeTime+dischargeTime)) == chargeTime+1
        useMode = 2;
        if includePCM
            totalEnergyRocksEndCharge = sum(intEnergyMatrixStorage(rowStartRock:numRows-2,1));
            storageEnergyMatrixChargeEndLog = intEnergyMatrixStorage(:,1);
            totalEnergyPCMEndCharge = sum(intEnergyMatrixStorage(3:rowCutoffPCM,1));
        else
            totalEnergyRocksEndCharge = sum(intEnergyMatrixStorage(3:numRows-2,1));
            storageEnergyMatrixChargeEndLog = intEnergyMatrixStorage(:,1);
        end
        % logInternalEnergy = intEnergyMatrix(:,2);

    % Switch back to charge mode for extra cycles
    elseif mod(t,(chargeTime+dischargeTime)) == 1
        useMode = 1;
    end
    
    vdot(t) = volumeFlowRate(t,chargeTime,timeStepInverse,bypassFractionInverse,dischargeTemp,getmdot,getReferenceDensity);
%     vdot(t) = volumeFlowRate(t,chargeTime,timeStepInverse,bypassFractionInverse,25,getmdot,getReferenceDensity);
%     mdot(t) = getmdot(t,chargeTime,timeStepInverse,bypassFractionInverse);
    mdot(t) = vdot(t) * getDensity(dischargeTemp);
    G=mdot(t)/tankArea;
    
    hvRocks = hvRocksCoeff*G^0.92;  % W/m2K
    hpRocks = hvRocks*rockAverageDiameter/6/(voidFractionRocksInverse);
    
    if useMode == 1 
        if strcmp(fluidName,'DryAirAtmo')
        tempMatrixFluid(2,t+1) = tempProfile(t,timeStep,0,tempMatrixFluid(2,t),includePCM,1);
%         tempMatrixFluid(2,t+1) = tempProfile(t,timeStep,0,tempMatrixFluid(2,t),false,1);
        else
            % Inlet temperature profile from PCM+ rocks case, increased so it's never
            % below 75 (for sCO2 purposes)
            tempMatrixFluid(2,t+1) = inletTempProfile(1,t+1);
        end
        
        tempMatrixFluid(1,t+1) = tempMatrixFluid(2,t+1);
        
        intEnergyMatrixFluid(1:2,1) = getEnthalpy(tempMatrixFluid(1:2,t));
        
        for i=1:1:numRows
            BiRocks = hpRocks * rockAverageRadius/kRocks(i);
            hvEffRocks = hvRocks/(1+0.25*BiRocks);
            
            switch i
                case 1 % top plate top
                    fluidStepConvection = 0;
                    fluidStepAdvection = 0;
                    fluidStepHeatLoss =0; 
                    storageStepConvection = 0;
                    storageStepConduction = 0;
                    storageStepRadiation = 0;
                    storageStepWall = 0;
                    storageStepPlate = 0;
                case 2 % Tinlettop / top plate bottom
                    fluidStepConvection = 0;
                    fluidStepAdvection = 0;
                    fluidStepHeatLoss =0; 
                    storageStepConvection = 0;
                    storageStepConduction = 0;
                    storageStepRadiation = 0;
                    storageStepWall = 0;
                    storageStepPlate = 0;
                case numRows-1 %bottom plate top
                    fluidStepConvection = 0;
                    fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
                    fluidStepHeatLoss = 0; %getFluidHeatLoss(0,lossInsulationTerms,densityFluid(i),tempMatrixFluid(i,t),ambientTemp,timeStep,rowWallSurfaceAreaRock,voidFractionRocks,rowVolumeRock,i);
                    storageStepConvection = 0;
                    storageStepConduction = 0;
                    storageStepRadiation = 0;
                    storageStepWall = 0;
                    storageStepPlate = 0;
                case numRows %Toutletbottom / bottom plate top
                    fluidStepConvection = 0;
                    fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
                    fluidStepHeatLoss = 0; %getFluidHeatLoss(0,lossInsulationTerms,densityFluid(i),tempMatrixFluid(i,t),ambientTemp,timeStep,rowWallSurfaceAreaRock,voidFractionRocks,rowVolumeRock,i);
                    storageStepConvection = 0;
                    storageStepConduction = 0;
                    storageStepRadiation = 0;
                    storageStepWall = 0;
                    storageStepPlate = 0;
                case 3              % PCM or rocks
                    if includePCM
                        Re = 3.4658*G*dEncap/viscosityFluid(3);
                        % Based on umax not u0, thus increased by
                        % umax=ST/(SD-dEncap)*u0 = 3.4658*u0
                        
                        Pr = specificHeatFluid(3)*viscosityFluid(3)/kFluid(3);
                        PrS = getSpecificHeat(tempMatrixStorage(3))*getViscosity(tempMatrixStorage(3))/getConductivity(tempMatrixStorage(3));

                        NuPCM = 0.4845 * sqrt(Re)*Pr^0.37*(Pr/PrS)^0.25;
                                                
                        hpPCM = hvPCMModifier*NuPCM*kFluid(3)/dEncap;
                        hvPCM = hpPCM*surfaceAreaPCM;
                        BiPCM = hpPCM * rEncap / kPCMSolid;

                        hvEffPCM = hvPCM/(1+0.25*BiPCM);

                        fluidStepConvection = getFluidConvection(hvEffPCM,densityFluid(3),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionPCM);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(3),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(3),timeStep,rowVolumePCM,voidFractionPCM);
                        fluidStepHeatLoss = 0;% getFluidHeatLoss(includeHeatLoss,i,heffWall(3),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(3),kFluid(3),Re,Pr,voidFractionPCM,dEncap,rowVolumePCM,rowWallSurfaceAreaPCM,timeStep);

                        storageStepConvection = getRockConvection(hvEffPCM,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionPCMInverse,densityPCMEff);
                        storageStepConduction = getPCMConduction(includeConduction,tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,tankArea,fcontPCM,voidFractionPCMInverse,densityPCMEff,rowVolumePCM);
                        storageStepRadiation = getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixStorage(2,t),tankArea,voidFractionPlateInverse,emissivityEncapsulation,emissivityPlate,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeRadiation)+ ...
                            getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixStorage(4,t),tankArea,voidFractionPCMInverse,emissivityEncapsulation,emissivityEncapsulation,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeRadiation);
                        storageStepWall = getQWallPCM(tempMatrixStorage(3,t),min(tempMatrixFluid(3,t),tempMatrixStorage(3,t)+wallTempIncreaseCharge),rowWallSurfaceAreaPCM,fcontWall,tankRadius,thicknessEncapsulation,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,timeStep,includeWallLoss);
%                         storageStepPlate = getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixStorage(3,t)+plateTempIncreaseCharge,plateArea,voidFractionPlateInverse,emissivityEncapsulation,emissivityPlate,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeCapRadiation) ;
                        storageStepPlate = getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixFluid(2,t)+plateTempIncreaseCharge,plateArea,voidFractionPlateInverse,emissivityEncapsulation,emissivityPlate,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeCapRadiation) ;
                        
                        logQWallPCM(t) = logQWallPCM(t) + storageStepWall*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logQRadPlate(t) = logQRadPlate(t) + storageStepPlate*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionPCM*densityFluid(3)*rowVolumePCM;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;

                        logPCMLoss(t) = logPCMLoss(t) + fluidStepHeatLoss*voidFractionPCM*densityPCMEff*rowVolumePCM; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionPCM*rowVolumePCM*densityFluid(3);
                    else
                        Re = G*rockAverageDiameter/viscosityFluid(3);
                        Pr = specificHeatFluid(3)*viscosityFluid(3)/kFluid(3);

                        fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(3),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(3),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(3),timeStep,rowVolumeRock,voidFractionRocks);

                        fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(3),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(3),kFluid(3),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                        storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                        storageStepConduction = getRockConduction(includeConduction,keffRocks(i-1),keffRocks(3),tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);

                        storageStepRadiation =0;% getRockRadiation(tempMatrixStorage(i,t));
                        storageStepWall = 0;
                        storageStepPlate = getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixFluid(2,t)+plateTempIncreaseCharge,plateArea,voidFractionPlateInverse,emissivityRocks,emissivityPlate,timeStep,voidFractionRocksInverse,densityRocks,rowVolumeRock,includeCapRadiation) ;

    %                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
    %                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(3),G);
                        
                        logQRadPlate(t) = logQRadPlate(t) + storageStepPlate*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(3)*rowVolumeRock;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(3)*rowVolumeRock; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(3);
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + logRockLoss(t);
                    end

                case num2cell(4:rowCutoffPCM-1)              % PCM or rocks
                    if includePCM
                        Re = 3.4658*G*dEncap/viscosityFluid(i);
                        % Based on umax not u0, thus increased by
                        % umax=ST/(SD-dEncap)*u0 = 3.4658*u0
                        
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                        PrS = getSpecificHeat(tempMatrixStorage(i))*getViscosity(tempMatrixStorage(i))/getConductivity(tempMatrixStorage(i));

                        NuPCM = 0.4845 * sqrt(Re)*Pr^0.37*(Pr/PrS)^0.25;
                                                
                        hpPCM = hvPCMModifier*NuPCM*kFluid(3)/dEncap;
                        hvPCM = hpPCM*surfaceAreaPCM;
                        BiPCM = hpPCM * rEncap / kPCMSolid;

                        hvEffPCM = hvPCM/(1+0.25*BiPCM);

                        fluidStepConvection = getFluidConvection(hvEffPCM,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionPCM);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(i),timeStep,rowVolumePCM,voidFractionPCM);
                        fluidStepHeatLoss = 0; %getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionPCM,dEncap,rowVolumePCM,rowWallSurfaceAreaPCM,timeStep);

                        storageStepConvection = getRockConvection(hvEffPCM,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionPCMInverse,densityPCMEff);
                        storageStepConduction = getPCMConduction(includeConduction,tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,tankArea,fcontPCM,voidFractionPCMInverse,densityPCMEff,rowVolumePCM);
                        storageStepRadiation =getPCMRadiation(tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),tankArea,densityPCMEff,rowVolumePCM,emissivityEncapsulation,timeStep,includeRadiation);
                        storageStepWall = getQWallPCM(tempMatrixStorage(i,t),min(tempMatrixFluid(i,t),tempMatrixStorage(i,t)+wallTempIncreaseCharge),rowWallSurfaceAreaPCM,fcontWall,tankRadius,thicknessEncapsulation,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,timeStep,includeWallLoss);
                        storageStepPlate = 0;
                        
                        logQWallPCM(t) = logQWallPCM(t) + storageStepWall*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionPCM*densityFluid(i)*rowVolumePCM;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;

                        logPCMLoss(t) = logPCMLoss(t) + fluidStepHeatLoss*voidFractionPCM*densityPCMEff*rowVolumePCM; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionPCM*rowVolumePCM*densityFluid(i);
                    else
                        Re = G*rockAverageDiameter/viscosityFluid(i);
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);

                        fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);

                        fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                        storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                        storageStepConduction = getRockConduction(includeConduction,keffRocks(i-1),keffRocks(i),tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);

                        storageStepRadiation =0;% getRockRadiation(tempMatrixStorage(i,t));
                        storageStepWall = 0;
                        storageStepPlate = 0;

    %                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
    %                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
                        
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + logRockLoss(t);
                    end
                case rowCutoffPCM              % PCM or rocks
                    if includePCM
                        Re = 3.4658*G*dEncap/viscosityFluid(i);
                        % Based on umax not u0, thus increased by
                        % umax=ST/(SD-dEncap)*u0 = 3.4658*u0
                        
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                        PrS = getSpecificHeat(tempMatrixStorage(i))*getViscosity(tempMatrixStorage(i))/getConductivity(tempMatrixStorage(i));

                        NuPCM = 0.4845 * sqrt(Re)*Pr^0.37*(Pr/PrS)^0.25;
                                                
                        hpPCM = hvPCMModifier*NuPCM*kFluid(3)/dEncap;
                        hvPCM = hpPCM*surfaceAreaPCM;
                        BiPCM = hpPCM * rEncap / kPCMSolid;

                        hvEffPCM = hvPCM/(1+0.25*BiPCM);

                        fluidStepConvection = getFluidConvection(hvEffPCM,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionPCM);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(i),timeStep,rowVolumePCM,voidFractionPCM);
                        fluidStepHeatLoss = 0; %getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionPCM,dEncap,rowVolumePCM,rowWallSurfaceAreaPCM,timeStep);

                        storageStepConvection = getRockConvection(hvEffPCM,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionPCMInverse,densityPCMEff);
                        storageStepConduction = getPCMConduction(includeConduction,tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i,t),timeStep,tankArea,fcontPCM,voidFractionPCMInverse,densityPCMEff,rowVolumePCM);
                        storageStepRadiation = getInterfaceRadiation(tempMatrixStorage(i,t),tempMatrixStorage(i-1,t),tankArea,voidFractionPCMInverse,emissivityEncapsulation,emissivityEncapsulation,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeRadiation)+ ...
                            getInterfaceRadiation(tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),tankArea,voidFractionPCMInverse,emissivityEncapsulation,emissivityRocks,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeRadiation);
                        storageStepWall = getQWallPCM(tempMatrixStorage(i,t),min(tempMatrixFluid(i,t),tempMatrixStorage(i,t)+wallTempIncreaseCharge),rowWallSurfaceAreaPCM,fcontWall,tankRadius,thicknessEncapsulation,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,timeStep,includeWallLoss);
                        storageStepPlate = 0;
                        
                        logQWallPCM(t) = logQWallPCM(t) + storageStepWall*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionPCM*densityFluid(i)*rowVolumePCM;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;

                        logPCMLoss(t) = logPCMLoss(t) + fluidStepHeatLoss*voidFractionPCM*densityPCMEff*rowVolumePCM; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionPCM*rowVolumePCM*densityFluid(i);
                    else
                        Re = G*rockAverageDiameter/viscosityFluid(i);
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);

                        fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);

                        fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                        storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                        storageStepConduction = getRockConduction(includeConduction,keffRocks(i-1),keffRocks(i),tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);

                        storageStepRadiation =0;% getRockRadiation(tempMatrixStorage(i,t));
                        storageStepWall = 0;
                        storageStepPlate = 0;

    %                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
    %                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
                        
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + logRockLoss(t);
                    end

                case rowStartRock
                    Re = G*rockAverageDiameter/viscosityFluid(i);
                    Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                    
                    fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                    fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);

                    fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                    %                     getFluidHeatLoss(lossInsulationTerms,densityFluid(i),tempMatrixFluid(i,t),ambientTemp,timeStep,rowWallSurfaceAreaRock,voidFractionRocks,rowVolumeRock,i);
                    storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                    storageStepConduction = getRockConduction(includeConduction,0,keffRocks(i),tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
%                     storageStepConduction = getRockConduction(includeConduction,keffRocks(i-1),keffRocks(i),tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
                    if includePCM
                        storageStepRadiation = getInterfaceRadiation(tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tankArea,voidFractionPCMInverse,emissivityEncapsulation,emissivityRocks,timeStep,voidFractionRocksInverse,densityRocks,rowVolumeRock,includeRadiation);
                    else
                        storageStepRadiation = 0;
                    end
                    storageStepWall = 0;
                    storageStepPlate = 0;
%                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
%                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
%                     energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
                    
                    logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                    logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
                otherwise % Rock section
                    Re = G*rockAverageDiameter/viscosityFluid(i);
                    Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                    
                    fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                    fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i-1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);

                    fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                    %                     getFluidHeatLoss(lossInsulationTerms,densityFluid(i),tempMatrixFluid(i,t),ambientTemp,timeStep,rowWallSurfaceAreaRock,voidFractionRocks,rowVolumeRock,i);
                    storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                    storageStepConduction = getRockConduction(includeConduction,keffRocks(i-1),keffRocks(i),tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);

                    storageStepRadiation = 0;%getRockRadiation(tempMatrixStorage(i,t));
                    storageStepWall = 0;
                    storageStepPlate = 0;
                   
%                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
%                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
%                     energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
                    
                    logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                    logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
            end
            intEnergyMatrixStorage(i,2) = intEnergyMatrixStorage(i,1) + storageStepConvection + storageStepConduction + storageStepRadiation + storageStepWall + storageStepPlate;
            intEnergyMatrixFluid(i,2) = intEnergyMatrixFluid(i,1) + fluidStepConvection + fluidStepAdvection - fluidStepHeatLoss;
            
            
            if t == logTime
                if includePCM && i < rowStartRock
                    fluidStepConvectionJ = fluidStepConvection*voidFractionPCM*densityFluid(i)*rowVolumePCM;
                    fluidStepAdvectionJ = fluidStepAdvection*voidFractionPCM*densityFluid(i)*rowVolumePCM;
                    fluidStepHeatLossJ =fluidStepHeatLoss*voidFractionPCM*densityFluid(i)*rowVolumePCM; 

                    storageStepConvectionJ = storageStepConvection*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                    storageStepConductionJ = storageStepConduction*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                    storageStepRadiationJ = storageStepRadiation*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;

                    rowEnergyStorageJ = intEnergyMatrixStorage(i,2)*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                    rowEnergyFluidJ = intEnergyMatrixFluid(i,2)*voidFractionPCM*densityFluid(i)*rowVolumePCM;
                else                
                    fluidStepConvectionJ = fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                    fluidStepAdvectionJ = fluidStepAdvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                    fluidStepHeatLossJ =fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 

                    storageStepConvectionJ = storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    storageStepConductionJ = storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    storageStepRadiationJ = storageStepRadiation*voidFractionRocksInverse*densityRocks*rowVolumeRock;

                    rowEnergyStorageJ = intEnergyMatrixStorage(i,2)*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    rowEnergyFluidJ = intEnergyMatrixFluid(i,2)*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                end
                
                logEnergyFlow(i,1:13) = [tempMatrixStorage(i,t),tempMatrixFluid(i,t),intEnergyMatrixStorage(i,2),intEnergyMatrixFluid(i,2),...
                    storageStepConvection, fluidStepConvection,storageStepConduction,fluidStepAdvection,...
                    storageStepRadiation,fluidStepHeatLoss,hvEffRocks,keffRocks(i),BiRocks];      
                logEnergyFlowJ(i,1:13) = [tempMatrixStorage(i,t),tempMatrixFluid(i,t),rowEnergyStorageJ,intEnergyMatrixFluid(i,2),...
                    storageStepConvectionJ, fluidStepConvectionJ,storageStepConductionJ,fluidStepAdvectionJ,...
                    storageStepRadiationJ,fluidStepHeatLossJ,hvEffRocks,kRocks(i),BiRocks];  
            end
        end
        logEnergyIn(t) = timeStep*mdot(t)*(intEnergyMatrixFluid(2,1)-intEnergyMatrixFluid(numRows-2,1)); %-timeStep*mdot(t)*(sum(intEnergyMatrixFluid(3:numRows-3,2))-sum(intEnergyMatrixFluid(3:numRows-3,1)));
%         energyInput = energyInput + timeStep*mdot(t)*(intEnergyMatrixFluid(2,1)-intEnergyMatrixFluid(numRows-2,1));
        energyInput = energyInput + logEnergyIn(t);
    else
        tempMatrixFluid(numRows-1:numRows,t) = dischargeTemp;
        intEnergyMatrixFluid(numRows-1:numRows,1) = getEnthalpy(tempMatrixFluid(numRows-1:numRows,t));
%     mdot(t) = vdot(t) * getDensity(dischargeTemp);

        for i=1:1:numRows
            BiRocks = hpRocks * rockAverageRadius/kRocks(i);
            hvEffRocks = hvRocks/(1+0.25*BiRocks);
            
            switch i
               case 1 % top plate
                    fluidStepConvection = 0;
                    fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
                    fluidStepHeatLoss =0; 
                    storageStepConvection = 0;
                    storageStepConduction = 0;
                    storageStepRadiation = 0;
                    storageStepWall = 0;
                    storageStepPlate = 0;
                case 2 % Tinlettop
                    fluidStepConvection = 0;
                    fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
                    fluidStepHeatLoss =0; 
                    storageStepConvection = 0;
                    storageStepConduction = 0;
                    storageStepRadiation = 0;
                    storageStepWall = 0;
                    storageStepPlate = 0;
                case numRows-1 %bottom plate
                    fluidStepConvection = 0;
                    fluidStepAdvection = 0; %getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i));
                    fluidStepHeatLoss = 0; %getFluidHeatLoss(0,densityFluid(i),tempMatrixFluid(i,t),ambientTemp,timeStep,rowWallSurfaceAreaRock,voidFractionRocks,rowVolumeRock);
                    storageStepConvection = 0;
                    storageStepConduction = 0;
                    storageStepRadiation = 0;
                    storageStepWall = 0;
                    storageStepPlate = 0;
                case numRows %Toutletbottom
                    fluidStepConvection = 0;
                    fluidStepAdvection = 0; %getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
                    fluidStepHeatLoss = 0;%getFluidHeatLoss(0,densityFluid(i),tempMatrixFluid(i,t),ambientTemp);
                    storageStepConvection = 0;
                    storageStepConduction = 0;
                    storageStepRadiation = 0;
                    storageStepWall = 0;
                    storageStepPlate = 0;
                case 3              % PCM or rocks
                    if includePCM
                        Re = 3.4658*G*dEncap/viscosityFluid(i);
                        % Based on umax not u0, thus increased by
                        % umax=ST/(SD-dEncap)*u0 = 3.4658*u0
                        
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                        PrS = getSpecificHeat(tempMatrixStorage(i))*getViscosity(tempMatrixStorage(i))/getConductivity(tempMatrixStorage(i));

                        NuPCM = 0.4845 * sqrt(Re)*Pr^0.37*(Pr/PrS)^0.25;
                                                
                        hpPCM = hvPCMModifier*NuPCM*kFluid(3)/dEncap;
                        hvPCM = hpPCM*surfaceAreaPCM;
                        BiPCM = hpPCM * rEncap / kPCMSolid;

                        hvEffPCM = hvPCM/(1+0.25*BiPCM);

                        fluidStepConvection = getFluidConvection(hvEffPCM,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionPCM);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumePCM,voidFractionPCM);
                        fluidStepHeatLoss = 0; %getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionPCM,dEncap,rowVolumePCM,rowWallSurfaceAreaPCM,timeStep);
                        storageStepConvection = getRockConvection(hvEffPCM,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionPCMInverse,densityPCMEff);
                        storageStepConduction = getPCMConduction(includeConduction,tempMatrixStorage(i,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,tankArea,fcontPCM,voidFractionPCMInverse,densityPCMEff,rowVolumePCM);
                        storageStepRadiation = getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixStorage(2,t),tankArea,voidFractionPlateInverse,emissivityEncapsulation,emissivityPlate,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeRadiation)+ ...
                            getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixStorage(4,t),tankArea,voidFractionPCMInverse,emissivityEncapsulation,emissivityEncapsulation,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeRadiation);
                        storageStepWall = getQWallPCM(tempMatrixStorage(3,t),tempMatrixStorage(3,t)+wallTempIncreaseDischarge,rowWallSurfaceAreaPCM,fcontWall,tankRadius,thicknessEncapsulation,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,timeStep,includeWallLoss);
                        storageStepPlate = getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixStorage(3,t)+plateTempIncreaseDischarge,plateArea,voidFractionPlateInverse,emissivityEncapsulation,emissivityPlate,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeCapRadiation) ;
                        
                        logQWallPCM(t) = logQWallPCM(t) + storageStepWall*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logQRadPlate(t) = logQRadPlate(t) + storageStepPlate*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionPCM*densityFluid(i)*rowVolumePCM;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;

                        logPCMLoss(t) = logPCMLoss(t) + fluidStepHeatLoss*voidFractionPCM*densityPCMEff*rowVolumePCM; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionPCM*rowVolumePCM*densityFluid(i);
                    else
                        Re = G*rockAverageDiameter/viscosityFluid(i);
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);

                        fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
                        fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                        storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                        storageStepConduction = getRockConduction(includeConduction,keffRocks(i+1),keffRocks(i),tempMatrixStorage(i+1,t),tempMatrixStorage(i,t),tempMatrixStorage(i-1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
                        storageStepRadiation =0;% getRockRadiation(tempMatrixStorage(i,t));
                        storageStepWall = 0;
                        storageStepPlate = getInterfaceRadiation(tempMatrixStorage(3,t),tempMatrixStorage(3,t)+plateTempIncreaseDischarge,plateArea,voidFractionPlateInverse,emissivityRocks,emissivityPlate,timeStep,voidFractionRocksInverse,densityRocks,rowVolumeRock,includeCapRadiation) ;

    %                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
    %                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
                        
                        logQRadPlate(t) = logQRadPlate(t) + storageStepPlate*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + logRockLoss(t);
                    end
                    
                case num2cell(4:rowCutoffPCM-1)              % PCM or rocks
                    if includePCM
                        Re = 3.4658*G*dEncap/viscosityFluid(i);
                        % Based on umax not u0, thus increased by
                        % umax=ST/(SD-dEncap)*u0 = 3.4658*u0
                        
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                        PrS = getSpecificHeat(tempMatrixStorage(i))*getViscosity(tempMatrixStorage(i))/getConductivity(tempMatrixStorage(i));

                        NuPCM = 0.4845 * sqrt(Re)*Pr^0.37*(Pr/PrS)^0.25;
                                                
                        hpPCM = hvPCMModifier*NuPCM*kFluid(3)/dEncap;
                        hvPCM = hpPCM*surfaceAreaPCM;
                        BiPCM = hpPCM * rEncap / kPCMSolid;

                        hvEffPCM = hvPCM/(1+0.25*BiPCM);

                        fluidStepConvection = getFluidConvection(hvEffPCM,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionPCM);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumePCM,voidFractionPCM);
                        fluidStepHeatLoss = 0; % getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionPCM,dEncap,rowVolumePCM,rowWallSurfaceAreaPCM,timeStep);
                        storageStepConvection = getRockConvection(hvEffPCM,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionPCMInverse,densityPCMEff);
                        storageStepConduction = getPCMConduction(includeConduction,tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,tankArea,fcontPCM,voidFractionPCMInverse,densityPCMEff,rowVolumePCM);
                        storageStepRadiation =getPCMRadiation(tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),tankArea,densityPCMEff,rowVolumePCM,emissivityEncapsulation,timeStep,includeRadiation);
                        storageStepWall = getQWallPCM(tempMatrixStorage(i,t),tempMatrixStorage(i,t)+wallTempIncreaseDischarge,rowWallSurfaceAreaPCM,fcontWall,tankRadius,thicknessEncapsulation,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,timeStep,includeWallLoss);
                        storageStepPlate = 0;
                        
                        logQWallPCM(t) = logQWallPCM(t) + storageStepWall*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionPCM*densityFluid(i)*rowVolumePCM;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;

                        logPCMLoss(t) = logPCMLoss(t) + fluidStepHeatLoss*voidFractionPCM*densityPCMEff*rowVolumePCM; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionPCM*rowVolumePCM*densityFluid(i);
                    else
                        Re = G*rockAverageDiameter/viscosityFluid(i);
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);

                        fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
                        fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                        storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                        storageStepConduction = getRockConduction(includeConduction,keffRocks(i+1),keffRocks(i),tempMatrixStorage(i+1,t),tempMatrixStorage(i,t),tempMatrixStorage(i-1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
                        storageStepRadiation =0;% getRockRadiation(tempMatrixStorage(i,t));
                        storageStepWall = 0;
                        storageStepPlate = 0;

    %                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
    %                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
                        
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + logRockLoss(t);
                    end
                    
                    case rowCutoffPCM              % PCM or rocks
                    if includePCM
                        Re = 3.4658*G*dEncap/viscosityFluid(i);
                        % Based on umax not u0, thus increased by
                        % umax=ST/(SD-dEncap)*u0 = 3.4658*u0
                        
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                        PrS = getSpecificHeat(tempMatrixStorage(i))*getViscosity(tempMatrixStorage(i))/getConductivity(tempMatrixStorage(i));

                        NuPCM = 0.4845 * sqrt(Re)*Pr^0.37*(Pr/PrS)^0.25;
                                                
                        hpPCM = hvPCMModifier*NuPCM*kFluid(3)/dEncap;
                        hvPCM = hpPCM*surfaceAreaPCM;
                        BiPCM = hpPCM * rEncap / kPCMSolid;

                        hvEffPCM = hvPCM/(1+0.25*BiPCM);

                        fluidStepConvection = getFluidConvection(hvEffPCM,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionPCM);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumePCM,voidFractionPCM);
                        fluidStepHeatLoss = 0; % getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionPCM,dEncap,rowVolumePCM,rowWallSurfaceAreaPCM,timeStep);
                        storageStepConvection = getRockConvection(hvEffPCM,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionPCMInverse,densityPCMEff);
                        storageStepConduction = getPCMConduction(includeConduction,tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i,t),timeStep,tankArea,fcontPCM,voidFractionPCMInverse,densityPCMEff,rowVolumePCM);
                        storageStepRadiation = getInterfaceRadiation(tempMatrixStorage(i,t),tempMatrixStorage(i-1,t),tankArea,voidFractionPCMInverse,emissivityEncapsulation,emissivityEncapsulation,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeRadiation)+ ...
                            getInterfaceRadiation(tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),tankArea,voidFractionPCMInverse,emissivityEncapsulation,emissivityRocks,timeStep,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,includeRadiation);
                        storageStepWall = getQWallPCM(tempMatrixStorage(i,t),tempMatrixStorage(i,t)+wallTempIncreaseDischarge,rowWallSurfaceAreaPCM,fcontWall,tankRadius,thicknessEncapsulation,voidFractionPCMInverse,densityPCMEff,rowVolumePCM,timeStep,includeWallLoss);
                        storageStepPlate = 0;
                        
                        logQWallPCM(t) = logQWallPCM(t) + storageStepWall*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionPCM*densityFluid(i)*rowVolumePCM;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionPCMInverse*densityPCMEff*rowVolumePCM;

                        logPCMLoss(t) = logPCMLoss(t) + fluidStepHeatLoss*voidFractionPCM*densityPCMEff*rowVolumePCM; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionPCM*rowVolumePCM*densityFluid(i);
                    else
                        Re = G*rockAverageDiameter/viscosityFluid(i);
                        Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);

                        fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                        fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
                        fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                        storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                        storageStepConduction = getRockConduction(includeConduction,keffRocks(i+1),keffRocks(i),tempMatrixStorage(i+1,t),tempMatrixStorage(i,t),tempMatrixStorage(i-1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
                        storageStepRadiation =0;% getRockRadiation(tempMatrixStorage(i,t));
                        storageStepWall = 0;
                        storageStepPlate = 0;

    %                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
    %                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
                        
                        logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                        logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                        logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
%                         energyWallLossRocksCharging = energyWallLossRocksCharging + logRockLoss(t);
                    end
                    
                case rowStartRock % Rock section
                    Re = G*rockAverageDiameter/viscosityFluid(i);
                    Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                    
                    fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                    fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);

                    fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                    storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                    storageStepConduction = getRockConduction(includeConduction,keffRocks(i+1),0,tempMatrixStorage(i+1,t),tempMatrixStorage(i,t),tempMatrixStorage(i-1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
                    if includePCM
                        storageStepRadiation = getInterfaceRadiation(tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tankArea,voidFractionPCMInverse,emissivityEncapsulation,emissivityRocks,timeStep,voidFractionRocksInverse,densityRocks,rowVolumeRock,includeRadiation);
                    else
                        storageStepRadiation = 0;
                    end
                    storageStepWall = 0;
                    storageStepPlate = 0;
%                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
%                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
%                     energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
                    
                    logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                    logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 

                otherwise % Rock section
                    Re = G*rockAverageDiameter/viscosityFluid(i);
                    Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
                    
                    fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
                    fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);

                    fluidStepHeatLoss = getFluidHeatLoss(includeHeatLoss,i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
                    %                     getFluidHeatLoss(lossInsulationTerms,densityFluid(i),tempMatrixFluid(i,t),ambientTemp,timeStep,rowWallSurfaceAreaRock,voidFractionRocks,rowVolumeRock,i);
                    storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
                    storageStepConduction = getRockConduction(includeConduction,keffRocks(i+1),keffRocks(i),tempMatrixStorage(i+1,t),tempMatrixStorage(i,t),tempMatrixStorage(i-1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
                    storageStepRadiation = 0;%getRockRadiation(tempMatrixStorage(i,t));
                    storageStepWall = 0;
                    storageStepPlate = 0;

                   
%                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
%                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluidF,densityFluid(i),G);
%                     energyWallLossRocksCharging = energyWallLossRocksCharging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
                    
                    logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
                    logConduction(t) = logConduction(t) + storageStepConduction*voidFractionRocksInverse*densityRocks*rowVolumeRock;
                    logRockLoss(t) = logRockLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
%                 otherwise %case 3:numRows-2  % Rock section
%                     fluidStepConvection = getFluidConvection(hvEffRocks,densityFluid(i),tempMatrixStorage(i,t),tempMatrixFluid(i,t),timeStep,voidFractionRocks);
%                     fluidStepAdvection = getFluidAdvection(mdot(t),densityFluid(i),intEnergyMatrixFluid(i+1),intEnergyMatrixFluid(i),timeStep,rowVolumeRock,voidFractionRocks);
%                     fluidStepHeatLoss = getFluidHeatLoss(i,heffWall(i),tempMatrixFluid(i,t),ambientTemp,lossInsulationTerms,densityFluid(i),kFluid(i),Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep);
%                     storageStepConvection = getRockConvection(hvEffRocks,tempMatrixFluid(i,t),tempMatrixStorage(i,t),timeStep,voidFractionRocksInverse,densityRocks);
%                     storageStepConduction = getRockConduction(keffRocks(i),tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
%                     storageStepRadiation = getRockRadiation(tempMatrixStorage(i,t));
%                    
% %                     viscosityFluidF = getViscosity(tempMatrixFluid(i,t));
% %                     pressureDrop(t) = pressureDrop(t) + getPressureDrop(viscosityFluid(i),densityFluid(i),G);
%                     energyWallLossRocksDischarging = energyWallLossRocksDischarging + fluidStepHeatLoss*voidFractionRocks*rowVolumeRock*densityFluid(i);
% 
%                     logConvectionStorage(t) = logConvectionStorage(t) + storageStepConvection*voidFractionRocksInverse*densityRocks*rowVolumeRock;
%                     logConvectionFluid(t) = logConvectionFluid(t) + fluidStepConvection*voidFractionRocks*densityFluid(i)*rowVolumeRock;
%                     logLoss(t) = logLoss(t) + fluidStepHeatLoss*voidFractionRocks*densityFluid(i)*rowVolumeRock; 
            end
            
            intEnergyMatrixStorage(i,2) = intEnergyMatrixStorage(i,1) + storageStepConvection + storageStepConduction + storageStepRadiation + storageStepWall + storageStepPlate;
            intEnergyMatrixFluid(i,2) = intEnergyMatrixFluid(i,1) + fluidStepConvection + fluidStepAdvection - fluidStepHeatLoss;
            
%             totalEnergyInput = totalEnergyInput + timeStep*mdot(t)*(intEnergyMatrixFluid(2,1)-intEnergyMatrixFluid(numRows-2,1));
            
            if t == logTime
                logEnergyFlow(i,1:13) = [tempMatrixStorage(i,t),tempMatrixFluid(i,t),intEnergyMatrixStorage(i,2),intEnergyMatrixFluid(i,2),storageStepConvection, ...
                    fluidStepConvection,storageStepConduction,fluidStepAdvection,storageStepRadiation,fluidStepHeatLoss,...
                    hvEffRocks,keffRocks(i),BiRocks]; 
            end
            logEnergyIn(t) = timeStep*mdot(t)*(-intEnergyMatrixFluid(2,1)+intEnergyMatrixFluid(numRows-2,1))-timeStep*mdot(t)*(sum(intEnergyMatrixFluid(3:numRows-3,2))-sum(intEnergyMatrixFluid(3:numRows-3,1)));
        end
        
    end
    
    if useMode == 1
        intEnergyMatrixFluid(1:2,2) = getEnthalpy(tempMatrixFluid(1,t+1));
    else
        intEnergyMatrixFluid(numRows,2) = getEnthalpy(dischargeTemp);
%         tempMatrixFluid(numRows-1,t+1) = tempMatrixFluid(,t+1);
    end
    
    intEnergyMatrixStorage(1:2,2) = intEnergyMatrixStorage(3,2);
    intEnergyMatrixStorage(numRows-1:numRows,2) = intEnergyMatrixStorage(numRows-2,2);
    
    %  
    if includePCM
    logEnergyStored(t) = (sum(intEnergyMatrixStorage(3:rowCutoffPCM,2)) - sum(intEnergyMatrixStorage(3:rowCutoffPCM,1)))*densityPCMEff*rowVolumePCM*voidFractionPCMInverse ...
        + (sum(intEnergyMatrixStorage(rowStartRock:numRows-2,2)) - sum(intEnergyMatrixStorage(rowStartRock:numRows-2,1)))*densityRocks*rowVolumeRock*voidFractionRocksInverse;
    else
        logEnergyStored(t) = (sum(intEnergyMatrixStorage(3:numRows-2,2)) - sum(intEnergyMatrixStorage(3:numRows-2,1)))*densityRocks*rowVolumeRock*voidFractionRocksInverse;
    end
    
    try
    [tempMatrixStorage(:,t+1),tempMatrixFluid(:,t+1)] = calculateTemps(intEnergyMatrixStorage(:,2),intEnergyMatrixFluid(:,2),tempMatrixStorage(:,t+1),tempMatrixFluid(:,t+1),pcmLowEnergyCutoff, pcmHighEnergyCutoff, TmeltStart, TmeltEnd, ceffPCMSolid, ceffPCMPhaseChange, ceffPCMLiquid,  rowCutoffPCM, rowStartRock, getEnthalpy,getEnergyDerivative);
%     disp(tempMatrixStorage(1:rowCutoffPCM,t+1))
    %     tempMatrixStorage(2,t+1) = tempMatrixStorage(3,t+1)+abs(tempMatrixFluid(3,t+1) - tempMatrixStorage(3,t+1))*0.5;
    catch MEx
        disp(t)
        disp(MEx.identifier)
        disp(MEx.message)
        return
    end
    if includePCM
        for i = 3:rowCutoffPCM
            pressureDrop(t) = pressureDrop(t) + getPressureDrop(rowLengthPCM,dEncap, voidFractionPCM, voidFractionPCMInverse, tempMatrixFluid(i,t+1)-tempMatrixFluid(i,t), tempMatrixFluid(i,t),viscosityFluid(i),densityFluid(i),G);
        end
        for i = rowStartRock:numRows-2
            pressureDrop(t) = pressureDrop(t) + getPressureDrop(rowLengthRock,rockAverageDiameter, voidFractionRocks, voidFractionRocksInverse, tempMatrixFluid(i,t+1)-tempMatrixFluid(i,t), tempMatrixFluid(i,t),viscosityFluid(i),densityFluid(i),G);
        end
    else
        for i = 3:numRows-2
            pressureDrop(t) = pressureDrop(t) + getPressureDrop(rowLengthRock,rockAverageDiameter, voidFractionRocks, voidFractionRocksInverse, tempMatrixFluid(i,t+1)-tempMatrixFluid(i,t), tempMatrixFluid(i,t),viscosityFluid(i),densityFluid(i),G);
        end
    end
    
%     intEnergyMatrixStorage(:,1) = intEnergyMatrixStorage(:,2);
%     intEnergyMatrixFluid(:,1) = intEnergyMatrixFluid(:,2);
        

    
    if mod(t,progressBarPercentTime) == 0
        clc
        timeSoFar = toc;
        fprintf('Progress: %g %%\n',round(t/totalTime*100,0))
        fprintf('Time elapsed so far: %g min\n',timeSoFar/60)
        fprintf('Estimated time to completion: %g min\n',timeSoFar*(totalTime/t-1)/60)
        if t < 2*onePercentTime && timeSoFar<5
            progressBarPercentTime = 10*onePercentTime;
        end
    end
    if t == logTime
        fprintf('Updated energyflowlog at t=%g\n',t)
        if breakAfterLogTime
            return
        end
        
    end
    t=t+1;
end

disp('Transient analysis complete')
toc
%% Energy logger calculations
% {
energyInputKWH = energyInput/3.6e6;
if totalTime > chargeTime+1
    energyStoredRocks = (totalEnergyRocksEndCharge - totalEnergyRocksStart)*densityRocks*rowVolumeRock*voidFractionRocksInverse;
else
    energyStoredRocks = (sum(intEnergyMatrixStorage(3:numRows-2)) - totalEnergyRocksStart)*densityRocks*rowVolumeRock*voidFractionRocksInverse;
end
    % energyStoredRocks2 = (sum(storageEnergyMatrixChargeEndLog(3:70)) - totalEnergyRocksStart)*densityRocks*rowVolumeRock*voidFractionRocksInverse;
energyStoredRocksKWH = energyStoredRocks/3.6e6;
if includePCM
    energyStoredPCM = (totalEnergyPCMEndCharge - totalEnergyPCMStart)*densityPCMEff*rowVolumePCM*voidFractionPCMInverse;
    energyStoredPCMKWH = energyStoredPCM/3.6e6;
    energyWallLossPCMChargingKWH = sum(logPCMLoss(1:chargeTime))/3.6e6;
    energyWallLossPCMDischargingKWH = sum(logPCMLoss(chargeTime+1:totalTime))/3.6e6;
    energyWallPCMKWHCharge = sum(logQWallPCM(1:chargeTime))/3.6e6;
    energyWallPCMKWHDischarge = sum(logQWallPCM(chargeTime+1:totalTime))/3.6e6;
end
energyWallLossRocksChargingKWH = sum(logRockLoss(1:chargeTime))/3.6e6;
energyWallLossRocksDischargingKWH = sum(logRockLoss(chargeTime+1:totalTime))/3.6e6;

energyRadPlateKWHCharge = sum(logQRadPlate(1:chargeTime))/3.6e6;
energyRadPlateKWHDischarge = sum(logQRadPlate(chargeTime+1:totalTime))/3.6e6;
%% Written Output

fprintf('Total energy input: %g kWh -- target 18.4 for PCM or 19.3 rocks\n',energyInputKWH)
fprintf('Total energy stored in rocks: %g kWh\n',energyStoredRocksKWH)
if includePCM
    fprintf('Total energy stored in PCM: %g kWh\n',energyStoredPCMKWH)
end
fprintf('Energy exchange with top plate during charge: %g kWh -- target 0.95 or 1.01 rocks\n',energyRadPlateKWHCharge);
fprintf('Energy exchange with top plate during discharge: %g kWh -- target -0.14 or +0.42 rocks\n',energyRadPlateKWHDischarge);
fprintf('Energy lost from rocks section during charging: %g kWh -- target 0.24 or 0.3 rocks\n',energyWallLossRocksChargingKWH)
fprintf('Energy lost from rocks section during discharging: %g kWh -- target 0.25 or 0.3 rocks\n',energyWallLossRocksDischargingKWH)
if includePCM
    fprintf('Energy exchange between PCM section and wall during charging: %g kWh -- target 1.27\n',energyWallPCMKWHCharge)
    fprintf('Energy exchange between PCM section and wall during discharging: %g kWh -- target 0.43\n',energyWallPCMKWHDischarge)
end
%}

%% Figure setup
close all
% Row selection based on PCM use
switch conditionSelection
    case 'Zanganeh2015'
        if includePCM
            R1 = 2+19;     R2 = 2+36;    R3 = 2+52;    R4 = 2+67; R0 = 2+3;
            legendEntry0 = '0 mm';
            legendEntry1 = '379 mm';
            legendEntry2 = '720 mm';
            legendEntry3 = '1045 mm';
            legendEntry4 = '1345 mm';
            legendEntry5 = ['tank top ('  fluidName  ')'];
            legendEntry6 = [fluidName ' mass flow (-->)'];
        %     I = imread('..\Model Outputs\Source Figures\Fig12aAverageTempsPCMpng');
        else
            R1 = 11+2;    R2 = 29+2;    R3 = 43+2;    R4 = 67+2;  R0 = 3;
            legendEntry0 = '0 mm';
            legendEntry1 = '225 mm';
            legendEntry2 = '574 mm';
            legendEntry3 = '865 mm';
            legendEntry4 = '1335 mm';
            legendEntry5 = ['tank top ('  fluidName  ')'];
            legendEntry6 = [fluidName ' mass flow (-->)'];
        %     I = imread('..\Model Outputs\Source Figures\Fig11aAverageTempsRocksOnly.png');
        end
    case 'Beasley1989-low'
        R1 = 12;
        R2 = 25;
        R3 = 44;
        R4 = 44;
        legendEntry1 = 'x/L = 0.25';
        legendEntry2 = 'x/L = 0.575';
        legendEntry3 = '--';
        legendEntry4 = '--';
    case 'Beasley1989-high'
        R1 = 12;
        R2 = 22;
        R3 = 44;
        R4 = 44;
        legendEntry1 = 'x/L = 0.25';
        legendEntry2 = 'x/L = 0.5';
        legendEntry3 = '--';
        legendEntry4 = '--';
end

plotTimeStart = 1;
plotTimeEnd = totalTime;

fileNameTime = clock;

mem = memory;
if mem.MemUsedMATLAB > 5e9
    return
end

%% Pressure Drop
%{
figure
hold on

plot(pressureDrop,'LineWidth',3)

hold off
ylabel('Pressure Drop (Pa)')
    
xticks(1*timeStepInverse*3600*[1 2 3 4 5 6 7 8])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
%ylabel('Temperature (C)')

legend('Pressure drop across tank (Pa)')

xlim(3600*timeStepInverse*[0 8])
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

% saveas(gcf,figureName)
figureNamePressureDropPNG = ['..\Model Outputs\PressureDrop\pressureDrop' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
% saveas(fig,figureNamePressureDropPNG)
% fprintf('Wrote to file: %s\n',figureNamePressureDropPNG)
pause(3)
close gcf
%}

%% Mass and volume flow rate
%{
figure
hold on

ylim([0 0.02]);
% uistack(h,'bottom');

% plot(mdot2,'LineWidth',6)
plot(mdot,'LineWidth',2)
ylabel('Mass flow (kg/s)')


yyaxis right
plot(vdot*2118.88,'LineWidth',2)
ylabel('Volume flow (cfm)')
ylim([0 30]);

hold off
    
xticks(1*timeStepInverse*3600*[1 2 3 4 5 6 7 8])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
%ylabel('Temperature (C)')

legend('mdot','vdot')

xlim(3600*timeStepInverse*[0 8])
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

% saveas(gcf,figureName)
figureNameFlowRatePNG = ['..\Model Outputs\FlowRate\flowRate' datestr(fileNameTime,'mmm-dd-yyyy HH-MM-SS') '.png'];
saveas(fig,figureNameFlowRatePNG)
fprintf('Wrote to file: %s\n',figureNameFlowRatePNG)
pause(3)
close gcf
%}
%% Average temperature in tank
if plotAverage
figure
hold on

% plot((tempMatrixStorage(R0,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R0,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)
plot((tempMatrixStorage(R1,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R1,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)
plot((tempMatrixStorage(R2,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R2,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)
plot((tempMatrixStorage(R3,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R3,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)
plot((tempMatrixStorage(R4,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R4,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)

% plot((tempMatrixFluid(2,plotTimeStart:plotTimeEnd)+tempMatrixStorage(2,plotTimeStart:plotTimeEnd))/2,'LineWidth',3)
plot(tempMatrixFluid(2,plotTimeStart:plotTimeEnd),'LineWidth',3)
plot(TmeltStart*ones(1,plotTimeEnd),'LineWidth',2)

ylim([20 750]);
xlabel('Time (hours)')
ylabel(['Temperature (' char(176) 'C)'])
% uistack(h,'bottom');

yyaxis right
plot(mdot,'LineWidth',2)
if includePCM
    ylim([0 0.02])
else
ylim([0 0.02]);
end


hold off
ylabel('Mass Flow (kg/s)')
    
xticks(1*timeStepInverse*3600*[1 2 3 4 5 6 7 8])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
%ylabel('Temperature (C)')

legendEntry0 = '20 mm';
% legend(legendEntry1,legendEntry2,legendEntry3,legendEntry4,legendEntry5,legendEntry6,legendEntry0)

if includePCM
    title('Average temperature in tank - PCM + rocks')
else
%     title('Average temperature in tank - Rocks only')
    title([num2str(energyRadPlateKWHCharge) '---' num2str(energyRadPlateKWHDischarge) '---' num2str(energyRadPlateKWHDischarge) '---' num2str(energyWallLossRocksDischargingKWH)])
end
   
xlim(3600*timeStepInverse*[0 8])
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

figureNameAveragePNG = ['..\Model Outputs\Averages\averageTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
saveas(fig,figureNameAveragePNG)
fprintf('Wrote to file: %s\n',figureNameAveragePNG)
end
%}

%% Storage temperature in tank
%{
figure
hold on

% plot((tempMatrixStorage(R0,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R0,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)
plot(tempMatrixStorage(R1,plotTimeStart:plotTimeEnd),'LineWidth',2)
plot(tempMatrixStorage(R2,plotTimeStart:plotTimeEnd),'LineWidth',2)
plot(tempMatrixStorage(R3,plotTimeStart:plotTimeEnd),'LineWidth',2)
plot(tempMatrixStorage(R4,plotTimeStart:plotTimeEnd),'LineWidth',2)


ylim([20 750]);
xlabel('Time (hours)')
ylabel(['Temperature (' char(176) 'C)'])
% uistack(h,'bottom');



yyaxis right
plot(mdot,'LineWidth',2)
ylim([0 0.02]);


hold off
ylabel('Mass Flow (kg/s)')
    
xticks(1*timeStepInverse*3600*[1 2 3 4 5 6 7 8])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
%ylabel('Temperature (C)')

legendEntry0 = '20 mm';
% legend(legendEntry1,legendEntry2,legendEntry3,legendEntry4,legendEntry5,legendEntry6,legendEntry0)

if includePCM
    title('Storage temperature in tank - PCM + rocks')
else
    title('Storage temperature in tank - Rocks only')
end
xlim(3600*timeStepInverse*[0 8])
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 600, 630, 395]);
fig.PaperPositionMode = 'auto';

figureNameStoragePNG = ['..\Model Outputs\Storage\storageTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
% saveas(fig,figureNameStoragePNG)
fprintf('Wrote to file: %s\n',figureNameStoragePNG)
%}

%% PCM temperature in tank
% {
if plotPCM
if includePCM
    figure
    hold on

    % plot((tempMatrixStorage(R0,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R0,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)
    plot(tempMatrixStorage(3,plotTimeStart:plotTimeEnd),'LineWidth',2)
    plot(tempMatrixFluid(3,plotTimeStart:plotTimeEnd),'LineWidth',2)
    plot(tempMatrixStorage(5,plotTimeStart:plotTimeEnd),'LineWidth',2)
    plot(tempMatrixFluid(5,plotTimeStart:plotTimeEnd),'LineWidth',2)


    ylim([550 650]);
    xlabel('Time (hours)')
    ylabel(['Temperature (' char(176) 'C)'])
    % uistack(h,'bottom');





    hold off

    xticks(1*timeStepInverse*3600*[1 2 3 4 5 6 7 8])
    xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
    %ylabel('Temperature (C)')

    legendEntry0 = '20 mm';
    legend('PCM 1','Fluid 1','PCM 3','Fluid 3')

    title('PCM/Fluid temperature in tank')
    xlim(3600*timeStepInverse*[1.5 5])
    fig = gcf;
    set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 600, 630, 395]);
    fig.PaperPositionMode = 'auto';

    figureNamePCMPNG = ['..\Model Outputs\PCM\PCMTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
    saveas(fig,figureNamePCMPNG)
    fprintf('Wrote to file: %s\n',figureNamePCMPNG)
end
end
%}

%%
% {
t=1;
while t*timeStepInverse < totalTime
    logTempFluid(:,t) = tempMatrixFluid(:,t*timeStepInverse);
    logTempStorage(:,t) = tempMatrixStorage(:,t*timeStepInverse);
    logPressure(t) = pressureDrop(t*timeStepInverse);
    logmdot(t) = mdot(t*timeStepInverse);
    logvdot(t) = vdot(t*timeStepInverse);
    logAverageTemp(:,t) = (logTempFluid(:,t)+logTempStorage(:,t))/2;
    
    t=t+1;
end
if includePCM
    parameters = {'Fluid: ',fluidName; 'Reference fluid: ',referenceFluid; 'PCM Included?', includePCM; 'PCM Height', heightPCM;'Ein',energyInputKWH;'Erocks',energyStoredRocksKWH;'EPCM',energyStoredPCMKWH;'Eradplatech',energyRadPlateKWHCharge;'Eradplatedisch',energyRadPlateKWHDischarge;'EwallPCMch',energyWallPCMKWHCharge;'EwallPCMdisch',energyWallPCMKWHDischarge;'Ewallrocksch',energyWallLossRocksChargingKWH;'Ewallrocksdisch',energyWallLossRocksDischargingKWH};
    else
    parameters = {'Fluid: ',fluidName; 'Reference fluid: ',referenceFluid; 'PCM Included?', includePCM; 'PCM Height', 'n/a'; 'Ein',energyInputKWH; 'Erocks',energyStoredRocksKWH;'EPCM', 'n/a';'Eradplatech',energyRadPlateKWHCharge;'Eradplatedisch',energyRadPlateKWHDischarge;'EwallPCMch', 'n/a';'EwallPCMdisch','n/a'; 'Ewallrocksch',energyWallLossRocksChargingKWH;'Ewallrocksdisch',energyWallLossRocksDischargingKWH};
end
workSpaceName = ['..\Model Outputs\Temp Histories\temps' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.xlsx'];
% writematrix(A,workSpaceName)
% type(workSpaceName)
warning('off','MATLAB:xlswrite:AddSheet');
try
    disp('Attempting to write temperature histories to file...')
    xlswrite(workSpaceName,parameters,'Sheet1','A1:B26');
    disp('Wrote model parameters to file')
    xlswrite(workSpaceName,logTempFluid','Fluid',['A1:BT' num2str(length(logTempFluid))]);
    disp('Wrote fluid temperatures to file')
    xlswrite(workSpaceName,logTempStorage','Storage',['A1:BT' num2str(length(logTempStorage))]);
    disp('Wrote storage temperatures to file')
    xlswrite(workSpaceName,logPressure','Pressure',['A1:BT' num2str(length(logPressure))]);
    disp('Wrote pressure drop to file')
    xlswrite(workSpaceName,logmdot','mdot',['A1:A' num2str(length(logmdot))]);
    disp('Wrote mass flow rate to file')
    xlswrite(workSpaceName,logvdot','vdot',['A1:A' num2str(length(logvdot))]);
    disp(['---> Successfully wrote temperature histories to file: ' workSpaceName])
catch
    disp('Error writing to file')
end
%}


% transVis = true; %transVis;
% transVis = input('Figure loaded. Press 1 to run transient visualization, or any other number to skip: ');
if transVis
%% Plot temp gradient at each time sequentially
    figure
    tic
    startRow = 3;

    endRow = numRows-2;
    for t = 1:timeStepInverse*60*1:chargeTime
        % % f
%         t= logTime

        clf
%         figure
        hold on
        yyaxis right
        ylim([0 100])
        plot(abs(tempMatrixFluid(startRow:endRow,t)-tempMatrixStorage(startRow:endRow,t)),'LineWidth',3)

        yyaxis left
        ylim([ambientTemp 705])
        plot(tempMatrixFluid(startRow:endRow,t),'LineWidth',3)
        plot(tempMatrixStorage(startRow:endRow,t),'LineWidth',3)
%         plot([2 2], [0 655])
%         plot([3 3], [0 655])
        plot([rowCutoffPCM rowCutoffPCM], [0 705])
%         plot([rowCutoffPCM-2 rowCutoffPCM-2], [0 705])
        xlim([1 numRows])


        title(['Time: ' num2str(round(t*timeStep/60,1)) ' min'])
        pause(0.05)
    %     if mod(t,timeStepInverse*600)==1
    %         pause()
    %     end
        hold off
    end
    disp('Done')
    toc
else
    close gcf
end