%% top
clc
clear all
close all
addpath(genpath('fluid properties'))

useTopLogTime = false;
% useTopLogTime = true;
logTime = 4632000;
breakAfterLogTime = false;
% breakAfterLogTime = true;

tempControlledTime = false;
	cycleNumber = 1;
	chargeFluidTempThreshold = 500;
	dischargeFluidTempThreshold = 600;


conditionSelection = 'test1'; % options are Zanganeh2015, Trahan, TrahanPCM, Beasley1989-low, Beasley1989-high, ...
shsChoice = 'Zanganeh2015';    % 'NIST' 'Smetana2013' 'Zanganeh2015' 'Trahan'
pcmChoice = 'KCl-NaCl-Tm=657';           % 'AlSi12', 'KCl-NaCl', 'KCl-NaCl-Tm=680', 'KCl-NaCl-Tm=670'
encapsChoice = 'steel';         % steel, ptfe, ceramic
hvFunction = 'Zanganeh2015';    % 'Zanganeh2015' 'Trahan' 'Beasley'
pcmFraction = 0;
cMinConstantCharge = false;
cMinConstantAbsMin = false;

% ********** Modifiers applied here for sensitivity analysis
mdotModifierChargeAir = 1;			% 0.9, 1, 1.1
mdotModifierDischargeAir = 1;
mdotModifierAir = mdotModifierChargeAir;
mdotModifierChargeCO2 = 1.1;
mdotModifierDischargeCO2 = 1;
mdotModifierCO2 = mdotModifierChargeCO2;
specificHeatModifier = 0;			% +/- 0.2
enthalpyFusionModifier = 0;			% +/- 0.2
effectivenessHX = 0.95; % 0.95 assumed based on HT recuperators in table 3 of https://doi.org/10.1016/j.energy.2016.06.013
% ********** End modifiers used for sensitivity analysis

useHX = true;
efficiencySCO2 = @(temp) 0.50; % NYI - based on Quadrennial Technology Review 2015, chapter 4 (Advancing Clean Electric Power Technologies), pg 9

fluidTempAlgorithm = 'nr'; % options are 'cp' or 'nr'
SHSTempAlgorithm = 'quadratic'; % options are 'quadratic' or 'cp' or 'nr'
minimumTemp = 400;	maximumTemp = 700;	avgTemp = (400+700)/2;

continueFromLastRun = true;
saveForNextRun = true;

constantFluidProperties = false;
includeConduction = true;
includeBulkHeatLoss = true;
TrahanWallLoss = false;

includeRadiation = false;           % treating PCM as packed bed = radiation included in conduction term
includePCMWallExchange = false;     % Not used because no reference temps available 
includeCapRadiation = false;        % Not used because no reference temps available 

useInternalEnergy = false;

tankFluidName = 'DryAirAtmo';           % options are CO2HP, CO2LP, DryAir, DryAirAtmo, MoistAir
inputMode = 1;                      % options are 1 (mass flow) or 2 (volume flow). NYI
bypassFraction = 0.0;					% 0 because large particle to tank ratio used

increaseFluidTemps = false;         % Always true (via override later) for CO2, 
increaseBaseTemp = false;           % can be true or false for air

runTransient = true;

if pcmFraction == 0
	includePCM = false;
	includeSHS = true;
elseif pcmFraction == 1
	includePCM = true;
	includeSHS = false;
else
	includePCM = true;
	includeSHS = true;
end

%% Plot selector etc -- moved to individual conditions
transVis= false;				transVisMdot = false;
saveFigures = false;			dumpToExcel = true;
closeFiguresAfterPlotting = false;
compareToValidation = false;	printEnergyFlows = true;
printLNorms = false;

plotFigures = false;
plotFluid = true;      
plotStorage = true;				plotmdot = false;
calculatePressureDrop = false;

%% Condition set
switch conditionSelection
	case 'test1'
		voidFractionPCM = 0.4; voidFractionPCMInverse = 1-voidFractionPCM;
        voidFractionSHS = 0.4; voidFractionSHSInverse = 1-voidFractionSHS;
        
        
        chargeTimeHours =  8;					% 8
        dischargeTimeHours = 8; %8;					% 5 standard
		chargeTimePartialHour = 0;				% Unneeded, preserved to avoid changing other code
        totalTimeHours = chargeTimeHours+dischargeTimeHours;
        rowCountModifier = 1;
        hvModifier = 1; 
		plateTempIncreaseCharge = 0;			plateTempIncreaseDischarge = 0;
		wallTempIncreaseCharge = zeros(1,52);	wallTempIncreaseDischarge = zeros(1,52);
		
        useTrahanTempProfile = false;
        trackError = false;
        trackFluidError = false;
    otherwise
        disp('Invalid condition selection')
        return
end
plateTempIncrease = plateTempIncreaseCharge;
wallTempIncrease = wallTempIncreaseCharge;

if dischargeTimeHours > 0
    includeDischarge = true;
else
    includeDischarge = false;
end

switch fluidTempAlgorithm
    case 'cp'
        calculateTempFluid = @calculateTempFluidCP;
		calculateTempsFluid = @calculateTempsFluidCP;
    case 'nr'
        calculateTempFluid = @calculateTempFluidNR;
		calculateTempsFluid = @calculateTempsFluidNR;
end
switch SHSTempAlgorithm
    case 'quadratic'
        calculateTempSHS = @calculateTempSHSQuadratic;
%         (energy)(-747.0995+sqrt(747.0995^2+4*0.2838*energy))/2/0.2838;
    case 'cp'
        calculateTempSHS = @calculateTempSHSCP;
    case 'nr'
        calculateTempSHS = @calculateTempSHSNR;
end

switch conditionSelection
	case 'test1'
		tempProfile = @tempProfileTest1;
		getmdot = @(x,y,z) 2000; % kg/s
		getStorageTopConduction = @getPackedBedConduction;
end

%% ============== Time setup
switch conditionSelection
	case 'test1'
        timeStepInverse = rowCountModifier*100;  %getmdot(1,1,1)/2*              % 41        % samples per second
	otherwise
        timeStepInverse = rowCountModifier*100;  %getmdot(1,1,1)/2*              % 41        % samples per second
end

timeStep = 1/timeStepInverse;                       % seconds per sample
chargeTime = 0; dischargeTime = 0;
chargeTime = floor((chargeTimeHours+chargeTimePartialHour)*3600*timeStepInverse);        % 3+1/6 hours in seconds
dischargeTime = floor(dischargeTimeHours*3600*timeStepInverse);        % Normally 5h, shortened for quicker runtime while testing charge mode

totalTime = floor((chargeTime+dischargeTime)) + 2;
onePercentTime = round(0.01*totalTime);
progressBarPercentTime = onePercentTime;
t = 1;

if ~useTopLogTime
    logTime = chargeTime-1; 
end

%% Encapsulation properties

[densityEncaps,cEncaps,getkEncaps,emissivityEncaps] = getPropertiesEncaps(encapsChoice);

%% PCM properties
diameterPCMParticle = 0.03; % m
encapsulationThickness = 1.0e-03;  % m

[TmeltStart, TmeltEnd, dTMelt, cPCMSolid, cPCMLiquid, enthalpyFusion, kPCMSolid, kPCMLiquid,densityPCMSolid,densityPCMLiquid] = getPropertiesPCM(pcmChoice);


cPCMSolid = (1+specificHeatModifier)*cPCMSolid;
cPCMLiquid = (1+specificHeatModifier)*cPCMLiquid;
enthalpyFusion = (1+enthalpyFusionModifier)*enthalpyFusion;

Tmelt = TmeltStart + dTMelt/2;
cPCMPhaseChange = (cPCMSolid+cPCMLiquid)/2 + enthalpyFusion/dTMelt; % J/kg

diameterPCM = diameterPCMParticle - 2*encapsulationThickness;
radiusPCMParticle = diameterPCMParticle/2;
radiusPCM = diameterPCM/2;
volumePCMInParticle = 4/3*pi*radiusPCM^3;                     % PCM only
volumePCMParticle = 4/3*pi*(radiusPCMParticle)^3;   % includes encapsulation
surfaceAreaPCMParticle = 4*pi*radiusPCMParticle^2;

massPCMinParticle = densityPCMLiquid * volumePCMInParticle;
massEncapsinParticle = densityEncaps * (volumePCMParticle-volumePCMInParticle);
massParticle = massPCMinParticle+massEncapsinParticle;

cEffPCMSolid = (massPCMinParticle*cPCMSolid + massEncapsinParticle*cEncaps)/(massParticle); % J/kgK
cEffPCMLiquid = (massPCMinParticle*cPCMLiquid + massEncapsinParticle*cEncaps)/(massParticle); % J/kgK
cEffPCMPhaseChange = (massPCMinParticle*cPCMPhaseChange + massEncapsinParticle*cEncaps)/(massParticle);
densityPCMEff = massParticle / volumePCMParticle;

pcmSensibleEnergyPerKG = ((TmeltStart-minimumTemp)*cEffPCMSolid + (maximumTemp - TmeltEnd)*cEffPCMLiquid); % J/kg
pcmLatentEnergyPerKG = dTMelt*cEffPCMPhaseChange; % J/kg
pcmEnergyPerKG = (pcmSensibleEnergyPerKG+pcmLatentEnergyPerKG);

%% SHS properties
diameterSHS = 0.03;

[getEnergySHS,getSpecificHeatSHS,getkSHS,densitySHS,emissivitySHS] = getPropertiesSHS(shsChoice);

SHSsensibleEnergyPerKG = (getEnergySHS(maximumTemp) - getEnergySHS(minimumTemp)); % J/kg
radiusSHS = diameterSHS/2;
volumeSHSParticle = 4/3*pi*radiusSHS^3;
surfaceAreaSHSParticle = 4*pi*radiusSHS^2;

%% ============== Tank dimensions
designPowerElectric = 100e6;        % W
designPowerCycleEfficiency = 0.45;  % cycle efficiency
designStorageHours = 12;            % hours
% designEnergyThermal = designPowerElectric * designStorageHours / designPowerCycleEfficiency; % Wh-th
% 
% totalTankEnergyPerM3KWH = (pcmEnergyPerM3*tankFractionPCM + SHSEnergyPerM3*(1-tankFractionPCM)); % J

% designTankVolume = designEnergyThermalJ/totalTankEnergyPerM3KWH; % m3

tankAspectRatio = 0.5; % H/D
tankDiameter = 40; %(4*designTankVolume/pi / tankAspectRatio)^(1/3);
tankRadius = tankDiameter/2;
tankHeight = 20; %tankDiameter*tankAspectRatio;
tankArea = pi*tankRadius^2;
tankVolume = pi*tankRadius^2*tankHeight;
plateArea = tankArea;
voidFractionPlate = 1;	voidFractionPlateInverse = 0;	emissivityPlate = 1;


volumeSHSTank = (1-pcmFraction) * tankVolume;
volumePCMTank = pcmFraction * tankVolume;

pcmMassKG = volumePCMTank*densityPCMEff*voidFractionPCMInverse;
shsMassKG = volumeSHSTank*densitySHS*voidFractionSHSInverse;

energyTankPCM = pcmEnergyPerKG*pcmMassKG; % J
energyTankSHS = shsMassKG*SHSsensibleEnergyPerKG; % J
energyTankTotal = energyTankPCM + energyTankSHS;

designEnergyThermalKWH = energyTankTotal/3.6e6; % GWh-th
designEnergyThermalGWH = designEnergyThermalKWH/1e6; % GWh-th

rowsAtTopAndBottom = 1;

numRowsStorage = rowCountModifier*50;
numRows = numRowsStorage + 2*rowsAtTopAndBottom;
numRowsPCM = round(numRowsStorage*pcmFraction,0);
numRowsSHS = numRowsStorage - numRowsPCM;

if includePCM
	if pcmFraction == 1
		rowCutoffPCM = numRows-1;
		rowStartSHS = numRows; 
	else
		rowCutoffPCM = rowsAtTopAndBottom+numRowsPCM;
		rowStartSHS = rowCutoffPCM+1; 
	end
else
	rowCutoffPCM = 0;
	rowStartSHS = 1; 
end

chargeMode = true;

%% ============== Discretization
rowLength(1:numRows) = tankHeight/numRows;                     % All model rows same length
rowVolume = rowLength*pi*(tankRadius)^2;
rowWallSurfaceArea = pi*2*tankRadius*rowLength;

%% ============== Insulation parameters
switch conditionSelection
    case 'test1'
		thicknessTank = 0.0127/2; % 0.25 inch
		thicknessInsulation = 0.1524; % 6 inches
		kTank = 16.2; % W /mK, using AISI 304 http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MQ304A
		kInsulation = 0.038; % W/mK, using rockwool from Zang 2015
		
        lossInsulationTerms(1:numRows) = tankRadius*...
			(1/kTank*log((tankRadius + thicknessTank)/tankRadius) ...
			+1/kInsulation*log((tankRadius + thicknessTank+thicknessInsulation)/(tankRadius + thicknessTank)));
    case 'Zanganeh2015'
        % Adjust these values for actual conditions. Whatever that means?
        lossInsulationTermsbackup = lossInsulationTerms;
        lossInsulationTerms(1:26) = lossInsulationTermsbackup(1);
        lossInsulationTerms(27) = (lossInsulationTermsbackup(1)+lossInsulationTermsbackup(2))/2;
        lossInsulationTerms(28:47) = lossInsulationTermsbackup(2);
        lossInsulationTerms(48:numRows) = lossInsulationTermsbackup(3);
        clear lossInsulationTermsbackup
end

%% ============== Row properties
particlesPerRowPCM = rowVolume*voidFractionPCMInverse / volumePCMParticle;
rowSurfaceAreaPCM = particlesPerRowPCM*surfaceAreaPCMParticle;

particlesPerRowSHS = rowVolume*voidFractionSHSInverse / volumeSHSParticle;
rowSurfaceAreaSHS = particlesPerRowSHS*surfaceAreaSHSParticle;   % m2/row

%% Fluid selection
[getEnthalpy,getSpecificHeat,getConductivity,getViscosity,getDensity,fluidPressure,getEntropy] = getPropertiesFluid(tankFluidName);

if constantFluidProperties
    fluidEnthalpy = getEnthalpy(avgTemp);
    getEnthalpy = @(x) getSpecificHeat(1)*x;
    fluidConductivity = getConductivity(avgTemp);
    getConductivity = @(x) fluidConductivity;
    fluidViscosity = getViscosity(avgTemp);
    getViscosity = @(x) fluidViscosity;
    fluidDensity = getDensity(avgTemp);
    getDensity = @(x) fluidDensity;
    fluidSpecificHeat = getSpecificHeat(avgTemp);
    getSpecificHeat = @(x) fluidSpecificHeat;
    getEnergyDerivative = getSpecificHeat;
end

%% ============= Temperature conditions
switch conditionSelection
    case 'test1'
        chargeTemp = 700;
        dischargeTemp = 400;      % C - temperature of HTF at bottom inlet during discharge
        baseTemp = 400;           % C - starting temperature of rock bed
end
ambientTemp = 25;

%% HX Setup
if useHX
cpCO2HXCharge = 1276.735826; % J/kg, from refprop: 25 MPa, 700 C
cpAirHXCharge = 1068.759482; % J/kg, from refprop: 0.1013 MPa, 400 C
cpCO2HXDischarge = 1248.461149; % J/kg, from refprop: 25 MPa, 700 C
cpAirHXDischarge = 1136.094885; % J/kg, from refprop: 0.1013 MPa, 700 C

CminCharge = getmdot() * min(mdotModifierChargeAir,mdotModifierChargeCO2) * min(cpAirHXCharge,cpCO2HXCharge);
if effectivenessHX == 1
	CminDischarge = getmdot() * min(mdotModifierDischargeAir,mdotModifierDischargeCO2) * min(cpAirHXDischarge,cpCO2HXDischarge);
else
	CminDischarge = getmdot() * min(mdotModifierDischargeAir,mdotModifierDischargeCO2) * min(cpAirHXDischarge,cpCO2HXDischarge);
end
if cMinConstantCharge
	CminDischarge = CminCharge;
end
if cMinConstantAbsMin
	CminCharge = min(CminCharge,CminDischarge);
	CminDischarge = CminCharge;
end

dischargeTempEnthalpyAir = getEnthalpy(dischargeTemp);
dischargeTempEnthalpyCO2 = getEnthalpyCO2HP(dischargeTemp);
end

% Avoid computer crashes
% if totalTime > 5e6
%     disp('Array size too large. Shutting down to prevent computer crash.')
%     return
% end

%% ============== Temp/energy/step/void fraction Matrix setup

tempStorage = zeros(totalTime,numRows);			% C
tempFluid = zeros(totalTime,numRows);			% C
if tempControlledTime && continueFromLastRun
	if cycleNumber > 1
		if pcmFraction > 0
			workSpaceName = ['.\Outputs\Excel Dumps\Tm=' num2str(Tmelt) ' -- X=' num2str(pcmFraction) ' -- cycle' num2str(cycleNumber-1) '.xlsx'];
		else
			workSpaceName = ['.\Outputs\Excel Dumps\rockOnly -- cycle ' num2str(cycleNumber-1) '.xlsx'];
		end
		prevRunTemps = xlsread(workSpaceName,'TankTempsAtEnd');
		tempFluid(1,:) = prevRunTemps(:,2)';		% C
		tempStorage(1,:) = prevRunTemps(:,1)' ;		% C
	else
		tempStorage(1,:) = baseTemp;				% C
	end
else
    tempStorage(1,:) = baseTemp;				% C
end
tempFluid(1,:) = tempStorage(1,:);				% C / assume initial condition of equilibrium within tank
% Above assumption is reasonable if tank sits idle overnight, because fluid
% and storage should equilibrate in temperature.

intEnergyStorage = zeros(2,numRows);
if includePCM
	for i = 1:rowCutoffPCM
		if tempStorage(1,i) < TmeltStart
			intEnergyStorage(1,i) = cEffPCMSolid * tempStorage(1,i);
		elseif tempStorage(1,i) < TmeltEnd
			intEnergyStorage(1,i) = cEffPCMSolid*TmeltStart + cEffPCMPhaseChange*(tempStorage(1,i)-TmeltStart);
		else
			intEnergyStorage(1,i) = cEffPCMSolid*TmeltStart + cEffPCMPhaseChange*(TmeltEnd-TmeltStart)+ cEffPCMLiquid*(tempStorage(1,i)-TmeltEnd);
		end
	end
end
intEnergyStorage(1,rowStartSHS:numRows) = getEnergySHS(tempStorage(1,rowStartSHS:numRows));
if includePCM
energySHSStart = (intEnergyStorage(1,rowStartSHS:numRows-rowsAtTopAndBottom));   % J/kg, kg stays constant
else
energySHSStart = (intEnergyStorage(1,2:numRows-rowsAtTopAndBottom));   % J/kg, kg stays constant
end
if includePCM
	energyPCMStart = (intEnergyStorage(1,2:rowCutoffPCM));          % J/kg, kg stays constant
	voidFraction(1:rowCutoffPCM) = voidFractionPCM;
	densityStorage(1:rowCutoffPCM) = densityPCMEff;
	storageParticleDiameter(1:rowCutoffPCM) = diameterPCMParticle;
	emissivityStorage(1:rowCutoffPCM) = emissivityEncaps;
end

voidFraction(rowStartSHS:numRows) = voidFractionSHS;
voidFractionInverse(:) = 1-voidFraction(:);
densityStorage(rowStartSHS:numRows) = densitySHS;
storageParticleDiameter(rowStartSHS:numRows) = diameterSHS;
storageParticleRadius = storageParticleDiameter/2;
emissivityStorage(rowStartSHS:numRows) = emissivitySHS;

rowMassStorage(1:numRows) = rowVolume.*densityStorage.*voidFractionInverse;
totalStorageMass = sum(rowMassStorage(rowsAtTopAndBottom+1:numRows-rowsAtTopAndBottom));


intEnergyFluid = zeros(2,numRows);
enthalpyFluid = zeros(2,numRows);
tempHXCO2 = zeros(totalTime,2);
enthalpyHXCO2 = zeros(totalTime,2);
entropyHX = zeros(totalTime,2);
exergy = zeros(totalTime,3);
logExergyInputRecirculating = zeros(totalTime,1);
logExergyInputTopTempOnly = zeros(totalTime,1);

enthalpyFluid(1,:) = getEnthalpy(tempFluid(1,:));
intEnergyFluid(1,:) = getEnthalpy(tempFluid(1,:));

intEnergyStorage(2,:) = intEnergyStorage(1,:);
intEnergyFluid(2,:) = intEnergyFluid(1,:);
enthalpyFluid(2,:) = enthalpyFluid(1,:);


%% heat transfer multipliers setup
storageRadiationPlateCoefficient = timeStep*5.67e-8*plateArea*voidFractionPlateInverse/(1/emissivityStorage(2)+1/emissivityPlate-1);
storageRadiationPCMCoefficient = timeStep*5.67e-8*tankArea*voidFractionPCMInverse/(2/emissivityEncaps-1);
storageRadiationInterfaceCoefficient = timeStep*5.67e-8*tankArea*voidFractionPCMInverse/(1/emissivityEncaps+1/emissivitySHS-1);

% Enthalpy coefficient matrix
enthalpyCoefficientCharge = zeros(numRows);
enthalpyCoefficientDischarge = zeros(numRows);
for i = 2:1:numRows-1
    enthalpyCoefficientCharge(i,i) = -1;
    enthalpyCoefficientCharge(i,i-1) = 1;
    enthalpyCoefficientDischarge(i,i) = -1;
    enthalpyCoefficientDischarge(i,i+1) = 1;
end
enthalpyCoefficientCharge(numRows,numRows) = -1;
enthalpyCoefficientCharge(numRows,numRows-1) = 1;
enthalpyCoefficientDischarge(1,1) = -1;
enthalpyCoefficientDischarge(1,2) = 1;

enthalpyCoefficient = enthalpyCoefficientCharge;

% Conduction coefficient matrix
conductionCoefficient = zeros(numRows);
for i = 2:1:numRows-1
    conductionCoefficient(i,i) = -2;
    conductionCoefficient(i,i-1) = 1;
    conductionCoefficient(i,i+1) = 1;
end

%% ============== Fluid Base properties and flow rates
densityFluid = getDensity(tempFluid(1,:));
viscosityFluid = getViscosity(tempFluid(1,:));
kFluid = getConductivity(tempFluid(1,:));
specificHeatFluid = getSpecificHeat(tempFluid(1,:));


mdotAir = zeros(1,totalTime);
mdotCO2 = zeros(1,totalTime);
vdot = zeros(1,totalTime);
bypassFractionInverse = 1-bypassFraction;               % bypass fraction 

%% Stability calculations
chargeMode = true;

switch conditionSelection
    case 'test1'
        mdotmax = max(mdotModifierDischargeAir,mdotModifierChargeAir)*getmdot(); %0.014; % Needs to be updated
        tempMax = 700;
        VFmin = 0.4;
        VFmax = 0.4;
        minrowlength = min(rowLength);
        rhomin = getDensity(tempMax);
        rhominsolid = min(densityPCMEff,densitySHS);
        kminsolid = 0.5;
        specificheatminsolid = 788.6150;
end

vmax = mdotmax / VFmin / tankArea / rhomin;

CFL = vmax * timeStep / minrowlength;
maxdTfluid = minrowlength/vmax;
mintimestepinversefluid = 1/maxdTfluid;

maxdTsolid = 0.5*minrowlength^2*(1-VFmax)*rhominsolid * specificheatminsolid / kminsolid;
mintimestepinversesolid = 1/maxdTsolid;
fprintf('dt = %g -- maxdTfluid = %g -- maxdTsolid = %g\n',timeStep,maxdTfluid,maxdTsolid)
fprintf('timestepinverse = %g -- minTSinversefluid = %g -- minTSinversesolid = %g\n',timeStepInverse,mintimestepinversefluid,mintimestepinversesolid)
if CFL >= 1
    disp('fluid stability check failed. Quitting.')
    return
end
if timeStep >= maxdTsolid
    disp('solid stability check failed. Quitting.')
    return
end

%% ============== Energy and heat transfer variables preallocation

breakFlag = false;

fluidAdvection = zeros(1,numRows);
fluidHeatLoss = zeros(1,numRows);
convectionExchange = zeros(1,numRows);
storageConduction = zeros(1,numRows);
storageRadiation = zeros(1,numRows);
storagePlate  = zeros(1,numRows);
QWallPCM = zeros(1,numRows);

pressureDrop(totalTime) = 0;

q(totalTime) = 0;

energyInput = 0;
energyExtracted = 0;
energyRadiatedPlateCharging = 0;
energyRadiatedPlateDischarging = 0;
energyWallLossPCMCharging = 0;
energyWallLossPCMDischarging = 0;
energyWallLossRocksCharging = 0;
energyWallLossRocksDischarging = 0;
energyStoredRocks = 0;
energyStoredPCM = 0;

%% Loggers
logConvectionStorage(totalTime) = 0;
logConduction(totalTime) = 0;
logConvectionFluid(totalTime) = 0;
logRadiationPlate(totalTime) = 0;
logFluidHeatLoss(totalTime) = 0;
logPCMLoss(totalTime) = 0; 
logEnergyInput(totalTime) = 0;
logEnergyExtracted(totalTime) = 0;
logEnergyStored(totalTime) = 0;
logQRadPlate(totalTime) = 0;
logQWallPCM(totalTime) = 0;

logFluidEnergyChange(totalTime) = 0;
logFluidTempChangeTheoretical(totalTime) = 0;
logFluidTempChangeActual(totalTime) = 0;

maxConduction = 0;
minConduction = 0;

%% Figure setup and excel dump
close all
% Row selection based on PCM use
R0 = rowsAtTopAndBottom;    
R1 = rowsAtTopAndBottom+round(numRowsStorage*0.1,0);     
R2 = rowsAtTopAndBottom+round(numRowsStorage*0.2,0);     
R3 = rowsAtTopAndBottom+round(numRowsStorage*0.3,0); 
R4 = rowsAtTopAndBottom+round(numRowsStorage*0.4,0); 
R5 = rowsAtTopAndBottom+round(numRowsStorage*0.5,0); 
R6 = rowsAtTopAndBottom+round(numRowsStorage*0.6,0);     
R7 = rowsAtTopAndBottom+round(numRowsStorage*0.7,0); 
R8 = rowsAtTopAndBottom+round(numRowsStorage*0.8,0); 
R9 = rowsAtTopAndBottom+round(numRowsStorage*0.9,0); 
Rbottom = numRows;

xlimHours = [0 totalTimeHours];
ylimFlow = [0 2];
ylimTemps = [400 700];

color0 = [255 50 50]/255;
color1 = [255 154 52]/255;
color2 = [255 255 5]/255;
color3 = [80 148 52]/255;
color4 = [52 105 255]/255;
color5 = [0 50 154]/255;
color6 = [102 102 204]/255;
color7 = [102 0 153]/255;
color8 = [255 103 206]/255;
color9 = [0 0 0]/255;
            
plotTimeStart = 1;
plotTimeEnd = totalTime;

%% Memory check for computer crash prevention
mem = memory;
if mem.MemUsedMATLAB > 15e9
    fprintf('\nMemory used by matlab is %g gigabytes, stopping to minimize risk of computer crash.',mem.MemUsedMATLAB/1e9)
    return
end

%% Transient
if ~runTransient
    disp('Setup complete. Transient analysis not performed.')
    return
end

tic
disp('Starting transient analysis...')
startTime = clock;
fprintf('Time: %g:%g:%g\n',startTime(4),startTime(5),floor(startTime(6)))

totalTime = totalTime - 1;
while t < totalTime && breakFlag == false % && cycleNumber <= numCycles
    % Swap energy parameters to set up for new time step
    intEnergyStorage(1,:) = intEnergyStorage(2,:);
    intEnergyFluid(1,:) = intEnergyFluid(2,:);
        
    % Fluid properties for time step
    densityFluid(:) = getDensity(tempFluid(t,:));
    viscosityFluid(:) = getViscosity(tempFluid(t,:));
    kFluid(:) = getConductivity(tempFluid(t,:));
    specificHeatFluid(:) = getSpecificHeat(tempFluid(t,:));
	Pr = specificHeatFluid.*viscosityFluid./kFluid;
            
	enthalpyFluid = getEnthalpy(tempFluid(t,:));
    rowMassFluid(:) = rowVolume(:).*densityFluid(:).*voidFraction(:);
    
	
    % Storage conductivity
    if includePCM
        % [rowCutoffPCM rowStartSHS]
        kStorage (rowStartSHS:numRows) = getkSHS(tempStorage(t,rowStartSHS:numRows));
        for i=1:1:rowCutoffPCM
            if tempStorage(t,i) < TmeltStart
                kStorage(i) = 17.7;
            elseif tempStorage(t,i) < TmeltEnd
                kStorage(i) = 22;
            else
                kStorage(i) = 23.1;
            end
        end
    else
        for i=1:1:numRows
            kStorage(i) = getkSHS(tempStorage(t,i));
        end
    end
      
   % Switch to discharge mode, log energy before proceeding
	if tempControlledTime
		if chargeMode && t>1 && tempHXCO2(t-1,2) > chargeFluidTempThreshold
			chargeMode = false;
			chargeTime = t-1;
			if includePCM % log energy at end of charge before switching to discharge
				energySHSEndCharge = (intEnergyStorage(1,rowStartSHS:numRows-1));
				storageEnergyChargeEndLog = intEnergyStorage(1,:);
				if strcmp(conditionSelection,'TrahanPCM')
					energyPCMEndCharge = (intEnergyStorage(1,2:numRows-1));
				else
					energyPCMEndCharge = (intEnergyStorage(1,2:rowCutoffPCM));
				end
			else
				energySHSEndCharge = (intEnergyStorage(1,2:numRows-1));
				storageEnergyChargeEndLog = intEnergyStorage(1,:);
			end
			% logInternalEnergy = intEnergyMatrix(:,2);
			enthalpyCoefficient = enthalpyCoefficientDischarge;

			plateTempIncrease = plateTempIncreaseDischarge;
			wallTempIncrease = wallTempIncreaseDischarge;
			mdotModifierAir = mdotModifierDischargeAir;
			mdotModifierCO2 = mdotModifierDischargeCO2;
			
			disp('Switching to discharge mode')
				disp(t)
		elseif ~chargeMode && tempHXCO2(t-1,1) < dischargeFluidTempThreshold
			enthalpyCoefficient = enthalpyCoefficientCharge;

			chargeMode = true;
			dischargeTime = t-1;
			plateTempIncrease = plateTempIncreaseCharge;
			wallTempIncrease = wallTempIncreaseCharge;
			mdotModifierAir = mdotModifierChargeAir;
			mdotModifierCO2 = mdotModifierChargeCO2;
			
% 			if cycleNumber == numCycles
				totalTime = t;
% 			else
% 				disp('Switching to charge mode')
% 				disp(t)
% 			end
% 			cycleNumber = cycleNumber+1;
			breakFlag = true;
		end
		
	else
		if mod(t,(chargeTime+dischargeTime)) == chargeTime+1 || (t == totalTime && chargeMode)
			chargeMode = false;
			if includePCM % log energy at end of charge before switching to discharge
				energySHSEndCharge = (intEnergyStorage(1,rowStartSHS:numRows-1));
				storageEnergyChargeEndLog = intEnergyStorage(1,:);
				if strcmp(conditionSelection,'TrahanPCM')
					energyPCMEndCharge = (intEnergyStorage(1,2:numRows-1));
				else
					energyPCMEndCharge = (intEnergyStorage(1,2:rowCutoffPCM));
				end
			else
				energySHSEndCharge = (intEnergyStorage(1,2:numRows-1));
				storageEnergyChargeEndLog = intEnergyStorage(1,:);
			end
			% logInternalEnergy = intEnergyMatrix(:,2);
			enthalpyCoefficient = enthalpyCoefficientDischarge;

			plateTempIncrease = plateTempIncreaseDischarge;
			wallTempIncrease = wallTempIncreaseDischarge;
			mdotModifierAir = mdotModifierDischargeAir;
			mdotModifierCO2 = mdotModifierDischargeCO2;
			disp('Switching to discharge mode')
    % Switch back to charge mode for extra cycles
		elseif t > chargeTime && mod(t,(chargeTime+dischargeTime)) == 1 
			enthalpyCoefficient = enthalpyCoefficientCharge;

			chargeMode = true;
			plateTempIncrease = plateTempIncreaseCharge;
			wallTempIncrease = wallTempIncreaseCharge;
			mdotModifierAir = mdotModifierChargeAir;
			mdotModifierCO2 = mdotModifierChargeCO2;
			disp('Switching to charge mode')
			disp(t)
		end
	end
%     fprintf('===timestep = %g===\n',t)
	if chargeMode
		mdotAir(t) = mdotModifierAir*bypassFractionInverse*getmdot(t,chargeTime,timeStep);
		mdotCO2(t) = mdotModifierCO2*bypassFractionInverse*getmdot(t,chargeTime,timeStep);
	else
		mdotAir(t) = mdotModifierAir*bypassFractionInverse*getmdot(t,chargeTime,timeStep);
		mdotCO2(t) = mdotModifierCO2*bypassFractionInverse*getmdot(t,chargeTime,timeStep);
	end
    G=mdotAir(t)/tankArea;
	[hp,hv] = calculateHvTrahan(G,storageParticleDiameter,viscosityFluid,Pr,0,kFluid,0,0,voidFractionInverse);
%     hvRocks = hvRocksCoeff*G^0.92;  % W/m2K
    
    % Inlet fluid temp depends on flow direction
    if chargeMode
        if ~useHX
            tempFluid(t,2) = tempProfile(true); % different from Dec1 but in a good way
        else
            tempHXCO2(t,1) = tempProfile(true);
            enthalpyHXCO2(t,1) = getEnthalpyCO2HP(tempHXCO2(t,1)); 
            q(t) = effectivenessHX*CminCharge*(tempProfile(true)-tempFluid(t,numRows));
            enthalpyHXCO2(t,2) = enthalpyHXCO2(t,1) - q(t)/(mdotCO2(t));
            
            tempFluid(t+1,1) = calculateTempFluidNR(enthalpyFluid(numRows)+q(t)/(mdotAir(t)),true, chargeTemp, 0,getEnthalpy,getSpecificHeat); 
            tempHXCO2(t,2) = calculateTempFluidNR(enthalpyHXCO2(t,2),true, chargeTemp, 0,@getEnthalpyCO2HP,@getCpCO2HP); 
		end
        intEnergyFluid(1,1) = getEnthalpy(tempFluid(t,1)); % different from Dec1. 
		% But this makes me wonder: should tempProfile be updating inlet at
		% t+1 or inlet at t? Seems like t would make more sense, no?
    else
        if ~useHX
            tempFluid(t+1,numRows) = tempProfile(false);
		else
            q(t) = effectivenessHX*CminDischarge*(tempFluid(t,1)-tempProfile(false));
			enthalpyFluid(1,numRows) = enthalpyFluid(1,1) - q(t)/(mdotAir(t));
			tempFluid(t+1,numRows) = calculateTempFluidNR(enthalpyFluid(1,numRows),false, chargeTemp, 0,getEnthalpy,getSpecificHeat); 
		end
	end
		
    % Fluid advection
    fluidAdvection = timeStep*mdotAir(t).*(enthalpyCoefficient*enthalpyFluid(1,:)')';

    % Fluid convection
        Re = G*storageParticleDiameter./viscosityFluid;
        % Convection - rock section
        if includeSHS
            BiotNumber(rowStartSHS:numRows) = hp(rowStartSHS:numRows) .* storageParticleRadius(rowStartSHS:numRows)./kStorage(rowStartSHS:numRows);
            hvEff(rowStartSHS:numRows) = hv(rowStartSHS:numRows)./(1+0.25*BiotNumber(rowStartSHS:numRows));
		else
			if includePCM
				BiotNumber(numRows) = hp(numRows) .* storageParticleRadius(numRows)./ kPCMSolid; %kPCMSolid(numRows);
				hvEff(numRows) = hv(numRows)./(1+0.25*BiotNumber(numRows)); % 0.25 in zang, 0.2 in Trahan
			end
		end
        % Convection - PCM section
        if includePCM
            BiotNumber(1:rowCutoffPCM) = hp(1:rowCutoffPCM) .* storageParticleRadius(1:rowCutoffPCM)./ kPCMSolid; %kStorage(1:rowCutoffPCM);
            hvEff(1:rowCutoffPCM) = hv(1:rowCutoffPCM)./(1+0.25*BiotNumber(1:rowCutoffPCM)); % 0.25 in zang, 0.2 in Trahan
		end
        hvEff(1) = 0; hvEff(numRows) = 0;
		convectionExchange = timeStep.*hvEff.*rowVolume.*(tempFluid(t,:)-tempStorage(t,:));

	% Dispersion and loss coefficients
	[keffStorage,heffWall] = calculatekeff(kStorage,kFluid,tempFluid(t,:),tempStorage(t,:),voidFraction,storageParticleDiameter,numRows,emissivityStorage,Re,Pr);
    
    % Wall loss - spherical particle packed bed section
    if includeBulkHeatLoss
        if TrahanWallLoss
            alphaConv = kFluid./storageParticleDiameter.*((0.203*(Re.*Pr).^(0.333333)) + (0.220*Re.^0.8.*Pr.^0.4));
        else
            alphaConv = kFluid./storageParticleDiameter.*(3.22*(Re.*Pr).^(0.333333)+0.117*Re.^0.8.*Pr.^0.4);
        end
        uInside = alphaConv + heffWall;
        uWall = 1./(1./uInside + lossInsulationTerms);
        fluidHeatLoss = timeStep.*uWall.*rowWallSurfaceArea.*(ambientTemp-tempFluid(t,:));
		fluidHeatLoss(1) = 0; fluidheatLoss(numRows) = 0;
    end
    % Wall loss - PCM section
    if includePCMWallExchange && includePCM
        QWallPCM(1:rowCutoffPCM) = timeStep.*kStorage(1:rowCutoffPCM)*storageWallHXareaTerms.*(wallTempIncrease(1:rowCutoffPCM));
	end
    	
    % Storage radiation
    if includeRadiation && includePCM && includeSHS
        % between rows 2,3 (cap calculated separately)
        storageRadiation(2) = storageRadiationPCMCoefficient*((tempStorage(t,3)+273.15)^4-(tempStorage(t,2)+273.15)^4);
        % rows i-1, i, i+1 from 3 to rowCutoffPCM-1
        storageRadiation(3:rowCutoffPCM-1) = storageRadiationPCMCoefficient* ...
            ((tempStorage(t,2:rowCutoffPCM-2)+273.15).^4 + (tempStorage(t,4:rowCutoffPCM)+273.15).^4 - 2*(tempStorage(t,3:rowCutoffPCM-1)+273.15).^4);
        % row rowCutoffPCM-1 to rowCutoffPCM and rowCutoffPCM to rowStartSHS
        storageRadiation(rowCutoffPCM) = storageRadiationPCMCoefficient * ((tempStorage(t,rowCutoffPCM-1)+273.15).^4-(tempStorage(t,rowCutoffPCM)+273.15).^4) + ...
            storageRadiationInterfaceCoefficient * ((tempStorage(t,rowCutoffPCM+1)+273.15).^4-(tempStorage(t,rowCutoffPCM)+273.15).^4);
        % row rowStartSHSs to rowCutoffPCM
        storageRadiation(rowStartSHS) =  storageRadiationInterfaceCoefficient * ((tempStorage(t,rowCutoffPCM)+273.15).^4-(tempStorage(t,rowCutoffPCM+1)+273.15).^4);
    end
    
    % Cap radiation - only relevant in Zanganeh case
    if includeCapRadiation
        storagePlate(2) = storageRadiationPlateCoefficient*((tempStorage(t,2)+plateTempIncrease+273.15)^4-(tempStorage(t,2)+273.15)^4);
    end
    
    if includeConduction
    % Storage conduction    % may need to add caveats to conduction between sections.
% 		if chargeMode
			if includePCM
				for i = 2:1:rowCutoffPCM-1  % first row of PCM through second to last row of PCM
					conductionCoefficient(i,i) = -(keffStorage(i)+keffStorage(i-1));
					conductionCoefficient(i,i-1) = keffStorage(i-1);
					conductionCoefficient(i,i+1) = keffStorage(i);
				end
				for i = rowStartSHS+1:1:numRows-1  % Second row of rock to last row of rock
					conductionCoefficient(i,i) = -(keffStorage(i)+keffStorage(i-1));
					conductionCoefficient(i,i-1) = keffStorage(i-1);
					conductionCoefficient(i,i+1) = keffStorage(i);
				end
				% last row of PCM
					conductionCoefficient(rowCutoffPCM,rowCutoffPCM) = -(keffStorage(rowCutoffPCM-1));
					conductionCoefficient(rowCutoffPCM,rowCutoffPCM-1) = keffStorage(rowCutoffPCM-1);
					conductionCoefficient(rowCutoffPCM,rowCutoffPCM+1) = 0;
				% first row of rock
					conductionCoefficient(rowStartSHS,rowStartSHS) = -(keffStorage(rowStartSHS));
					conductionCoefficient(rowStartSHS,rowStartSHS-1) = 0;
					conductionCoefficient(rowStartSHS,rowStartSHS+1) = keffStorage(rowStartSHS);
% 				% pcm section
% 				storageConduction(1:rowCutoffPCM) = timeStep*tankArea*fcontPCM./rowLength(1:rowCutoffPCM).* ...
% 					(conductionCoefficient(1:rowCutoffPCM,1:rowCutoffPCM)*tempStorage(t,1:rowCutoffPCM)')';
				% rock section
				storageConduction(1:numRows) = timeStep*tankArea./rowLength(1:numRows).* ...
					(conductionCoefficient(1:numRows,1:numRows)*tempStorage(t,(1:numRows))')';
			else
				for i = 2:1:numRows-1
					conductionCoefficient(i,i) = -(keffStorage(i)+keffStorage(i-1));
					conductionCoefficient(i,i-1) = keffStorage(i-1);
					conductionCoefficient(i,i+1) = keffStorage(i);
				end
				storageConduction = timeStep*tankArea./rowLength.*(conductionCoefficient*tempStorage(t,:)')';
			end
	end
	if max(storageConduction) > maxConduction
		maxConduction=max(storageConduction);
	end
	if min(storageConduction) > minConduction
		minConduction=max(storageConduction);
	end
    intEnergyStorage(2,:) = intEnergyStorage(1,:) + (convectionExchange + storageConduction  + storageRadiation + storagePlate + QWallPCM)./rowMassStorage;
    intEnergyFluid(2,:) = intEnergyFluid(1,:) + (fluidAdvection - convectionExchange + fluidHeatLoss)./rowMassFluid;
    
    if t == logTime
%         QWallPCM = zeros(1,numRows);
%         storagePlate = zeros(1,numRows);
        logEnergyFlow = [tempStorage(t,:);	tempFluid(t,:);			intEnergyStorage(1,:);
			intEnergyStorage(2,:);			intEnergyFluid(1,:);	intEnergyFluid(2,:); 
            convectionExchange;				storageConduction;		fluidAdvection; 
			fluidHeatLoss;					hvEff;					keffStorage;
			BiotNumber]';      
    end
    logConvectionStorage(t) = sum(-convectionExchange.*voidFractionInverse.*densityStorage.*rowVolume);
    logQRadPlate(t) = storagePlate(2);
    logQWallPCM(t) = sum(QWallPCM);
    logFluidHeatLoss(t) = sum(fluidHeatLoss);
	
% 	if useHX
% 		if chargeMode
% 			logEnergyInput(t) = timeStep*mdotAir(t)*(enthalpyFluid(1)-enthalpyFluid(numRows));
% 		else
% 			logEnergyExtracted(t) = timeStep*mdotAir(t)*(enthalpyFluid(1)-enthalpyFluid(numRows));
% 		end
% 	else
% 		if chargeMode
% 			logEnergyInput(t) = timeStep*mdotAir(t)*(enthalpyFluid(1)-enthalpyFluid(numRows));
% 		else
% 			logEnergyExtracted(t) = timeStep*mdotAir(t)*(enthalpyFluid(1)-enthalpyFluid(numRows));
% 		end
% 	end
    
	if useHX
		if chargeMode
			intEnergyFluid(2,1) = getEnthalpy(tempFluid(t+1,1));
			intEnergyFluid(2,numRows) = intEnergyFluid(2,numRows-1);
		else
			intEnergyFluid(2,numRows) = getEnthalpy(tempFluid(t+1,numRows));
		end
	else
		if chargeMode
			intEnergyFluid(2,1) = getEnthalpy(tempFluid(t+1,1));
			intEnergyFluid(2,numRows) = intEnergyFluid(2,numRows-1);
		else
			intEnergyFluid(2,numRows) = getEnthalpy(tempFluid(t+1,numRows));
		end
	end
	
    intEnergyStorage(2,1) = intEnergyStorage(2,2);
    intEnergyStorage(2,numRows) = intEnergyStorage(2,numRows-1);
    
%     if includePCM
%     logEnergyStored(t) = (sum(intEnergyStorage(2:rowCutoffPCM,2)) - sum(intEnergyStorage(2:rowCutoffPCM,1)))*densityPCMEff*rowVolume(i)*voidFractionInverse(i) ...
%         + (sum(intEnergyStorage(rowStartSHS:numRows-1,2)) - sum(intEnergyStorage(rowStartSHS:numRows-1,1)))*densityRocks*rowVolume(i)*voidFractionInverse(i);
%     else
%         logEnergyStored(t) = (sum(intEnergyStorage(2:numRows-1,2)) - sum(intEnergyStorage(2:numRows-1,1)))*densityRocks*rowVolume(i)*voidFractionInverse(i);
%     end
    try             % Calculate temperatures or throw error
        for i = 1:1:rowCutoffPCM
            if tempStorage(t,i) > TmeltEnd
                cPCM(1,i) = cEffPCMLiquid;
            elseif tempStorage(t,i) > TmeltStart
                cPCM(1,i) = cEffPCMPhaseChange;
            else
                cPCM(1,i) = cEffPCMSolid;
            end
		end
        % storage temps don't depend on int energy vs enthalpy
        if includePCM
            tempStorage(t+1,1:rowCutoffPCM) = calculateTempPCMMatrix(tempStorage(t,1:rowCutoffPCM),intEnergyStorage(2,1:rowCutoffPCM),intEnergyStorage(1,1:rowCutoffPCM),cPCM(1,1:rowCutoffPCM));
		end
        tempStorage(t+1,rowStartSHS:numRows) = calculateTempSHS(tempStorage(t,rowStartSHS:numRows),intEnergyStorage(2,rowStartSHS:numRows),intEnergyStorage(1,rowStartSHS:numRows),getSpecificHeatSHS);
        tempFluid(t+1,:) = calculateTempsFluid(intEnergyFluid(2,:),intEnergyFluid(1,:), tempFluid(t,:), specificHeatFluid,getEnthalpy,getSpecificHeat);
	catch MEx
        disp(t)
        disp(MEx.identifier)
        disp(MEx.message)
        fprintf(['Ended at time ' datestr(clock,'yyyy/mm/dd--HH:MM:SS') '\n'])
        return
	end
	
	% HX temp for discharge
    if ~chargeMode && useHX
        tempHXCO2(t,2) = tempProfile(false);
        enthalpyHXCO2(t,2) = getEnthalpyCO2HP(tempHXCO2(t,2)); 
%         q(t) = effectivenessHX*CminCharge*(tempFluid(t,1)-tempProfile(false));
        enthalpyHXCO2(t,1) = enthalpyHXCO2(t,2) + q(t)/(mdotCO2(t));

        tempHXCO2(t,1) = calculateTempFluidNR(enthalpyHXCO2(t,1),false, chargeTemp, 0,@getEnthalpyCO2HP,@getCpCO2HP); 
	end
    
	if useHX
		exergy(t,1) = timeStep*mdotCO2(t)*getExergy(tempHXCO2(t,1),ambientTemp,@getEnthalpyCO2HP,@getEntropyCO2HP);
		exergy(t,2) = timeStep*mdotCO2(t)*getExergy(tempHXCO2(t,2),ambientTemp,@getEnthalpyCO2HP,@getEntropyCO2HP);
		exergy(t,3) = exergy(t,1)-exergy(t,2);
        logExergyInputRecirculating(t) = exergy(t,1)-exergy(t,2);
		logExergyInputTopTempOnly(t) = exergy(t,1);
% 		timeStep*mdotCO2(t)*getExergy(tempHXCO2(t,1),ambientTemp,@getEnthalpyCO2HP,@getEntropyCO2HP);

%         logEnergyInput(t) = timeStep*mdotCO2(t)*(getEnthalpyCO2HP(tempHXCO2(t,1)) - getEnthalpyCO2HP(tempHXCO2(t,2)));
        logEnergyInput(t) = timeStep*mdotCO2(t)*(enthalpyHXCO2(t,1) - enthalpyHXCO2(t,2));
	else
        logExergyInputRecirculating(t) = timeStep*mdotAir(t)*(getExergy(tempFluid(t,1),ambientTemp,getEnthalpy,getEntropy) ...
            - getExergy(tempFluid(t,numRows),ambientTemp,getEnthalpy,getEntropy));
		logExergyInputTopTempOnly(t) = timeStep*mdotAir(t)*getExergy(tempFluid(t,1),ambientTemp,getEnthalpy,getEntropy);
		
        logEnergyInput(t) = timeStep*mdotAir(t)*(enthalpyFluid(1,1) - enthalpyFluid(1,numRows));
    end
	
    logFluidTempChangeActual(t) = tempFluid(t+1,10) - tempFluid(t,10);
    if calculatePressureDrop
        if includePCM   % Pressure drop calculations
            pressureDrop(t) = pressureDrop(t) + getPressureDrop(rowLengthPCM,dEncap, voidFraction(2:rowCutoffPCM), voidFractionInverse(2:rowCutoffPCM), tempFluid(t+1,2:rowCutoffPCM)-tempFluid(t,2:rowCutoffPCM), tempFluid(t,2:rowCutoffPCM),viscosityFluid(2:rowCutoffPCM),densityFluid(2:rowCutoffPCM),G);
        end
        if includeRock
            pressureDrop(t) = pressureDrop(t) + getPressureDrop(rowLengthSHS,rockAverageDiameter, voidFraction(rowStartSHS:numRows-1), voidFractionInverse(rowStartSHS:numRows-1), tempFluid(t+1,rowStartSHS:numRows-1)-tempFluid(t,rowStartSHS:numRows-1), tempFluid(t,rowStartSHS:numRows-1),viscosityFluid(rowStartSHS:numRows-1),densityFluid(rowStartSHS:numRows-1),G);
        end
    end
    
%     intEnergyMatrixStorage(:,1) = intEnergyMatrixStorage(:,2);
%     intEnergyMatrixFluid(:,1) = intEnergyMatrixFluid(:,2);
        
    if mod(t,progressBarPercentTime) == 0
        clc
        timeSoFar = toc;
        fprintf('Progress: %g %%\n',round(t/totalTime*100,0))
		if chargeMode
			fprintf('Charge mode. sCO2 chargeTemp: %g / threshold %g\n',tempHXCO2(t,2),chargeFluidTempThreshold );
		else
			fprintf('Discharge mode. sCO2 dischargeTemp: %g / threshold %g\n',tempHXCO2(t,1),dischargeFluidTempThreshold);
		end
        fprintf('Time elapsed so far: %g min\n',timeSoFar/60)
        fprintf('Estimated time to completion: %g min\n',timeSoFar*(totalTime/t-1)/60)
        if t < 2*onePercentTime && timeSoFar<5
            progressBarPercentTime = 10*onePercentTime;
        end
    end
    if t == logTime
        fprintf('Updated logEnergyFlow at t=%g\n',t)
        if breakAfterLogTime
            break
        end
    end
    t=t+1;
end

disp('Transient analysis complete')
transientTime = toc;
fprintf('The transient analysis took %g minutes.\n',transientTime/60)

% [includeConduction,includeRadiation,includeBulkHeatLoss,includePCMWallExchange,includeCapRadiation]
if breakAfterLogTime
	return
end
%% Energy logger calculations
% {
% {
    energyInputKWH = sum(logEnergyInput(1:chargeTime))/3.6e6;
    exergyInputKWHRecirculating = sum(logExergyInputRecirculating(1:chargeTime))/3.6e6;
	
    exergyInputKWHVented = sum(logExergyInputTopTempOnly(1:chargeTime))/3.6e6;
	
%     if strcmp(conditionSelection,'Zanganeh2015') || strcmp(conditionSelection,'Trahan')
	if totalTime > chargeTime+1
		if includePCM
			energySHSEndDischarge = (intEnergyStorage(2,rowStartSHS:numRows-1));
			energyPCMEndDischarge =  (intEnergyStorage(2,2:rowCutoffPCM));
			energyStoredSHS = sum((energySHSEndCharge - energySHSStart).*densityStorage(rowStartSHS:numRows-1).*rowVolume(rowStartSHS:numRows-1).*voidFractionInverse(rowStartSHS:numRows-1));
			energyRemovedSHS = sum((energySHSEndCharge - energySHSEndDischarge).*densityStorage(rowStartSHS:numRows-1).*rowVolume(rowStartSHS:numRows-1).*voidFractionInverse(rowStartSHS:numRows-1));
		else
			energySHSEndDischarge = (intEnergyStorage(2,2:numRows-1));
			energyStoredSHS = sum((energySHSEndCharge - energySHSStart).*densityStorage(2:numRows-1).*rowVolume(2:numRows-1).*voidFractionInverse(2:numRows-1));
			energyRemovedSHS = sum((energySHSEndCharge - energySHSEndDischarge).*densityStorage(2:numRows-1).*rowVolume(2:numRows-1).*voidFractionInverse(2:numRows-1));
		end
		
disp('hi')
		energyExtractedKWH = sum(logEnergyInput(chargeTime+1:totalTime))/3.6e6;
		exergyOutputKWHRecirculating = sum(logExergyInputRecirculating(chargeTime+1:totalTime))/3.6e6;
		exergyOutputKWHVented = sum(logExergyInputTopTempOnly(chargeTime+1:totalTime))/3.6e6;
	else % prevents code from breaking if breakAfterLogTime=true
		totalEnergySHSEndCharge = sum(intEnergyStorage(1,3:numRows-2));
		energyStoredSHS = sum((energySHSEndCharge - energySHSStart).*densityStorage(rowStartSHS:numRows-1).*rowVolume(rowStartSHS:numRows-1).*voidFractionInverse(rowStartSHS:numRows-1));
		energyRemovedSHS = 0;
		exergyOutputKWHRecirculating = 0;
	end
	energyStoredSHSKWH = energyStoredSHS/3.6e6;
	energyRemovedSHSKWH = energyRemovedSHS/3.6e6;
	energyWallLossChargingKWH = sum(logFluidHeatLoss(1:chargeTime))/3.6e6;
	energyWallLossDischargingKWH = sum(logFluidHeatLoss(chargeTime+1:totalTime))/3.6e6;
%     else
%         energyStoredSHSKWH = 0;
%         energyWallLossRocksChargingKWH =sum(logFluidHeatLoss(1:chargeTime))/3.6e6;
%         energyWallLossRocksDischargingKWH=0;
%     end

    if includePCM
%         exist totalEnergyPCMEndCharge
        if ~exist('energyPCMEndDischarge','var')
			energyPCMEndCharge = (intEnergyStorage(1,2:rowCutoffPCM));
            energyPCMEndDischarge =  0; %sum(intEnergyStorage(3:rowCutoffPCM,1));
        end

        energyStoredPCM = sum((energyPCMEndCharge - energyPCMStart).*densityStorage(2:rowCutoffPCM).*rowVolume(2:rowCutoffPCM).*voidFractionInverse(2:rowCutoffPCM));
        energyRemovedPCM = sum((energyPCMEndCharge - energyPCMEndDischarge).*densityStorage(2:rowCutoffPCM).*rowVolume(2:rowCutoffPCM).*voidFractionInverse(2:rowCutoffPCM));
        energyStoredPCMKWH = energyStoredPCM/3.6e6;
        energyRemovedPCMKWH = energyRemovedPCM / 3.6e6;
        energyWallLossPCMChargingKWH = sum(logPCMLoss(1:chargeTime))/3.6e6;
        energyWallLossPCMDischargingKWH = sum(logPCMLoss(chargeTime+1:totalTime))/3.6e6;
        energyWallPCMKWHCharge = sum(logPCMLoss(1:chargeTime))/3.6e6;
        energyWallPCMKWHDischarge = sum(logPCMLoss(chargeTime+1:totalTime))/3.6e6;
    else
        energyStoredPCMKWH = 0;
        energyWallLossPCMChargingKWH = 0;
        energyWallLossPCMDischargingKWH = 0;
        energyRemovedPCMKWH = 0;
    end
    if dischargeTimeHours == 0
        energyRemovedKWH = 0;
    end


	totalEnergyStoredKWH = energyStoredSHSKWH+energyStoredPCMKWH;
    energyBalanceChargingKWH = energyInputKWH-energyStoredSHSKWH-energyStoredPCMKWH+energyWallLossChargingKWH;
    if dischargeTimeHours > 0
        energyBalanceDischargingKWH = -energyExtractedKWH+energyRemovedSHSKWH+energyRemovedPCMKWH+energyWallLossDischargingKWH;
    else
        energyBalanceDischargingKWH = 0;
    end
	utilizationRatio = totalEnergyStoredKWH/designEnergyThermalKWH*100;

%% Written Output
clc
fprintf('PCM volume fraction: %g\n',pcmFraction)
fprintf('Total energy input: %g MWh\n',energyInputKWH/1000)
fprintf('Utilization ratio: %g %%\n',utilizationRatio)
fprintf('Total energy stored in SHS: %g MWh\n',energyStoredSHSKWH/1000)
fprintf('Total energy stored in PCM: %g MWh\n',energyStoredPCMKWH/1000)
fprintf('Energy lost through walls during charging: %g MWh\n',energyWallLossChargingKWH/1000)
fprintf('Energy lost through walls during discharging: %g MWh\n',energyWallLossDischargingKWH/1000)
fprintf('Overall energy balance during charging: %g MWh = %g%% of input -- target 0\n',energyBalanceChargingKWH/1000,energyBalanceChargingKWH/(energyInputKWH)*100)
fprintf('Overall energy balance during discharging: %g MWh = %g%% of input -- target 0\n',energyBalanceDischargingKWH/1000,energyBalanceDischargingKWH/(energyRemovedSHSKWH+energyRemovedPCMKWH)*100)
fprintf('Total energy output: %g MWh\n',energyExtractedKWH/1000)
fprintf('Total energy removed from SHS: %g MWh\n',energyRemovedSHSKWH/1000)
fprintf('Total energy removed from PCM: %g MWh\n',energyRemovedPCMKWH/1000)
fprintf('Run energy efficiency: %g %%\n',(energyExtractedKWH)/(energyInputKWH)*100)
fprintf('Total exergy input: %g MWh\n',exergyInputKWHRecirculating/1000)
fprintf('Total exergy output: %g MWh\n',exergyOutputKWHRecirculating/1000)
fprintf('Run exergy efficiency: %g %%\n',(exergyOutputKWHRecirculating)/(exergyInputKWHRecirculating)*100)
fprintf('Run exergy efficiency assuming venting: %g %%\n',(exergyOutputKWHVented)/(exergyInputKWHVented)*100)



energyFlowOutputs = [energyInputKWH, energyExtractedKWH;
	energyStoredSHSKWH, energyRemovedSHSKWH;
	energyStoredPCMKWH, energyRemovedPCMKWH;
    energyWallLossChargingKWH,energyWallLossDischargingKWH;
	energyBalanceChargingKWH, energyBalanceDischargingKWH;
	((energyExtractedKWH)/(energyInputKWH)*100), 0];

    
% disp([energyInputKWH; (energyStoredPCMKWH+energyWallLossRocksChargingKWH)]')


%}

%% Excel dump and memory check before plotting
fileNameTime = clock;

mem = memory;
if mem.MemUsedMATLAB > 12e9
    fprintf('\nMemory used by matlab is %g gigabytes, stopping to minimize risk of computer crash.',mem.MemUsedMATLAB/1e9)
%     return
end

%% Excel output
% {
if dumpToExcel
    secondsPerRow = 60;
    for t = 1:1:floor(totalTime/timeStepInverse/secondsPerRow)
        logTempFluid(1,t) = tempHXCO2(t*timeStepInverse*secondsPerRow,1);
        logTempFluid(2,t) = tempHXCO2(t*timeStepInverse*secondsPerRow,2);
        logTempFluid(3:numRows+2,t) = tempFluid(t*timeStepInverse*secondsPerRow,:);
        logTempStorage(:,t) = tempStorage(t*timeStepInverse*secondsPerRow,:);
        logPressure(t) = pressureDrop(t*timeStepInverse*secondsPerRow);
        logmdotAir(t) = mdotAir(t*timeStepInverse*secondsPerRow);
        logmdotCO2(t) = mdotCO2(t*timeStepInverse*secondsPerRow);
        logvdot(t) = vdot(t*timeStepInverse*secondsPerRow);
		exergyLog(1:3,t) = [exergy(t*timeStepInverse*secondsPerRow,1),exergy(t*timeStepInverse*secondsPerRow,2),exergy(t*timeStepInverse*secondsPerRow,1)-exergy(t*timeStepInverse*secondsPerRow,2)];
	end
	
	parametersDescriptions = {'Condition: ';'Fluid: '; ...
        'PCM name'; ...
        'SHS name'; ...
        'PCM fraction'; ...
        'Tank diameter (m)'; ...
        'Tank height (m)'; ...
        'Design energy (MWh-th)'; ...
        'Ein (MWh)'; ...
        'Eout (MWh)'; ...
        'E stored SHS (MWh)'; ...
        'E stored PCM (MWh)'; ...
        'E stored total (MWh)'; ...
        'Wall loss charging (MWh)'; ...
        'Wall loss discharging (MWh)'; ...
        'Exergy input (MWh)';
        'Exergy output (MWh)';
        'Energy efficiency';
        'Exergy efficiency (%)';
        'Utilization ratio';
		'Transient Time (s)';
		'Charge energy balance %';
		'Discharge energy balance %';
		'HX effectiveness';
		'enthalpy fusion modifier';
		'specific heat modifier';
        };
	
	parametersValues = {conditionSelection; tankFluidName; ...
		pcmChoice; ...
        shsChoice; ...
        pcmFraction; ...
        tankDiameter; ...
        tankHeight; ...
        designEnergyThermalGWH*1000; ...
        energyInputKWH/1000; ...
        energyExtractedKWH/1000; ...
        energyStoredSHSKWH/1000; ...
        energyStoredPCMKWH/1000; ...
        (energyStoredSHSKWH+energyStoredPCMKWH)/1000; ...
        energyWallLossChargingKWH/1000; ...
        energyWallLossDischargingKWH/1000; ...
        exergyInputKWHRecirculating/1000;
        exergyOutputKWHRecirculating/1000;
        (energyExtractedKWH)/(energyInputKWH)*100;
        exergyOutputKWHRecirculating/exergyInputKWHRecirculating*100;
        energyInputKWH/designEnergyThermalKWH*100;
		transientTime; 
		energyBalanceChargingKWH/(energyInputKWH)*100;
		energyBalanceDischargingKWH/(energyRemovedSHSKWH+energyRemovedPCMKWH)*100;
		effectivenessHX;
		enthalpyFusionModifier;
		specificHeatModifier;
        };
	
	
% 	prevRunTemps = xlsread(['.\Outputs\Excel Dumps\' pcmChoice '--X=' num2str(pcmFraction) '-- cycle' num2str(cycleNumber-1) '.xlsx'],'TankTempsAtEnd');
	if tempControlledTime
		if pcmFraction > 0
			workSpaceName = ['.\Outputs\Excel Dumps\Tm=' num2str(Tmelt) ' -- X=' num2str(pcmFraction) ' -- cycle' num2str(cycleNumber) '.xlsx'];
		else
			workSpaceName = ['.\Outputs\Excel Dumps\rockOnly -- cycle ' num2str(cycleNumber) '.xlsx'];
		end
	else
		if pcmFraction > 0
			workSpaceName = ['.\Outputs\Excel Dumps\' pcmChoice '_X=' num2str(pcmFraction) '_' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.xlsx'];
		else
			workSpaceName = ['.\Outputs\Excel Dumps\rockOnly_' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.xlsx'];
		end
	end
	
	% writematrix(A,workSpaceName)
    % type(workSpaceName)
    warning('off','MATLAB:xlswrite:AddSheet');
    try
        disp('Attempting to write temperature histories to file...')
        
        % logEnergyFlow
        xlswrite(workSpaceName,logEnergyFlow,'logEnergyFlow','A1:M72');
        disp('Wrote logEnergyFlow to file')
        
        % parameters
        xlswrite(workSpaceName,parametersDescriptions,'Parameters','A1:A26');
        xlswrite(workSpaceName,parametersValues,'Parameters','B1:B26');
        disp('Wrote model parameters to file')
        
        % Main tank, both PCM and rock cases
		mainTankOutput = [logTempFluid',logTempStorage',logmdotAir',logmdotCO2',exergyLog'];
        if useHX
			mainTankOutputLabels = {'sCO2 high temp','sCO2 low temp','fluid then storage then mdot'};
		else
			mainTankOutputLabels = {'Fluid top','Fluid bottom','fR1','fR2','fR3','fR4','fR5','fR6','fR7','fR8','fR9',...
				's1','sR1','sR2','sR3','sR4','sR5','sR6','sR7','sR8','sR9','sN','mdot'};
        end
        xlswrite(workSpaceName,mainTankOutputLabels,'Tank Temps','A1:C1');
        xlswrite(workSpaceName,mainTankOutput,'Tank Temps',['A2:DG' num2str(1+length(logTempFluid))]);
        xlswrite(workSpaceName,[tempStorage(chargeTime,:);tempFluid(chargeTime,:)]','Tank temps at chargeTime',['A2:B' num2str(1+numRows)]);
		xlswrite(workSpaceName,[tempStorage(totalTime,:);tempFluid(totalTime,:)]','TankTempsAtEnd',['A2:B' num2str(1+numRows)]);
		disp('Wrote tank temperatures to file')
     

        
        disp(['---> Successfully wrote temperature histories to file: ' workSpaceName])
		clear mainTankOutput


    catch
        disp('Error writing to file')
	end
end
%}

if plotFigures
%% Pressure Drop Plot
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
%% Mass and volume flow rate Plots
% {
if plotmdot
figure
hold on

ylim([0 0.02]);
% uistack(h,'bottom');

% plot(mdot2,'LineWidth',6)
plot(mdotAir,'LineWidth',2)
plot(mdotCO2,'LineWidth',2)
if compareToValidation
    plot(timeStepInverse*3600*validationDataRocks(:,21),validationDataRocks(:,22),':','LineWidth',2,'Color',[0.4940 0.1840 0.5560])
end
ylabel('Mass flow (kg/s)')


% yyaxis right
% plot(vdot*2118.88,'LineWidth',2)
% ylabel('Volume flow (cfm)')
% ylim([0 30]);
% 
% hold off
    
xticks(1*timeStepInverse*3600*[1 2 3 4 5 6 7 8])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
%ylabel('Temperature (C)')

legend('mdot','vdot')

xlim(3600*timeStepInverse*[0 8])
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

if saveFigures
    % saveas(gcf,figureName)
    figureNameFlowRatePNG = ['..\Model Outputs\FlowRate\flowRate' datestr(fileNameTime,'mmm-dd-yyyy HH-MM-SS') '.png'];
    saveas(fig,figureNameFlowRatePNG)
    fprintf('Wrote to file: %s\n',figureNameFlowRatePNG)
end
if closeFiguresAfterPlotting
	close gcf
end
end
%}
%% Fluid temperature in tank Plot
if plotFluid
figure
hold on

if useHX
	plot(tempHXCO2(plotTimeStart:plotTimeEnd,1),'LineWidth',2,'Color',color0)
else
	plot(tempFluid(plotTimeStart:plotTimeEnd,1),'LineWidth',2,'Color',color0)
end
plot(tempFluid(plotTimeStart:plotTimeEnd,R1),'LineWidth',2,'Color',color1)
plot(tempFluid(plotTimeStart:plotTimeEnd,R2),'LineWidth',2,'Color',color2)
plot(tempFluid(plotTimeStart:plotTimeEnd,R3),'LineWidth',2,'Color',color3)
plot(tempFluid(plotTimeStart:plotTimeEnd,R4),'LineWidth',2,'Color',color4)
plot(tempFluid(plotTimeStart:plotTimeEnd,R5),'LineWidth',2,'Color',color5)
plot(tempFluid(plotTimeStart:plotTimeEnd,R6),'LineWidth',2,'Color',color6)
plot(tempFluid(plotTimeStart:plotTimeEnd,R7),'LineWidth',2,'Color',color7)
plot(tempFluid(plotTimeStart:plotTimeEnd,R8),'LineWidth',2,'Color',color8)
plot(tempFluid(plotTimeStart:plotTimeEnd,R9),'LineWidth',2,'Color',color9)
if useHX
	plot(tempHXCO2(plotTimeStart:plotTimeEnd,2),'LineWidth',2,'Color',color0)
else
	plot(tempFluid(plotTimeStart:plotTimeEnd,numRows),'LineWidth',2,'Color',color0)
end
    
if includePCM
plot(TmeltStart*ones(1,plotTimeEnd),'LineWidth',0.5)
end

ylim(ylimTemps);
ylabel(['Temperature (' char(176) 'C)'])

hold off


xlabel('Time (hours)')
xticks(1*timeStepInverse*3600*[1:20])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'})
%ylabel('Temperature (C)')

yticks([400:50:700])

if useHX
legend({'sCO2 high','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','sCO2 low','Tmelt'},'Location','south','NumColumns',2)
else
legend({'Air tank top','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','Air tank bottom','Tmelt'},'Location','south','NumColumns',2)
end
title(['Fluid temperature in tank - PCM Fraction = ' num2str(pcmFraction) ' Tmelt = ' num2str(TmeltStart+dTMelt/2)])

   
xlim(3600*timeStepInverse*xlimHours)
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

if saveFigures
    figureNameFluidPNG = ['..\Model Outputs\Fluid\fluidTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
    saveas(fig,figureNameFluidPNG)
    fprintf('Wrote to file: %s\n',figureNameFluidPNG)
end
if closeFiguresAfterPlotting
% 	close gcf
end
end
%}
%% Fluid Individual Row temperatures in tank Plot
if strcmp(conditionSelection,'TrahanPCM')
    
figure
hold on
plot((tempFluid(plotTimeStart:plotTimeEnd,R2)+tempFluid(plotTimeStart:plotTimeEnd,R2a))/2,'LineWidth',2,'Color',color3)
uppererrorR2 = validationDataRocks(:,2) - validationDataRocks(:,8);
lowererrorR2 = validationDataRocks(:,14) - validationDataRocks(:,2);
errorbar(timeStepInverse*60*validationDataRocks(:,1),validationDataRocks(:,2),lowererrorR2,uppererrorR2,'o','LineWidth',1,'Color',color3)
ylim(ylimTemps);
hold off
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time (min)')
xticks(timeStepInverse*60*[0:20:220])
xticklabels({'0' '20' '40' '60' '80' '100' '120' '140' '160' '180' '200' '220'})
title('Fluid temperature R2')
xlim(3600*timeStepInverse*xlimHours)
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

if saveFigures
    figureNameFluidPNG = ['..\Model Outputs\Fluid\fluidTempTrahanPCMR2-' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
    saveas(fig,figureNameFluidPNG)
    fprintf('Wrote to file: %s\n',figureNameFluidPNG)
end

figure
hold on
plot((tempFluid(plotTimeStart:plotTimeEnd,R3)+tempFluid(plotTimeStart:plotTimeEnd,R3a))/2,'LineWidth',2,'Color',color4)
uppererrorR3 = validationDataRocks(:,4) - validationDataRocks(:,10);
lowererrorR3 = validationDataRocks(:,16) - validationDataRocks(:,4);
errorbar(timeStepInverse*60*validationDataRocks(:,3),validationDataRocks(:,4),lowererrorR3,uppererrorR3,'o','LineWidth',1,'Color',color4)
ylim(ylimTemps);
hold off
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time (min)')
xticks(timeStepInverse*60*[0:20:220])
xticklabels({'0' '20' '40' '60' '80' '100' '120' '140' '160' '180' '200' '220'})
title('Fluid temperature R3')
xlim(3600*timeStepInverse*xlimHours)
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

if saveFigures
    figureNameFluidPNG = ['..\Model Outputs\Fluid\fluidTempTrahanPCMR3-' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
    saveas(fig,figureNameFluidPNG)
    fprintf('Wrote to file: %s\n',figureNameFluidPNG)
end      

figure
hold on
plot((tempFluid(plotTimeStart:plotTimeEnd,R4)+tempFluid(plotTimeStart:plotTimeEnd,R4a))/2,'LineWidth',2,'Color',color5)
uppererrorR4 = validationDataRocks(:,6) - validationDataRocks(:,12);
lowererrorR4 = validationDataRocks(:,18) - validationDataRocks(:,6);
errorbar(timeStepInverse*60*validationDataRocks(:,5),validationDataRocks(:,6),lowererrorR4,uppererrorR4,'o','LineWidth',1,'Color',color5)
ylim(ylimTemps);
hold off
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time (min)')
xticks(timeStepInverse*60*[0:20:220])
xticklabels({'0' '20' '40' '60' '80' '100' '120' '140' '160' '180' '200' '220'})
title('Fluid temperature R4')
xlim(3600*timeStepInverse*xlimHours)
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

if saveFigures
    figureNameFluidPNG = ['..\Model Outputs\Fluid\fluidTempTrahanPCMR4-' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
    saveas(fig,figureNameFluidPNG)
    fprintf('Wrote to file: %s\n',figureNameFluidPNG)
end
if closeFiguresAfterPlotting
	close gcf
end
end
%}
%% Storage temperature in tank
% {
if plotStorage
figure
hold on


plot(tempStorage(plotTimeStart:plotTimeEnd,2),'LineWidth',2,'Color',color1)
plot(tempStorage(plotTimeStart:plotTimeEnd,R1),'LineWidth',2,'Color',color1)
plot(tempStorage(plotTimeStart:plotTimeEnd,R2),'LineWidth',2,'Color',color2)
plot(tempStorage(plotTimeStart:plotTimeEnd,R3),'LineWidth',2,'Color',color3)
plot(tempStorage(plotTimeStart:plotTimeEnd,R4),'LineWidth',2,'Color',color4)
plot(tempStorage(plotTimeStart:plotTimeEnd,R5),'LineWidth',2,'Color',color5)
plot(tempStorage(plotTimeStart:plotTimeEnd,R6),'LineWidth',2,'Color',color6)
plot(tempStorage(plotTimeStart:plotTimeEnd,R7),'LineWidth',2,'Color',color7)
plot(tempStorage(plotTimeStart:plotTimeEnd,R8),'LineWidth',2,'Color',color8)
plot(tempStorage(plotTimeStart:plotTimeEnd,R9),'LineWidth',2,'Color',color9)
plot(tempStorage(plotTimeStart:plotTimeEnd,numRows),'LineWidth',2,'Color',color9)
    
if includePCM
plot(TmeltStart*ones(1,plotTimeEnd),'LineWidth',2)
end

ylim(ylimTemps);
ylabel(['Temperature (' char(176) 'C)'])

hold off


xlabel('Time (hours)')
xticks(1*timeStepInverse*3600*[1:20])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20'})
%ylabel('Temperature (C)')

legend({'Top','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','bottom','Tmelt'},'Location','south','NumColumns',2)

title(['Storage temperature in tank - PCM Fraction = ' num2str(pcmFraction) ' Tmelt = ' num2str(TmeltStart+dTMelt/2)])

   
xlim(3600*timeStepInverse*xlimHours)
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

if saveFigures
    figureNameStoragePNG = ['..\Model Outputs\Storage\storageTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
    saveas(fig,figureNameStoragePNG)
    fprintf('Wrote to file: %s\n',figureNameStoragePNG)
end
if closeFiguresAfterPlotting
	close gcf
end
end
%}
%}
%% Plot temp gradient at each time sequentially
% {
% transVis = true; %transVis;
% transVis = input('Figure loaded. Press 1 to run transient visualization, or any other number to skip: ');
if transVis

    figure
    tic
    startRow = 2;

    endRow = numRows-1;
    for t = 1:timeStepInverse*60*0.5:totalTime
        % % f
%         t= logTime

        clf
%         figure
        hold on
        yyaxis right
        ylim([0 100])
        plot(abs(tempFluid(t,startRow:endRow)-tempStorage(t,startRow:endRow)),'LineWidth',3)

        yyaxis left
        ylim(ylimTemps)
        plot(tempFluid(t,startRow:endRow),'LineWidth',3)
        plot(tempStorage(t,startRow:endRow),'LineWidth',3)
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
%}

end
%% notify of ending
fprintf(['\nEnded at time ' datestr(clock,'yyyy/mm/dd--HH:MM:SS') '\n'])