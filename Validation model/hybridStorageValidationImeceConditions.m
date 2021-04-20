%% top
clc
clear all
close all
addpath(genpath('fluid properties'))

useTopLogTime = false;
% useTopLogTime = true;

breakAfterLogTime = false;
% breakAfterLogTime = true; % false

logTime = 3; %floor(chargeTime/2); 

conditionSelection = 'Zanganeh2015'; % options are Zanganeh2015, Trahan, TrahanPCM, Beasley1989-low, Beasley1989-high, IMECE...
rockChoice = 'Zanganeh2015';    % 'NIST' 'Smetana2013' 'Zanganeh2015' 'Trahan'
fluidTempAlgorithm = 'cp'; % options are 'cp' or 'nr'
rockTempAlgorithm = 'quadratic'; % options are 'quadratic' or 'cp' or 'nr'

includePCM = true;

includeConduction = true;
includeRadiation = true;
includeBulkHeatLoss = true;
includePCMWallExchange = true;
includeCapRadiation = true;

useInternalEnergy = false;

fluidName = 'DryAirAtmo';           % options are CO2HP, CO2LP, DryAir, DryAirAtmo, MoistAir
fluidMethod = 'REFPROP';             % options are Trahan or refprop. Doesn't seem to make a difference :)
inputMode = 1;                      % options are 1 (mass flow) or 2 (volume flow). NYI
referenceFluid = 'DryAirAtmo';      % When using CO2HP for fluid at controlled vdot, ref fluid should be DryAir
                                    % for controlled mdot, ref fluid should be same as fluidName
bypassFraction = 0;%.15;           % 0.15 used in Zanganeh 2015 but unclear whether it was applied before or after figure plot. Seems to be applied before.

increaseFluidTemps = false;         % Always true (via override later) for CO2, 
increaseBaseTemp = false;           % can be true or false for air

runTransient = true;
conductionModifier = 1;
constantFluidProperties = false;

if strcmp(fluidName,'CO2HP') || strcmp(fluidName,'CO2LP')
    load('inletTempProfile40Hz3hrcharge75Cmin.mat')
end
useRockTempProfileForPCM = false;

%% Plot selector etc -- moved to individual conditions
transVis= false;        transVisMdot = false;
saveFigures = false;     dumpToExcel = true;
compareToValidation = false;
printEnergyFlows = true;
printLNorms = false;

%% Condition set
switch conditionSelection
    case 'Zanganeh2015'
		imeceTempIncrease = 0;
        zanganeh = true;
        includeRocks = true;
        nonsphericalPCM = true;
        avgTemp = (25+650)/2;
        plotAverage = false;     
        if includePCM plotPCM = false; else plotPCM = false; end
        plotFluid = false;      plotStorage = false;    plotmdot = false;
        calculatePressureDrop = false;
        useTrahanTempProfile = false;
        trackError = true;
        trackFluidError = true;
        
        load('inputsZanganeh2015.mat')
        chargeTimeHours = 3;            % 3 is standard
        chargeTimePartialHour = 1/6;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 5;             % 5 standard
        rowCountModifier = 1;
        numRows = rowCountModifier*68+2;
        
        if includePCM
%             tuned values for hvModifier= 1 are: 
%             plateTempIncreaseCharge = 49.94759434;  plateTempIncreaseDischarge = -7.187420888;
%             tuned values for hvModifier= 0.75 are: 
%             plateTempIncreaseCharge = 51.28005047;   plateTempIncreaseDischarge = -7.124672876;

            hvModifier = 1; %0.75;% 0.3;
            
            plateTempIncreaseCharge = 51.68;
            plateTempIncreaseDischarge = -7.22;
            plateTempIncrease = plateTempIncreaseCharge;
			wallTempIncreaseCharge=zeros(1,numRows); wallTempIncreaseDischarge=zeros(1,numRows);
			wallTempIncreaseCharge(1,2:5) = 9.28;
			wallTempIncreaseDischarge(1,2:5) = 2.02;
			
			wallTempIncrease = wallTempIncreaseCharge;
        else
            hvModifier = 1; %0.75;% 0.3;
            
            plateTempIncreaseCharge = 41.9;
            plateTempIncreaseDischarge = 17.7;
            plateTempIncrease = plateTempIncreaseCharge;
            wallTempIncreaseCharge = zeros(1,numRows);        % static
            wallTempIncreaseDischarge = zeros(1,numRows);
        end
        validationrows3plus = true;
		
    case 'Trahan'
        voidFractionPCMInverse = 1;
        emissivityEncapsulation = 0;
        zanganeh = false;
        includeRocks = true;
        nonsphericalPCM = false;
        plotAverage = false;    plotPCM = false;
        plotFluid = true;      plotStorage = true;    plotmdot = false;
        calculatePressureDrop = false;
        useTrahanTempProfile = true;
        trackError = false;
        trackFluidError = true;
        
        load('inputsTrahan.mat')
        includePCM = false;
        chargeTimeHours = 1.5;             % 3 is standard
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 1.5;             % 5 standard
        rowCountModifier = 1;
        numRows = rowCountModifier*100+2;
        avgTemp = (25+500)/2;
        
        plateTempIncreaseCharge = 0;
        plateTempIncreaseDischarge = 0;
        wallTempIncreaseCharge = [0 0 4 3 2 1]*0;        % static
        wallTempIncreaseDischarge = [0 0 4 3 2 1]*0;
        
        hvModifier = 1; %0.75;% 0.3;
        
        validationrows3plus = true;
    case 'TrahanPCM'
        includeRocks = false;
        nonsphericalPCM = false;
        emissivityEncapsulation = 0;
        zanganeh = false;
        plotAverage = false;    plotPCM = false;
        plotFluid = true;      plotStorage = false;    plotmdot = false;
        calculatePressureDrop = false;
        useTrahanTempProfile = false;
        trackError = false;
        trackFluidError = true;
        
        load('inputsTrahanPCM.mat')
        includePCM = true;
        chargeTimeHours = 2;             % only 2 hours really necessary
%         chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 3;             % 5 standard
        rowCountModifier = 1;
        numRows = rowCountModifier*50+2;
        
        plateTempIncreaseCharge = 0;
        plateTempIncreaseDischarge = 0;
        wallTempIncreaseCharge = [1:numRows]*0;        % static
        wallTempIncreaseDischarge = [1:numRows]*0;
        
        hvModifier = 1; %0.75;% 0.3;
        
        validationrows3plus = true;
    case 'Beasley1989-low'
        includeRadiation = false;
        includeRocks = false;
        nonsphericalPCM = false;
        zanganeh = false;
        avgTemp = (25+52)/2;
        plotAverage = false;    plotPCM = false;
        plotFluid = false;      plotStorage = true;    plotmdot = false;
        calculatePressureDrop = false;
        useTrahanTempProfile = false;
        trackFluidError = false;
        trackError = false;
        
        load('inputsBeasley1989.mat')
        plotAverage = false;
        plotStorage = true;
        hvModifier = 1;
        rowCountModifier = 1;
        numRows = rowCountModifier*40+2;
                
        chargeTimeHours = 3;                % 3 is standard\
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 0;             % 5 standard
        includePCM = true;
        includePCMWallExchange = false;
        plateTempIncreaseCharge = 0;
        plateTempIncreaseDischarge = 0;
        wallTempIncreaseCharge = zeros(1,numRows);        % static
        wallTempIncreaseDischarge = zeros(1,numRows);
        validationrows3plus = false;
    case 'Beasley1989-high'
        includeRadiation = false;
        includeRocks = false;
        zanganeh = false;
        nonsphericalPCM = false;
        avgTemp = (25+52)/2;
        plotAverage = false;    plotPCM = false;
        plotFluid = false;      plotStorage = true;    plotmdot = false;
        calculatePressureDrop = false;
        useTrahanTempProfile = false;
        trackFluidError = false;
        trackError = false;
        
        load('inputsBeasley1989.mat')
        hvModifier = 1;% 0.3;
        rowCountModifier = 1;
        numRows = rowCountModifier*20+2;
        
        chargeTimeHours = 2;                % 3 is standard\
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 0;             % 5 standard
        includePCM = true;
        includePCMWallExchange = false;
        plateTempIncreaseCharge = 0;
        plateTempIncreaseDischarge = 0;
        wallTempIncreaseCharge = zeros(1,numRows);        % static
        wallTempIncreaseDischarge = zeros(1,numRows);
        validationrows3plus = false;
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
    case 'nr'
        calculateTempFluid = @calculateTempFluidNR;
end
switch rockTempAlgorithm
    case 'quadratic'
        calculateTempRock = @calculateTempRockQuadratic;
%         (energy)(-747.0995+sqrt(747.0995^2+4*0.2838*energy))/2/0.2838;
    case 'cp'
        calculateTempRock = @calculateTempRockCP;
    case 'nr'
        calculateTempRock = @calculateTempRockNR;
end

%% ============== Time setup

if strcmp(fluidName,'DryAirAtmo')
    timeStepInverse = rowCountModifier*40; %50;                % 41        % samples per second
else
    timeStepInverse = 40;                % 41        % samples per second
end

switch conditionSelection
    case 'Trahan'
        if numRows == 102
            timeStepInverse = 145;
        else
            timeStepInverse = 100;
        end
    case 'TrahanPCM'
        timeStepInverse = 1000;
    case 'Beasley1989-high'
        timeStepInverse = rowCountModifier*200;
    case 'Beasley1989-low'
        timeStepInverse = rowCountModifier*125;
end

timeStep = 1/timeStepInverse;                       % seconds per sample
chargeTime = floor((chargeTimeHours+chargeTimePartialHour)*3600*timeStepInverse);        % 3+1/6 hours in seconds
dischargeTime = floor(dischargeTimeHours*3600*timeStepInverse);        % Normally 5h, shortened for quicker runtime while testing charge mode
numCycles = 1;                                      % #*(charge+discharge)
totalTime = floor(numCycles*(chargeTime+dischargeTime)) + 2;
onePercentTime = round(0.01*totalTime);
progressBarPercentTime = onePercentTime;
t = 1;

if ~useTopLogTime
    logTime = chargeTime-10; 
%     logTime = chargeTime+floor(dischargeTime/2); 
end

%% ============== Tank dimensions
switch conditionSelection
    case 'Zanganeh2015'
        if includePCM
            heightPCM = 0.09*1; %0.18; % 0.09
            rockHeight = tankHeight - heightPCM;
            numRowsPCM = round(4*heightPCM/0.09)*rowCountModifier;
            rowCutoffPCM = 1+numRowsPCM;
            rowStartRock = rowCutoffPCM+1; 

            numRowsRock = numRows-numRowsPCM-2;
            
            TmeltStart = 575;
            dTMelt = 4;
            TmeltEnd = TmeltStart + dTMelt;
%             ceffPCMSolid = ceffPCMSolid * 1.25;
%             ceffPCMPhaseChange = ceffPCMPhaseChange * 1.25;
            pcmLowEnergyCutoff = TmeltStart*ceffPCMSolid;
            pcmHighEnergyCutoff = pcmLowEnergyCutoff + ceffPCMPhaseChange*dTMelt;
            cPCM = zeros(1,rowCutoffPCM);
        else
            rowStartRock = 1; 
            rockHeight = tankHeight;      
            heightPCM = 0;
            numRowsPCM = 0;
            rowCutoffPCM = 2;
            numRowsRock = numRows-2;
        end
        beasley = false;
        TrahanPCM = false;
    case 'Trahan'
        rowStartRock = 1; 
        rockHeight = tankHeight;      
        heightPCM = 0;
        numRowsPCM = 0;
        rowCutoffPCM = 2;
        numRowsRock = numRows-4;
        beasley = false;
        TrahanPCM = false;
    case 'TrahanPCM'
        heightPCM = tankHeight;
        rockHeight = 0;
        numRowsPCM = (numRows-4);
        rowCutoffPCM = numRows;
        cPCM = zeros(1,rowCutoffPCM);
        rowStartRock = rowCutoffPCM+1; 
        numRowsRock = 0;
        beasley = false;
        TrahanPCM = true;
    case 'Beasley1989-high'
        heightPCM = tankHeight;
        rockHeight = 0;
        numRowsPCM = (numRows-4);
        rowCutoffPCM = numRows;
        cPCM = zeros(1,rowCutoffPCM);
        rowStartRock = rowCutoffPCM+1; 
        numRowsRock = 0;
        beasley = true;
        TrahanPCM = false;
    case 'Beasley1989-low'
        heightPCM = tankHeight;
        rockHeight = 0;
        numRowsPCM = (numRows-4);
        rowCutoffPCM = numRows;
        cPCM = zeros(1,rowCutoffPCM);
        rowStartRock = rowCutoffPCM+1; 
        numRowsRock = 0;
        beasley = true;
        TrahanPCM = false;
end
chargeMode = true;
firstRow = 1;
lastRow = numRows;
rowIncrement = 1;

%% ============== Discretization
rowLengthRock = rockHeight/numRowsRock;
rowVolumeRock = rowLengthRock*pi*(tankRadius)^2;
rowWallSurfaceAreaRock = pi*2*tankRadius*rowLengthRock;
if isnan(rowWallSurfaceAreaRock)
    rowWallSurfaceAreaRock = 0;
end

if includePCM
    rowLengthPCM = heightPCM/numRowsPCM;
    rowVolumePCM = rowLengthPCM*pi*(tankRadius)^2;
    rowWallSurfaceAreaPCM = pi*2*tankRadius*rowLengthPCM;
end

%% ============== Insulation parameters
switch conditionSelection
    case 'Zanganeh2015'
        % Adjust these values for actual conditions. Whatever that means?
        lossInsulationTermsbackup = lossInsulationTerms;
        lossInsulationTerms(1:26) = lossInsulationTermsbackup(1);
        lossInsulationTerms(27) = (lossInsulationTermsbackup(1)+lossInsulationTermsbackup(2))/2;
        lossInsulationTerms(28:47) = lossInsulationTermsbackup(2);
        lossInsulationTerms(48:numRows) = lossInsulationTermsbackup(3);
        clear lossInsulationTermsbackup
    case 'Trahan'
        tInsulationInside = 0.133; %m
        tInsulationOutside = 50.8/1000; %m
        tWall = 3.175/1000; %m
        
        kInsulation = 0.06; %W/mK
        kWall = 43; %W/mK based on 1% carbon steel
        
        lossInsulationTerms(1:numRows) = tankRadius * (1/kInsulation * log((tankRadius+tInsulationInside)/(tankRadius)) + ...
            1/kWall * log((tankRadius+tInsulationInside+tWall)/(tankRadius+tInsulationInside)) + ...
            1/kInsulation * log((tankRadius+tInsulationInside+tWall+tInsulationOutside)/(tankRadius+tInsulationInside+tWall)));
    case 'TrahanPCM'
        lossInsulationTerms(1:numRows) = tankRadius * (1/kInsulation * log((tankRadius+tInsulationInside)/(tankRadius)) + ...
            1/kWall * log((tankRadius+tInsulationInside+tWall)/(tankRadius+tInsulationInside)) + ...
            1/kInsulation * log((tankRadius+tInsulationInside+tWall+tInsulationOutside)/(tankRadius+tInsulationInside+tWall)));
        
    otherwise
        lossInsulationTerms(1:numRows) = 0;
end


%% ============== Rock properties
switch conditionSelection
    case 'Zanganeh2015'
        rocksPerRow = rowVolumeRock*voidFractionRocksInverse / rockAverageVolume;
        rowSurfaceAreaRocks = rocksPerRow*rockAverageSurfaceArea;   % m2/row
        getkRocks = @(x) 5+(1-5)*(x-25)/(700-25);					%@(x)5+(1-5)*(x-25)/(700-25);
    case 'Trahan'
        rocksPerRow = rowVolumeRock*voidFractionRocksInverse / rockAverageVolume;
        rowSurfaceAreaRocks = rocksPerRow*rockAverageSurfaceArea;   % m2/row
        getkRocks = @(x) 1.5;
    case 'TrahanPCM'
        rocksPerRow = rowVolumeRock*voidFractionRocksInverse / rockAverageVolume;
        rowSurfaceAreaRocks = rocksPerRow*rockAverageSurfaceArea;   % m2/row
        getkRocks = @(x) 0.52;
end

%% ============== Metal Plate properties

%% Fluid selection
switch fluidName
    case 'CO2LP'        
        fluidPressure=12.38e6;
        disp('CO2 properties set to Low Pressure Mode (P=12.38 MPa = 12.38e6 Pa).')
        getEnthalpy = @getEnthalpyCO2;
        getInternalEnergy = @getInternalEnergyCO2;
        getSpecificHeat = @getCpCO2;
        getConductivity = @getConductivityCO2;
        getViscosity = @getViscosityCO2;
        getDensity = @getDensityCO2;
        if useInternalEnergy
            getEnergyDerivative = @getCvCO2LP;
        else
            getEnergyDerivative = @getCpCO2LP;
        end
        useCO2 = true; 
        increaseTemps = true;
        fluidTempIncreaseAmount = 50;
    case 'CO2HP'
        fluidPressure=25e6;
        disp('CO2 properties set to High Pressure Mode (P=25 MPa = 25e6 Pa).')
        getEnthalpy = @getEnthalpyCO2HP;
        getInternalEnergy = @getInternalEnergyCO2HP;
        getSpecificHeat = @getCpCO2HP;
        getConductivity = @getConductivityCO2HP;
        getViscosity = @getViscosityCO2HP;
        getDensity = @getDensityCO2HP;
        if useInternalEnergy
            getEnergyDerivative = @getCvCO2HP;
        else
            getEnergyDerivative = @getCpCO2HP;
        end
        useCO2 = true; 
        increaseTemps = true;
        fluidTempIncreaseAmount = 50;
    case 'DryAir'
        fluidPressure = 25e6;
        getEnthalpy = @getEnthalpyDryAir;
        getInternalEnergy = @getInternalEnergyDryAir;
        getSpecificHeat = @getCpDryAir;
        getConductivity = @getConductivityDryAir;
        getViscosity = @getViscosityDryAir;
        getDensity = @getDensityDryAir;
        if useInternalEnergy
            getEnergyDerivative = @getCvDryAir;
        else
            getEnergyDerivative = @getCpDryAir;
        end
        useCO2 = false;
    case 'DryAirAtmo'
        fluidPressure = 101.3e3;
        if strcmp(conditionSelection,'Trahan') || strcmp(conditionSelection,'TrahanPCM') || strcmp(fluidMethod,'Trahan')
            getDensity = @(Tave) ((-5.75399E-16)*(Tave.^5))+((3.02846E-12)*(Tave.^4))-((6.18352E-9)*(Tave.^3))+((6.29927E-6)*(Tave.^2))-((3.5422E-3)*Tave)+1.25079; %density of air (kg/m3) 
            getViscosity = @(Tave) (((6.10504E-10)*(Tave.^3))-((2.13036E-6)*(Tave.^2))+((4.71398E-3)*(Tave))+1.67555)*(10.^-5); %Dynamic viscosity of air (kg/ms) 
            getSpecificHeat = @(Tave) (((1.28806E-13)*(Tave.^4))-((4.46054E-10)*(Tave.^3))+((4.8772E-7)*(Tave.^2))+((1.82754E-5)*Tave)+1.00651)*1000; %Cp of air (J/kg-K) 
            getConductivity = @(Tave) ((-4.44955E-15)*(Tave.^4))+((2.41702E-11)*(Tave.^3))-((4.09601E-8)*(Tave.^2))+((7.91034E-5)*Tave)+.0242006; % thermal conductivity of air (W/mK)
            getEnthalpy = @(Tave) (((1.28806E-13)*(Tave.^5)/5)-((4.46054E-10)*(Tave.^4)/4)+((4.8772E-7)*(Tave.^3)/3)+((1.82754E-5)*Tave.^2)/2+1.00651*Tave)*1000; %Cp of air (J/kg-K) 
        else
            getEnthalpy = @getEnthalpyDryAirAtmo;
            getInternalEnergy = @getInternalEnergyDryAirAtmo;
            getConductivity = @getConductivityDryAirAtmo;
            getViscosity = @getViscosityDryAirAtmo;
            getDensity = @getDensityDryAirAtmo;
            getSpecificHeat = @getCpDryAirAtmo;
            if useInternalEnergy
                getEnergyDerivative = @getCvDryAirAtmo;
            else
                getEnergyDerivative = @getCpDryAirAtmo;
            end
        end
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

if constantFluidProperties
    getEnthalpy = @(x) getSpecificHeat(1)*x;
%     fluidIntEnergy = getInternalEnergy(avgTemp);
%     getInternalEnergy = @(x) fluidIntEnergy;
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
%% PCM/encapsulation Properties

switch conditionSelection
    case 'Zanganeh2015'
        if includePCM
            %calculateTemps = @calculateTempsWithPCM;
            tempProfile = @tempProfileZanganeh;
            getmdot = @getmdotWithPCM;
            getStorageTopConduction = @getPCMConductionZanganeh;
            %getPCMConduction = @getPCMConductionZanganeh;
            storageWallHXareaTerms = rowWallSurfaceAreaPCM*fcontWall / (tankRadius-thicknessEncapsulation) / log(tankRadius/(tankRadius-thicknessEncapsulation));
        else
            pcmLowEnergyCutoff = 0;
            pcmHighEnergyCutoff = 0;
            TmeltStart = 0;
            TmeltEnd = 0;
            ceffPCMSolid = 0;
            ceffPCMPhaseChange = 0;
            ceffPCMLiquid = 0;

           % calculateTemps = @calculateTempsNoPCM;
            tempProfile = @tempProfileZanganeh;
            getmdot = @getmdotNoPCM;
            
            %getPCMConduction = @getPCMConductionZanganeh;
            getStorageTopConduction = @getRockConduction;
        end
        calculatehvPCM = @calculateHvPCMZanganeh;
    case 'Trahan'
        tempProfile = @tempProfileTrahan;
        getmdot = @getmdotTrahan;
        pcmLowEnergyCutoff = 0;
        pcmHighEnergyCutoff = 0;
        TmeltStart = 0;
        TmeltEnd = 0;
        ceffPCMSolid = 0;
        ceffPCMPhaseChange = 0;
        ceffPCMLiquid = 0;

        getStorageTopConduction = @getRockConduction;
    case 'TrahanPCM'
%         tempProfile = @tempProfileTrahan;
%         getmdot = @getmdotTrahan;
%         pcmLowEnergyCutoff = 0;
%         pcmHighEnergyCutoff = 0;
%         TmeltStart = 0;
%         TmeltEnd = 0;
%         ceffPCMSolid = 0;
%         ceffPCMPhaseChange = 0;
%         ceffPCMLiquid = 0;

        calculatehvPCM = @calculateHvTrahan;
    case 'Beasley1989-high'
            %calculateTemps = @calculateTempsNoRocks;
            tempProfile = @tempProfileBeasleyHigh;
            getmdot = @getmdotBeasleyHigh;
            getStorageTopConduction = @getPCMConductionBeasley;
            calculatehvPCM = @calculateHvTrahan; %@calculateHvPCMBeasley;
            %getPCMConduction = @getPCMConductionBeasley;
            wallTempIncreaseChargeVariable = 0;
            storageWallHXareaTerms = 0;
    case 'Beasley1989-low'
            %calculateTemps = @calculateTempsNoRocks;
            tempProfile = @tempProfileBeasleyLow;
            getmdot = @getmdotBeasleyLow;
            getStorageTopConduction = @getPCMConductionBeasley;
            calculatehvPCM = @calculateHvTrahan; %@calculateHvPCMBeasley;
            %getPCMConduction = @getPCMConductionBeasley;
            wallTempIncreaseChargeVariable = 0;
            storageWallHXareaTerms = 0;
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
    case 'Trahan'
        dischargeTemp = 192+baseTempIncreaseAmount;                                 % C - temperature of HTF at bottom inlet during discharge
        baseTemp = 25+baseTempIncreaseAmount;                                      % C - starting temperature of rock bed
        ambientTemp = 25;
    case 'TrahanPCM'
        baseTemp = 286;
    case 'Beasley1989-high'
        baseTemp = 25; %26.5;
        ambientTemp = 25;
    otherwise
        baseTemp = 26.5; %26.5;
        ambientTemp = 25;
end


% Avoid computer crashes
% if totalTime > 5e6
%     disp('Array size too large. Shutting down to prevent computer crash.')
%     return
% end

%% ============== Temp/energy/step/void fraction Matrix setup

tempStorage = ones(totalTime,numRows);    % C
tempStorage(1,:) = baseTemp;    % C
intEnergyStorage = zeros(2,numRows);
% fluidEnergyChange = zeros(numRows,1);
% storageEnergyChange = zeros(numRows,1);

% rock internal energy, J/kg
% intEnergyRock = 747.0995*tempRock+0.2838*tempRock.^2;       % J
switch conditionSelection
    case 'Trahan'
        % Jamie value
        % raising bed from 25 to 500 C takes input of 60.8 kWh, far higher
        % than the 36 kWh specified by Jamie's dissertation
%         getRockEnergy = @(temp) 1.2953375e-06*(temp+273).^4 - 0.002526691666667*(temp+273).^3 + 2.942740725*(temp+273).^2 - 429.234632*(temp+273);
%         getRockSpecificHeat = @(temp) 5.18135e-6*(temp+273).^3 - 7.580075e-3*(temp+273).^2 + 5.88548145*(temp+273) - 429.234632;
        
        % NIST value
%         getRockSpecificHeat = @(temp) 
%         getRockEnergy = @(temp) 

        % Zanganeh2015 value
        switch rockChoice
            case 'Trahan'
                getRockEnergy = @(temp) 1.2953375e-06*(temp+273).^4 - 0.002526691666667*(temp+273).^3 + 2.942740725*(temp+273).^2 - 429.234632*(temp+273);
                getRockSpecificHeat = @(temp) 5.18135e-6*(temp+273).^3 - 7.580075e-3*(temp+273).^2 + 5.88548145*(temp+273) - 429.234632;
            case 'NIST'
%                 ((temp+273)/1000) = t
                getRockSpecificHeat = @(temp) 6.262211312058515*(93.43834+108.3577*((temp+273)/1000)-50.86447*((temp+273)/1000).^2+25.58683*((temp+273)/1000).^3-1.61133*((temp+273)/1000).^(-2));
%                 ((temp+273)/1000)
                getRockEnergy = @(temp) (749.8980157*temp) + 10090488.95*((temp + 273).^(-1)) + (0.270235027*temp.^2) - (0.0000624319*temp.^3) + (0.0000000400575*temp.^4);
            case 'Zanganeh2015'
                getRockEnergy = @(temp) 747.0995*temp+0.2838*temp.^2;
                getRockSpecificHeat = @(temp) 747.0995+2*0.2838*temp;
            case 'Smetana2013'
                getRockSpecificHeat = @(temp) 2.42E-11*(temp+273).^4 - 6.71E-08*(temp+273).^3 + 6.96E-05*(temp+273).^2 - 3.15E-02*(temp+273) + 5.78E+00;
                getRockEnergy = @(temp) 2.42E-11/5*(temp+273).^5 - 6.71E-08/4*(temp+273).^4 + 6.96E-05/3*(temp+273).^3 - 3.15E-02/2*(temp+273)^2 + 5.78E+00*(temp);
        end
    case 'TrahanPCM'
        getRockEnergy = @(temp) 0;
        getRockSpecificHeat = @(temp) 0;
    case 'Zanganeh2015'
        getRockEnergy = @(temp) 747.0995*temp+0.2838*temp.^2;
        getRockSpecificHeat = @(temp) 747.0995+2*0.2838*temp;
    case 'Beasley1989-high'
        getRockEnergy = @(temp) 0;
        getRockSpecificHeat = @(temp) 0;
    case 'Beasley1989-low'
        getRockEnergy = @(temp) 0;
        getRockSpecificHeat = @(temp) 0;
end

if strcmp('Zanganeh2015',conditionSelection) || strcmp('Trahan',conditionSelection)
    if includePCM
        intEnergyStorage(1,1:rowCutoffPCM) = ceffPCMSolid * tempStorage(1,1:rowCutoffPCM);
        intEnergyStorage(1,rowStartRock:numRows) = getRockEnergy(tempStorage(1,rowStartRock:numRows));
        totalEnergyRocksStart = sum(intEnergyStorage(1,rowStartRock:numRows-1));
        totalEnergyPCMStart = sum(intEnergyStorage(1,2:rowCutoffPCM));
        emissivityStorage(1:rowCutoffPCM) = emissivityEncapsulation;
        emissivityStorage(rowStartRock:numRows) = emissivityRocks;
        voidFraction(1:rowCutoffPCM) = voidFractionPCM;
        voidFraction(rowStartRock:numRows) = voidFractionRocks;
        voidFractionInverse(:) = 1-voidFraction(:);
        rowVolume(1:rowCutoffPCM) = rowVolumePCM;
        rowVolume(rowStartRock:numRows) = rowVolumeRock;
        densityStorage(1:rowCutoffPCM) = densityPCMEff;
        densityStorage(rowStartRock:numRows) = densityRocks;
        rowLength(1:rowCutoffPCM) = rowLengthPCM;
        rowLength(rowStartRock:numRows) = rowLengthRock;
        
        for i = 1:rowCutoffPCM
            getStorageConduction{i} = @getPCMConductionZanganeh;
            calculateTempStorage{i} = @calculateTempPCM;
        end
        for i = rowStartRock:numRows
            getStorageConduction{i} = @getRockConduction;
            calculateTempStorage{i} = calculateTempRock;
        end
    else
        intEnergyStorage(1,1:numRows) = getRockEnergy(tempStorage(1,1:numRows));
        totalEnergyRocksStart = sum(intEnergyStorage(1,2:numRows-1));
        emissivityStorage(1:numRows) = emissivityRocks;
        voidFraction(:) = voidFractionRocks;
        voidFraction(1:numRows) = voidFractionRocks;
        voidFractionInverse(:) = 1-voidFraction(:);
        rowVolume(1:numRows) = rowVolumeRock;
        densityStorage(1:numRows) = densityRocks;
        rowLength(1:numRows) = rowLengthRock;
        for i = 1:numRows
            getStorageConduction{i} = @getRockConduction;
            calculateTempStorage{i} = calculateTempRock;
        end
    end
else
    
    emissivityStorage(1:numRows) = 0.7;
    intEnergyStorage(1,1:rowCutoffPCM) = ceffPCMSolid * tempStorage(1,:);
    totalEnergyPCMStart = sum(intEnergyStorage(1,2:numRows-1));
    voidFraction(1:numRows) = voidFractionPCM;
    voidFractionInverse(1:numRows) = 1-voidFraction(:);
    rowVolume(1:numRows) = rowVolumePCM;
    densityStorage(1:numRows) = densityPCMEff;
    rowLength(1:numRows) = rowLengthPCM;
    for i = 1:numRows
        getStorageConduction{i} = @getPCMConductionBeasley;
        calculateTempStorage{i} = @calculateTempPCM;
    end
end
rowMassStorage(1:numRows) = rowVolume(:).*densityStorage(:).*voidFractionInverse(:);
totalMass = sum(rowMassStorage(2:numRows-1));

storageRadiationPlateCoefficient = timeStep*5.67e-8*plateArea*voidFractionPlateInverse/(1/emissivityStorage(2)+1/emissivityPlate-1);
storageRadiationPCMCoefficient = timeStep*5.67e-8*tankArea*voidFractionPCMInverse/(2/emissivityEncapsulation-1);
storageRadiationInterfaceCoefficient = timeStep*5.67e-8*tankArea*voidFractionPCMInverse/(1/emissivityEncapsulation+1/emissivityRocks-1);

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

% fluid matrices
tempFluid = ones(totalTime,numRows);    % C
tempFluid(1,:) = baseTemp;    % C
intEnergyFluid = zeros(2,numRows);
enthalpyFluid = zeros(1,numRows);
logenthalpyFluid = zeros(totalTime,1);
logIntEnergyFluid = zeros(totalTime,1);

enthalpyFluid(1,:) = getEnthalpy(tempFluid(1,:));
if useInternalEnergy
    intEnergyFluid(1,1:numRows) = getInternalEnergy(tempFluid(1,1:numRows));
else
    intEnergyFluid(1,1:numRows) = getEnthalpy(tempFluid(1,1:numRows));
end

intEnergyStorage(2,:) = intEnergyStorage(1,:);
intEnergyFluid(2,:) = intEnergyFluid(1,:);

%% ============== Fluid Base properties and flow rates
% densityFluidDischarge = getDensity(dischargeTemp);
    densityFluid(numRows) = 0;
    densityFluidPrevious(numRows) = 0;
    viscosityFluid(numRows) = 0;
    kFluid(numRows) = 0;
    specificHeatFluid(numRows) = 0;

for i = 1:numRows
    densityFluid(i) = getDensity(tempFluid(t,i));
    densityFluidPrevious = densityFluid;
    viscosityFluid(i) = getViscosity(tempFluid(t,i));
    kFluid(i) = getConductivity(tempFluid(t,i));
    specificHeatFluid(i) = getSpecificHeat(tempFluid(t,i));
end

mdot = zeros(1,totalTime);
vdot = zeros(1,totalTime);
bypassFractionInverse = 1-bypassFraction;               % bypass fraction 
mdot(1)= 0.008*bypassFractionInverse;

G = mdot(1)/tankArea;

%% Stability calculations
switch conditionSelection
    case 'Zanganeh2015'
        mdotmax = 0.014;
        tempMax = 650;
        VFmin = 0.4;
        if includePCM
            minrowlength = min(rowLengthRock,rowLengthPCM);
            VFmax = 0.549;
        else
            minrowlength = rowLengthRock;
            VFmax = 0.4;
        end
        rhomin = getDensity(tempMax);
        rhominsolid = densityRocks;
        kminsolid = 1;
        specificheatminsolid = 788.6150;
        
    case 'Trahan'
        mdotmax = 0.04484;
        tempMax = 500;
        VFmin = 0.51;
        minrowlength = rowLengthRock;
        VFmax = 0.4;
        rhomin = getDensity(tempMax);
        rhominsolid = densityRocks;
        kminsolid = 1.5;
        specificheatminsolid = 747;
    case 'TrahanPCM'
        mdotmax = 0.05;
        tempMax = 326;
        VFmin = 0.348;
        VFmax = VFmin;
        minrowlength = rowLengthPCM;
        rhomin = getDensity(tempMax);
        rhominsolid = densityPCMEff;
        kminsolid = 0.50;
        specificheatminsolid = 1655;
    case 'Beasley1989-high'
        mdotmax = 0.0328;
        tempMax = 60;
        rhomin = getDensity(tempMax);
        rhominsolid = 812;
        VFmax = 0.369;
        VFmin = 0.369;
        kminsolid = 0.24;
        specificheatminsolid = 2093;
        minrowlength = min(rowLengthRock,rowLengthPCM);
    case 'Beasley1989-low'
        mdotmax = 0.008792;
        tempMax = 60;
        rhomin = getDensity(tempMax);
        rhominsolid = 812;
        VFmin = 0.369;
        VFmax = 0.369;
        kminsolid = 0.24;
        specificheatminsolid = 2093;
        minrowlength = min(rowLengthRock,rowLengthPCM);
end
% vmax = mdotmax / VFmin/tankArea/rhomin
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
%     return
end
if timeStep >= maxdTsolid
    disp('solid stability check failed. Quitting.')
%     return
end

%% ============== Energy and heat transfer variables preallocation

fluidAdvection = zeros(1,numRows);
fluidHeatLoss = zeros(1,numRows);

convectionExchange = zeros(1,numRows);
storageConduction = zeros(1,numRows);
storageRadiation = zeros(1,numRows);

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
pressureDrop(totalTime) = 0;
energyIn = 0;
logConvectionStorage(totalTime) = 0;
logConduction(totalTime) = 0;
logConvectionFluid(totalTime) = 0;
logRadiationPlate(totalTime) = 0;
logRockLoss(totalTime) = 0;
logPCMLoss(totalTime) = 0; 
logEnergyInput(totalTime) = 0;
logEnergyExtracted(totalTime) = 0;
logEnergyStored(totalTime) = 0;
logQRadPlate(totalTime) = 0;
logQWallPCM(totalTime) = 0;
QWallPCM = zeros(1,numRows);
storagePlate  = zeros(1,numRows);

logFluidEnergyChange(totalTime) = 0;
logFluidTempChangeTheoretical(totalTime) = 0;
logFluidTempChangeActual(totalTime) = 0;

%% Figure setup and excel dump
close all
% Row selection based on PCM use
switch conditionSelection
    case 'Zanganeh2015'
        if includePCM
            R1 = 1+19;     R2 = 1+36;    R3 = 1+52;    R4 = 1+67; R0 = 1;
            legendEntry0 = '0 mm';
            legendEntry1 = '379 mm';
            legendEntry2 = '720 mm';
            legendEntry3 = '1045 mm';
            legendEntry4 = '1345 mm';
            legendEntry5 = ['tank top ('  fluidName  ')'];
            legendEntry6 = [fluidName ' mass flow (-->)'];
        %     I = imread('..\Model Outputs\Source Figures\Fig12aAverageTempsPCMpng');
        else
            R1 = 11+1;    R2 = 29+1;    R3 = 43+1;    R4 = 67+1;  R0 = 1;
            legendEntry0 = '0 mm';
            legendEntry1 = '225 mm';
            legendEntry2 = '574 mm';
            legendEntry3 = '865 mm';
            legendEntry4 = '1335 mm';
            legendEntry5 = ['tank top ('  fluidName  ')'];
            legendEntry6 = [fluidName ' mass flow (-->)'];
        %     I = imread('..\Model Outputs\Source Figures\Fig11aAverageTempsRocksOnly.png');
        end
        xlimHours = [0 8];
        ylimFlow = [0 0.02];
        ylimTemps = [25 675];
    case 'Trahan'
        % row numbers TBD
        if numRows == 70
            R0 = 2;
            R1 = 3;     R1a = 3; 
            R2 = 13;    R2a = 14;
            R3 = 31;    R3a = 32;
            R4 = 48;    R4a = 49;
            R5 = 66;    R5a = 67;
        elseif numRows == 102
            R0 = 2;
            R1 = 3;    R1a = 3;
            R2 = 19;   R2a = 19;
            R3 = 44;   R3a = 46;
            R4 = 70;   R4a = 71;
            R5 = 96;   R5a = 97;
        end
        legendEntry0 = '0 mm';
        legendEntry1 = '225 mm';
        legendEntry2 = '574 mm';
        legendEntry3 = '865 mm';
        legendEntry4 = '1335 mm';
        legendEntry5 = ['tank top ('  fluidName  ')'];
        legendEntry6 = [fluidName ' mass flow (-->)'];
        %     I = imread('..\Model Outputs\Source Figures\Fig11aAverageTempsRocksOnly.png');
        xlimHours = [0 (chargeTimeHours+dischargeTimeHours)];
        ylimTemps = [20 520];
        ylimTempsDischarge = [150 550];
    case 'TrahanPCM'
        % row numbers TBD
        R0 = 2;
        R4 = (floor(0.2*(numRows-2)))+1;    R4a = R4;
        R3 = (floor(0.4*(numRows-2)))+1;    R3a = R3;
        R2 = (floor(0.6*(numRows-2)))+1;    R2a = R2;
        R1 = (floor(0.8*(numRows-2)))+1;    R1a = R1;
%         R5 = 50;    R5a = 50;
        legendEntry0 = '0 mm';
        legendEntry1 = '225 mm';
        legendEntry2 = '574 mm';
        legendEntry3 = '865 mm';
        legendEntry4 = '1335 mm';
        legendEntry5 = ['tank top ('  fluidName  ')'];
        legendEntry6 = [fluidName ' mass flow (-->)'];
        %     I = imread('..\Model Outputs\Source Figures\Fig11aAverageTempsRocksOnly.png');
        xlimHours = [0 3.6667];
        ylimTemps = [270 340];
%         ylimTempsDischarge = [150 550];
    case 'Beasley1989-low'
        R1 = round(0.242*(numRows-2))+1;
        R2 = round(0.579*(numRows-2))+1;
        R3 = 1;
        R4 = numRows;
        legendEntry1 = 'x/L = 0.25';
        legendEntry2 = 'x/L = 0.575';
        legendEntry3 = 'Air at x/L = 0.25';
        legendEntry4 = 'Air at x/L = 0.575';
        xlimHours = [0 3];
        ylimTemps = [25 65];
        ylimFlow = [0 0.04];
    case 'Beasley1989-high'
        R1 = round(0.242*(numRows-2))+1;
        R2 = round(0.494*(numRows-2))+1;
        R3 = 1;
        R4 = numRows;
        legendEntry1 = 'x/L = 0.25';
        legendEntry2 = 'x/L = 0.5';
        legendEntry3 = 'Air at x/L = 0.25';
        legendEntry4 = 'Air at x/L = 0.5';
        xlimHours = [0 2];
        ylimTemps = [25 65];
        ylimFlow = [0 0.04];
end

color1 = [0.4660 0.6740 0.1880];
color2 = [0 0.4470 0.7410];
color3 = [0.8500 0.3250 0.0980];
color4 = [0.9290 0.6940 0.1250];
color5 = [0.4940 0.1840 0.5560];
color6 = [0.750 0 0.25];
            
plotTimeStart = 1;
plotTimeEnd = totalTime;

%% Validation setup
if compareToValidation
    disp('Reading validation data file...')
    switch conditionSelection
        case 'Zanganeh2015'
            if includePCM
                validationDataSource = 'Data Files\Validation\ValidationPCMRedo.xlsx';
                validationDataRocks = xlsread(validationDataSource,'Rocks Section');
                validationDataPCM   = xlsread(validationDataSource,'PCM Section');
            else
                validationDataSource = 'Data Files\Validation\ValidationRocksOnly.xlsx';
                validationDataRocks = xlsread(validationDataSource,'Rocks Section');
            end
            row1validationModel=1; row2validationModel=1; row3validationModel=1; row4validationModel=1; fluidvalidationModel=1;  % used for error calc
            row1validationExp=1; row2validationExp=1; row3validationExp=1; row4validationExp=1; fluidvalidationExp=1;  % used for error calc
            PCM1validationModel=1; PCM3validationModel=1; fluid1validationModel=1; fluid3validationModel=1; % used for error calc
            PCM1validationExp=1; PCM3validationExp=1; fluid1validationExp=1; fluid3validationExp=1; % used for error calc
        case 'Trahan'
            validationDataSource = 'Data Files\Validation\ValidationTrahan.xlsx';
            validationDataRocks = xlsread(validationDataSource,'Exp 1 data');
            row1validationModel=1; row2validationModel=1; row3validationModel=1; row4validationModel=1; fluidvalidationModel=1;  % used for error calc
            row1validationExp=1; row2validationExp=1; row3validationExp=1; row4validationExp=1; fluidvalidationExp=1;  % used for error calc
            
            PCM1validationModel=1; PCM3validationModel=1; fluid1validationModel=1; fluid3validationModel=1; % Not used but kept to avoid breaking stuff
            PCM1validationExp=1; PCM3validationExp=1; fluid1validationExp=1; fluid3validationExp=1; % Not used but kept to avoid breaking stuff
        case 'TrahanPCM'
            validationDataSource = 'Data Files\Validation\ValidationTrahanPCM.xlsx';
            validationDataRocks = xlsread(validationDataSource,'LHS data');
            row1validationModel=1; row2validationModel=1; row3validationModel=1; row4validationModel=1; fluidvalidationModel=1;  % used for error calc
            row1validationExp=1; row2validationExp=1; row3validationExp=1; row4validationExp=1; fluidvalidationExp=1;  % used for error calc
            
            PCM1validationModel=1; PCM3validationModel=1; fluid1validationModel=1; fluid3validationModel=1; % Not used but kept to avoid breaking stuff
            PCM1validationExp=1; PCM3validationExp=1; fluid1validationExp=1; fluid3validationExp=1; % Not used but kept to avoid breaking stuff
        case 'Beasley1989-low'
            validationDataSource = 'Data Files\Validation\ValidationBeasleyLow.xlsx';
            validationDataRocks = xlsread(validationDataSource);    % Rocks to make code easier later; Beasley is pure PCM
                row1validationModel=1; row2validationModel=1;   % used for error calc
                row1validationExp=1; row2validationExp=1;   % used for error calc
        case 'Beasley1989-high'
            validationDataSource = 'Data Files\Validation\ValidationBeasleyHigh.xlsx';
            validationDataRocks = xlsread(validationDataSource);
                row1validationModel=1; row2validationModel=1;   % used for error calc
                row1validationExp=1; row2validationExp=1;   % used for error calc
    end
    disp('Done reading validation data')
    fluiderrorModel(length(validationDataRocks),2) = 0;
    row1errorModel(length(validationDataRocks),2) = 0;
    row2errorModel(length(validationDataRocks),2) = 0;
    row3errorModel(length(validationDataRocks),2) = 0;
    row4errorModel(length(validationDataRocks),2) = 0;
    
    fluiderrorExp(length(validationDataRocks),2) = 0;
    row1errorExp(length(validationDataRocks),2) = 0;
    row2errorExp(length(validationDataRocks),2) = 0;
    row3errorExp(length(validationDataRocks),2) = 0;
    row4errorExp(length(validationDataRocks),2) = 0;
    if includePCM && strcmp(conditionSelection,'Zanganeh2015')
        PCM1errorModel(length(validationDataPCM),2) = 0;
        PCM3errorModel(length(validationDataPCM),2) = 0;
        fluid1errorModel(length(validationDataPCM),2) = 0;
        fluid3errorModel(length(validationDataPCM),2) = 0;

        PCM1errorExp(length(validationDataPCM),2) = 0;
        PCM3errorExp(length(validationDataPCM),2) = 0;
        fluid1errorExp(length(validationDataPCM),2) = 0;
        fluid3errorExp(length(validationDataPCM),2) = 0;
    end
    
    
    errorRowMax = length(validationDataRocks)+1;
end

%% Memory check for computer crash prevention
mem = memory;
if mem.MemUsedMATLAB > 15e9
%     fprintf('\nMemory used by matlab is %g gigabytes, stopping to minimize risk of computer crash.',mem.MemUsedMATLAB/1e9)
%     return
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
while t < totalTime
    % Swap energy parameters to set up for new run
    intEnergyStorage(1,:) = intEnergyStorage(2,:);
    intEnergyFluid(1,:) = intEnergyFluid(2,:);
    
    if compareToValidation && trackError
        if beasley
            if row1validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row1validationModel,1)) < 1
                row1errorModel(row1validationModel,1) = validationDataRocks(row1validationModel,1);
                row1errorModel(row1validationModel,2) = tempStorage(t,R1) - validationDataRocks(row1validationModel,2);
                row1errorModel(row1validationModel,3) = tempStorage(t,R1);
                row1errorModel(row1validationModel,4) = validationDataRocks(row1validationModel,2);
                row1validationModel = row1validationModel + 1;
            end
%{             if row1validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row1validationExp,1)) < 1
%                 row1errorExp(row1validationExp,1) = validationDataRocks(row1validationExp,1);
%                 row1errorExp(row1validationExp,2) = tempStorage(R1,t) - validationDataRocks(row1validationExp,2);
%                 row1errorExp(row1validationExp,3) = tempStorage(R1,t);
%                 row1errorExp(row1validationExp,4) = validationDataRocks(row1validationExp,2);
%                 row1validationExp = row1validationExp + 1;
%             end
            if row2validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row2validationModel,3)) < 1
                row2errorModel(row2validationModel,1) = validationDataRocks(row2validationModel,3);
                row2errorModel(row2validationModel,2) = tempStorage(t,R2) - validationDataRocks(row2validationModel,4);
                row2validationModel = row2validationModel + 1;
            end
%             if row2validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row2validationExp,3)) < 1
%                 row2errorExp(row2validationExp,1) = validationDataRocks(row2validationExp,3);
%                 row2errorExp(row2validationExp,2) = tempStorage(R2,t) - validationDataRocks(row2validationExp,4);
%                 row2validationExp = row2validationExp + 1;
%             end
        else           
            if trackFluidError
                if fluidvalidationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(fluidvalidationModel,1)) < 1
                    fluiderrorModel(fluidvalidationModel,1) = validationDataRocks(fluidvalidationModel,1);
                    fluiderrorModel(fluidvalidationModel,2) = tempFluid(t,2) - validationDataRocks(fluidvalidationModel,2);
                    fluidvalidationModel = fluidvalidationModel + 1;
                end
                if fluidvalidationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(fluidvalidationExp,11)) < 1
                    fluiderrorExp(fluidvalidationExp,1) = validationDataRocks(fluidvalidationExp,11);
                    fluiderrorExp(fluidvalidationExp,2) = tempFluid(t,2) - validationDataRocks(fluidvalidationExp,12);
                    fluidvalidationExp = fluidvalidationExp + 1;
                end
            end
            if row1validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row1validationModel,3)) < 1
                row1errorModel(row1validationModel,1) = validationDataRocks(row1validationModel,3);
                row1errorModel(row1validationModel,2) = (tempFluid(t,R1)+tempStorage(t,R1))/2 - validationDataRocks(row1validationModel,4);
                row1validationModel = row1validationModel + 1;
            end
            if row1validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row1validationExp,13)) < 1
                row1errorExp(row1validationExp,1) = validationDataRocks(row1validationExp,13);
                row1errorExp(row1validationExp,2) = (tempFluid(t,R1)+tempStorage(t,R1))/2 - validationDataRocks(row1validationExp,14);
                row1validationExp = row1validationExp + 1;
            end
            if row2validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row2validationModel,5)) < 1
                row2errorModel(row2validationModel,1) = validationDataRocks(row2validationModel,5);
                row2errorModel(row2validationModel,2) = (tempFluid(t,R2)+tempStorage(t,R2))/2 - validationDataRocks(row2validationModel,6);
                row2validationModel = row2validationModel + 1;
            end
            if row2validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row2validationExp,15)) < 1
                row2errorExp(row2validationExp,1) = validationDataRocks(row2validationExp,15);
                row2errorExp(row2validationExp,2) = (tempFluid(t,R2)+tempStorage(t,R2))/2 - validationDataRocks(row2validationExp,16);
                row2validationExp = row2validationExp + 1;
            end
            if validationrows3plus
                if row3validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row3validationModel,7)) < 1
                    row3errorModel(row3validationModel,1) = validationDataRocks(row3validationModel,7);
                    row3errorModel(row3validationModel,2) = (tempFluid(t,R3)+tempStorage(t,R3))/2 - validationDataRocks(row3validationModel,8);
                    row3validationModel = row3validationModel + 1;
                end
                if row3validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row3validationExp,17)) < 1
                    row3errorExp(row3validationExp,1) = validationDataRocks(row3validationExp,17);
                    row3errorExp(row3validationExp,2) = (tempFluid(t,R3)+tempStorage(t,R3))/2 - validationDataRocks(row3validationExp,18);
                    row3validationExp = row3validationExp + 1;
                end
                if row4validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row4validationModel,9)) < 1
                    row4errorModel(row4validationModel,1) = validationDataRocks(row4validationModel,9);
                    row4errorModel(row4validationModel,2) = (tempFluid(t,R4)+tempStorage(t,R4))/2 - validationDataRocks(row4validationModel,10);
                    row4validationModel = row4validationModel + 1;
                end
                if row4validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataRocks(row4validationExp,19)) < 1
                    row4errorExp(row4validationExp,1) = validationDataRocks(row4validationExp,19);
                    row4errorExp(row4validationExp,2) = (tempFluid(t,R4)+tempStorage(t,R4))/2 - validationDataRocks(row4validationExp,20);
                    row4validationExp = row4validationExp + 1;
                end
            end
            
            % PCM error trackers need to be edited still: 
            % rows from source
            if includePCM
                if PCM1validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataPCM(PCM1validationModel,1)) < 1
                    PCM1errorModel(PCM1validationModel,1) = validationDataPCM(PCM1validationModel,1);
                    PCM1errorModel(PCM1validationModel,2) = tempStorage(t,3) - validationDataPCM(PCM1validationModel,2);
                    PCM1validationModel = PCM1validationModel + 1;
                end
                if PCM1validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataPCM(PCM1validationExp,9)) < 1
                    PCM1errorExp(PCM1validationExp,1) = validationDataPCM(PCM1validationExp,9);
                    PCM1errorExp(PCM1validationExp,2) = tempStorage(t,3) - validationDataPCM(PCM1validationExp,10);
                    PCM1validationExp = PCM1validationExp + 1;
                end
                if PCM3validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataPCM(PCM3validationModel,5)) < 1
                    PCM3errorModel(PCM3validationModel,1) = validationDataPCM(PCM3validationModel,5);
                    PCM3errorModel(PCM3validationModel,2) = tempStorage(t,5) - validationDataPCM(PCM3validationModel,6);
                    PCM3validationModel = PCM3validationModel + 1;
                end
                if PCM3validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataPCM(PCM3validationExp,13)) < 1
                    PCM3errorExp(PCM3validationExp,1) = validationDataPCM(PCM3validationExp,13);
                    PCM3errorExp(PCM3validationExp,2) = tempStorage(t,5) - validationDataPCM(PCM3validationExp,14);
                    PCM3validationExp = PCM3validationExp + 1;
                end
                
                if fluid1validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataPCM(fluid1validationModel,3)) < 1
                    fluid1errorModel(fluid1validationModel,1) = validationDataPCM(fluid1validationModel,3);
                    fluid1errorModel(fluid1validationModel,2) = tempFluid(t,3) - validationDataPCM(fluid1validationModel,4);
                    fluid1validationModel = fluid1validationModel + 1;
                end
                if fluid1validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataPCM(fluid1validationExp,11)) < 1
                    fluid1errorExp(fluid1validationExp,1) = validationDataPCM(fluid1validationExp,11);
                    fluid1errorExp(fluid1validationExp,2) = tempFluid(t,3) - validationDataPCM(fluid1validationExp,12);
                    fluid1validationExp = fluid1validationExp + 1;
                end
                if fluid3validationModel < errorRowMax && abs(t - timeStepInverse*3600*validationDataPCM(fluid3validationModel,7)) < 1
                    fluid3errorModel(fluid3validationModel,1) = validationDataPCM(fluid3validationModel,7);
                    fluid3errorModel(fluid3validationModel,2) = tempFluid(t,5) - validationDataPCM(fluid3validationModel,8);
                    fluid3validationModel = fluid3validationModel + 1;
                end
                if fluid3validationExp < errorRowMax && abs(t - timeStepInverse*3600*validationDataPCM(fluid3validationExp,15)) < 1
                    fluid3errorExp(fluid3validationExp,1) = validationDataPCM(fluid3validationExp,15);
                    fluid3errorExp(fluid3validationExp,2) = tempFluid(t,5) - validationDataPCM(fluid3validationExp,16);
                    fluid3validationExp = fluid3validationExp + 1;
                end
            end
        end
        
    end
    
    % Fluid properties for time step
    densityFluidPrevious = densityFluid;
    densityFluid(:) = getDensity(tempFluid(t,:));
    viscosityFluid(:) = getViscosity(tempFluid(t,:));
    kFluid(:) = getConductivity(tempFluid(t,:));
    specificHeatFluid(:) = getSpecificHeat(tempFluid(t,:));
            
    if useInternalEnergy
        enthalpyFluid = getEnthalpy(tempFluid(t,:));
    else
%         enthalpyFluid = intEnergyFluid(1,:);
        enthalpyFluid = getEnthalpy(tempFluid(t,:));
    end
    rowMassFluid(:) = rowVolume(:).*densityFluid(:).*voidFraction(:);
    logenthalpyFluid(t,1) = getEnthalpy(tempFluid(t,2));
    logIntEnergyFluid(t,1) = intEnergyFluid(1,2);
    
	
    % Storage conductivity
    if strcmp('Zanganeh2015',conditionSelection) && includePCM
        % [rowCutoffPCM rowStartRock]
        kStorage (rowStartRock:numRows) = getkRocks(tempStorage(t,rowStartRock:numRows));
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
            kStorage(i) = getkRocks(tempStorage(t,i));
        end
    end
    
    % Dispersion and loss coefficients
    % Add a Trahan version of this function since her correlations might be
    % different than Zanganeh's
   [keffRocks,heffWall] = calculatekeff(kStorage,kFluid,tempFluid(t,:),tempStorage(t,:),voidFraction,rockAverageDiameter,numRows,emissivityRocks);
      
   % Switch to discharge mode, log energy before proceeding
    if mod(t,(chargeTime+dischargeTime)) == chargeTime+1 || (t == totalTime && chargeMode)
        chargeMode = false;
        if includePCM % log energy at end of charge before switching to discharge
            totalEnergyRocksEndCharge = sum(intEnergyStorage(1,rowStartRock:numRows-1));
            storageEnergyChargeEndLog = intEnergyStorage(1,:);
            if strcmp(conditionSelection,'TrahanPCM')
                totalEnergyPCMEndCharge = sum(intEnergyStorage(1,2:numRows-1));
            else
                totalEnergyPCMEndCharge = sum(intEnergyStorage(1,2:rowCutoffPCM));
            end
        else
            totalEnergyRocksEndCharge = sum(intEnergyStorage(1,2:numRows-1));
            storageEnergyChargeEndLog = intEnergyStorage(1,:);
        end
        % logInternalEnergy = intEnergyMatrix(:,2);
        enthalpyCoefficient = enthalpyCoefficientDischarge;
        
        plateTempIncrease = plateTempIncreaseDischarge;
        wallTempIncrease = wallTempIncreaseDischarge;
        disp('Switching to discharge mode')
    % Switch back to charge mode for extra cycles
    elseif t > chargeTime && mod(t,(chargeTime+dischargeTime)) == 1 
        enthalpyCoefficient = enthalpyCoefficientCharge;
        
        chargeMode = true;
        plateTempIncrease = plateTempIncreaseCharge;
        wallTempIncrease = wallTempIncreaseCharge;
        disp('Switching to charge mode')
        disp(t)
    end
    
%     fprintf('===timestep = %g===\n',t)
    mdot(t) = bypassFractionInverse*getmdot(t,chargeTime,timeStep);
    G=mdot(t)/tankArea;
    hvRocks = hvRocksCoeff*G^0.92;  % W/m2K
    
    % Inlet fluid temp depends on flow direction
    if chargeMode
        if strcmp(fluidName,'DryAirAtmo')
            tempFluid(t+1,1) = tempProfile(t,timeStep,0,tempFluid(t,1),includePCM,1,includeDischarge);
% 			fprintf('Fluid temp from profile function: %g\n',tempFluid(t+1,1))
        else
            tempFluid(t+1,1) = inletTempProfile(1,t+1);
        end
        if useInternalEnergy
            intEnergyFluid(1,1) = getInternalEnergy(tempFluid(t,1));
        else
% 			fprintf('Fluid temp used for HT calcs inlet enthalpy: %g\n',tempFluid(t,1))
            intEnergyFluid(1,1) = getEnthalpy(tempFluid(t,1));
% 			disp(intEnergyFluid(1,1))
            enthalpyFluid(1,1) = intEnergyFluid(1,1);
		end
	else
        if useTrahanTempProfile
            try             % Calculate temperatures or throw error
                tempFluid(t+1,numRows) = dischargeTemp;
            catch MEx
                disp(t)
                disp(MEx.identifier)
                disp(MEx.message)
                return
            end
        elseif strcmp(conditionSelection,'TrahanPCM')
            tempFluid(t+1,numRows) = tempProfileTrahanPCMDischarge(t,timeStep,chargeTime);
        else
            tempFluid(t+1,numRows) = dischargeTemp;
        end
        intEnergyFluid(1,numRows) = getEnthalpy(tempFluid(t,numRows));
    end
    
    % Fluid advection
    fluidAdvection = timeStep*mdot(t).*(enthalpyCoefficient*enthalpyFluid(1,:)')';

    % Fluid convection
        Re = G*rockAverageDiameter./viscosityFluid;
        Pr = specificHeatFluid.*viscosityFluid./kFluid;
        % Convection - rock section
        if includeRocks
            hpRocks = hvRocks*rockAverageDiameter/6./(voidFractionInverse);
            BiotNumber(rowStartRock:numRows) = hpRocks(rowStartRock:numRows) * rockAverageRadius./kStorage(rowStartRock:numRows);
            hvEffRocks(rowStartRock:numRows) = hvRocks./(1+0.25*BiotNumber(rowStartRock:numRows));
            hvEff(rowStartRock:numRows) = hvEffRocks(rowStartRock:numRows)*hvModifier;
        end
        % Convection - PCM section
        if includePCM
            PrS = Pr; %getSpecificHeat(tempStorage(t,1:rowCutoffPCM)).*getViscosity(tempStorage(t,1:rowCutoffPCM))./getConductivity(tempStorage(t,1:rowCutoffPCM));
            [hpPCM,hvPCM] = calculatehvPCM(G,dEncap,viscosityFluid(1:rowCutoffPCM),Pr(1:rowCutoffPCM),PrS(1:rowCutoffPCM),kFluid(1:rowCutoffPCM),kStorage(1:rowCutoffPCM),surfaceAreaPCM,voidFractionInverse(1:rowCutoffPCM));
            BiotNumber(1:rowCutoffPCM) = hpPCM * rEncap ./ kStorage(1:rowCutoffPCM);
            hvEffPCM = hvPCM./(1+0.25*BiotNumber(1:rowCutoffPCM)); % 0.25 in zang, 0.2 in Trahan
            hvEff(1:rowCutoffPCM) = hvEffPCM(1:rowCutoffPCM)*hvModifier;
		end
        hvEff(1) = 0; hvEff(numRows) = 0;
		convectionExchange = timeStep.*hvEff.*rowVolume.*(tempFluid(t,:)-tempStorage(t,:));

    % Wall loss - spherical particle packed bed section
    if includeBulkHeatLoss
        if TrahanPCM
            alphaConv = kFluid./rockAverageDiameter.*((0.203*(Re.*Pr).^(0.333333)) + (0.220*Re.^0.8.*Pr.^0.4));
        else
            alphaConv = kFluid./rockAverageDiameter.*(3.22*(Re.*Pr).^(0.333333)+0.117*Re.^0.8.*Pr.^0.4);
        end
        uInside = alphaConv + heffWall;
        uWall = 1./(1./uInside + lossInsulationTerms);
        fluidHeatLoss = timeStep.*uWall.*rowWallSurfaceAreaRock.*(ambientTemp-tempFluid(t,:));
		fluidHeatLoss(1) = 0; fluidHeatLoss(numRows) = 0;
        if nonsphericalPCM && includePCM
            fluidHeatLoss(1:rowCutoffPCM) = 0;
        end
    end
    % Wall loss - PCM section
    if includePCMWallExchange && includePCM
        QWallPCM(1:rowCutoffPCM) = timeStep.*kStorage(1:rowCutoffPCM)*storageWallHXareaTerms.*(wallTempIncrease(1:rowCutoffPCM));
	end
    
%%  TEST BLOCK
% 	TEMPSTORAGE = 50* [1:6];
% 	TEMPFLUID = TEMPSTORAGE + 5;
% 	NUMROWS = length(TEMPSTORAGE);
% 	for i=1:1:NUMROWS
% 		if TEMPSTORAGE(i) < TmeltStart
% 			KSTORAGE(i) = 17.7;
% 		elseif TEMPSTORAGE(i) < TmeltEnd
% 			KSTORAGE(i) = 22;
% 		else
% 			KSTORAGE(i) = 23.1;
% 		end
% 	end
% 	KFLUID = getConductivity(TEMPFLUID);
% 	VOIDFRACTION(1:NUMROWS) = 0.40;
% 	ROWLENGTH(1:NUMROWS) = rowLength(20);
% 	
% 	QWALLPCM(1:rowCutoffPCM) = timeStep.*KSTORAGE(1:rowCutoffPCM)*storageWallHXareaTerms.*(wallTempIncrease(1:rowCutoffPCM))
		
% 	[KEFFROCKS,HEFFWALL] = calculatekeff(KSTORAGE,KFLUID,TEMPFLUID,TEMPSTORAGE,VOIDFRACTION,rockAverageDiameter,NUMROWS,emissivityRocks);
% 	CONDUCTIONCOEFFICIENT = zeros(NUMROWS);
% 	for i = 2:1:NUMROWS-1
% 		CONDUCTIONCOEFFICIENT(i,i) = -(KEFFROCKS(i)+KEFFROCKS(i-1));
% 		CONDUCTIONCOEFFICIENT(i,i-1) = KEFFROCKS(i-1);
% 		CONDUCTIONCOEFFICIENT(i,i+1) = KEFFROCKS(i);
% 	end
% 	STORAGECONDUCTION = timeStep*tankArea./ROWLENGTH.*(CONDUCTIONCOEFFICIENT*TEMPSTORAGE')'
				
%% REAL CODE
	
    % Storage radiation
    % more complicated, add later
%     disp('finish radiation terms')
%     return
    if includeRadiation && includePCM && includeRocks
        % between rows 2,3 (cap calculated separately)
        storageRadiation(2) = storageRadiationPCMCoefficient*((tempStorage(t,3)+273.15)^4-(tempStorage(t,2)+273.15)^4);
        
        % rows i-1, i, i+1 from 3 to rowCutoffPCM-1
        storageRadiation(3:rowCutoffPCM-1) = storageRadiationPCMCoefficient* ...
            ((tempStorage(t,2:rowCutoffPCM-2)+273.15).^4 + (tempStorage(t,4:rowCutoffPCM)+273.15).^4 - 2*(tempStorage(t,3:rowCutoffPCM-1)+273.15).^4);
        
        % row rowCutoffPCM-1 to rowCutoffPCM and rowCutoffPCM to
        % rowStartRocks
        storageRadiation(rowCutoffPCM) = storageRadiationPCMCoefficient * ((tempStorage(t,rowCutoffPCM-1)+273.15).^4-(tempStorage(t,rowCutoffPCM)+273.15).^4) + ...
            storageRadiationInterfaceCoefficient * ((tempStorage(t,rowCutoffPCM+1)+273.15).^4-(tempStorage(t,rowCutoffPCM)+273.15).^4);
        
        % row rowStartRocks to rowCutoffPCM
        storageRadiation(rowStartRock) =  storageRadiationInterfaceCoefficient * ((tempStorage(t,rowCutoffPCM)+273.15).^4-(tempStorage(t,rowCutoffPCM+1)+273.15).^4);
    end
    
    % Cap radiation - only relevant in Zanganeh case
    if includeCapRadiation
        storagePlate(2) = storageRadiationPlateCoefficient*((tempStorage(t,2)+plateTempIncrease+273.15)^4-(tempStorage(t,2)+273.15)^4);
    end
    
    if includeConduction
    % Storage conduction    % need to add caveats to conduction between
    % sections.
		if chargeMode
			if zanganeh
				if includePCM
					for i = 2:1:rowCutoffPCM-1  % first row of PCM through second to last row of PCM
						conductionCoefficient(i,i) = -(kStorage(i)+kStorage(i-1));
						conductionCoefficient(i,i-1) = kStorage(i-1);
						conductionCoefficient(i,i+1) = kStorage(i);
					end
					for i = rowStartRock+1:1:numRows-1  % Second row of rock to last row of rock
						conductionCoefficient(i,i) = -(keffRocks(i)+keffRocks(i-1));
						conductionCoefficient(i,i-1) = keffRocks(i-1);
						conductionCoefficient(i,i+1) = keffRocks(i);
					end
					% last row of PCM
						conductionCoefficient(rowCutoffPCM,rowCutoffPCM) = -(kStorage(rowCutoffPCM-1));
						conductionCoefficient(rowCutoffPCM,rowCutoffPCM-1) = kStorage(rowCutoffPCM-1);
						conductionCoefficient(rowCutoffPCM,rowCutoffPCM+1) = 0;
					% first row of rock
						conductionCoefficient(rowStartRock,rowStartRock) = -(keffRocks(rowStartRock));
						conductionCoefficient(rowStartRock,rowStartRock-1) = 0;
						conductionCoefficient(rowStartRock,rowStartRock+1) = keffRocks(rowStartRock);
					% pcm section
					storageConduction(1:rowCutoffPCM) = timeStep*tankArea*fcontPCM./rowLength(1:rowCutoffPCM).* ...
						(conductionCoefficient(1:rowCutoffPCM,1:rowCutoffPCM)*tempStorage(t,1:rowCutoffPCM)')';
					% rock section
					storageConduction(rowStartRock:numRows) = timeStep*tankArea./rowLength(rowStartRock:numRows).* ...
						(conductionCoefficient(rowStartRock:numRows,rowStartRock:numRows)*tempStorage(t,(rowStartRock:numRows))')';
				else
					for i = 2:1:numRows-1
						conductionCoefficient(i,i) = -(keffRocks(i)+keffRocks(i-1));
						conductionCoefficient(i,i-1) = keffRocks(i-1);
						conductionCoefficient(i,i+1) = keffRocks(i);
					end
					storageConduction = timeStep*tankArea./rowLength.*(conductionCoefficient*tempStorage(t,:)')';
				end
			else
				for i = 2:1:numRows-1
					conductionCoefficient(i,i) = -(keffRocks(i)+keffRocks(i-1));
					conductionCoefficient(i,i-1) = keffRocks(i-1);
					conductionCoefficient(i,i+1) = keffRocks(i);
				end
				storageConduction = timeStep*tankArea./rowLength.*(conductionCoefficient*tempStorage(t,:)')';
			end
		else
% ************************************* Need to fill this in!!!!! *************************************
% ************************************* Need to fill this in!!!!! *************************************
% ************************************* Need to fill this in!!!!! *************************************
% ************************************* Need to fill this in!!!!! *************************************
% ************************************* Need to fill this in!!!!! *************************************
% ************************************* Need to fill this in!!!!! *************************************
			if zanganeh % Need to fill this in!!!!!
				if includePCM
					for i = 2:1:rowCutoffPCM-1  % first row of PCM through second to last row of PCM
						conductionCoefficient(i,i) = -(kStorage(i)+kStorage(i-1));
						conductionCoefficient(i,i-1) = kStorage(i);
						conductionCoefficient(i,i+1) = kStorage(i);
					end
					for i = rowStartRock+1:1:numRows-1  % Second row of rock to last row of rock
						conductionCoefficient(i,i) = -(keffRocks(i)+keffRocks(i-1));
						conductionCoefficient(i,i-1) = keffRocks(i-1);
						conductionCoefficient(i,i+1) = keffRocks(i);
					end
					% last row of PCM
						conductionCoefficient(rowCutoffPCM,rowCutoffPCM) = -(kStorage(rowCutoffPCM-1));
						conductionCoefficient(rowCutoffPCM,rowCutoffPCM-1) = kStorage(rowCutoffPCM-1);
						conductionCoefficient(rowCutoffPCM,rowCutoffPCM+1) = 0;
					% first row of rock
						conductionCoefficient(rowStartRock,rowStartRock) = -(keffRocks(rowStartRock));
						conductionCoefficient(rowStartRock,rowStartRock-1) = 0;
						conductionCoefficient(rowStartRock,rowStartRock+1) = keffRocks(rowStartRock);
					% pcm section
					storageConduction(1:rowCutoffPCM) = timeStep*tankArea*fcontPCM./rowLength(1:rowCutoffPCM).* ...
						(conductionCoefficient(1:rowCutoffPCM,1:rowCutoffPCM)*tempStorage(t,1:rowCutoffPCM)')';
					% rock section
					storageConduction(rowStartRock:numRows) = timeStep*tankArea./rowLength(rowStartRock:numRows).* ...
						(conductionCoefficient(rowStartRock:numRows,rowStartRock:numRows)*tempStorage(t,(rowStartRock:numRows))')';
				else
					for i = 2:1:numRows-1
						conductionCoefficient(i,i) = -(keffRocks(i)+keffRocks(i-1));
						conductionCoefficient(i,i-1) = keffRocks(i-1);
						conductionCoefficient(i,i+1) = keffRocks(i);
					end
					storageConduction = timeStep*tankArea./rowLength.*(conductionCoefficient*tempStorage(t,:)')';
				end
			else
				for i = 2:1:numRows-1
					conductionCoefficient(i,i) = -(keffRocks(i)+keffRocks(i-1));
					conductionCoefficient(i,i-1) = keffRocks(i-1);
					conductionCoefficient(i,i+1) = keffRocks(i);
				end
				storageConduction = timeStep*tankArea./rowLength.*(conductionCoefficient*tempStorage(t,:)')';
			end
		end
    end
    
    
    
%     disp(intEnergyFluid(1,1))
% 	convectionExchange(rowStartRock+5:numRows) = 0;
    intEnergyStorage(2,:) = intEnergyStorage(1,:) + (convectionExchange + conductionModifier*storageConduction  + storageRadiation + storagePlate + QWallPCM)./rowMassStorage;
    intEnergyFluid(2,:) = intEnergyFluid(1,:) + (fluidAdvection - convectionExchange + fluidHeatLoss)./rowMassFluid;

		
% 	Abba = (fluidAdvection - convectionExchange + fluidHeatLoss)./rowMassFluid(i);
% 	disp(intEnergyFluid(1,1))
	
    
    if t == logTime
%         QWallPCM = zeros(1,numRows);
%         storagePlate = zeros(1,numRows);
        logEnergyFlow = [tempStorage(t,:);	tempFluid(t,:);			intEnergyStorage(1,:);
			intEnergyStorage(2,:);			intEnergyFluid(1,:);	intEnergyFluid(2,:); 
            convectionExchange;				storageConduction;		fluidAdvection; 
            QWallPCM;						storagePlate;			storageRadiation;
			fluidHeatLoss;					hvEff;					keffRocks;
			BiotNumber]';      
    end
    logConvectionStorage(t) = sum(-convectionExchange.*voidFractionInverse(i).*densityStorage.*rowVolume(i));
    logQRadPlate(t) = storagePlate(2);
    logQWallPCM(t) = sum(QWallPCM);
    logRockLoss(t) = sum(fluidHeatLoss);
        
% all of this needs to get processed into matrix format
    
    if chargeMode
        logEnergyInput(t) = timeStep*mdot(t)*(enthalpyFluid(1)-enthalpyFluid(numRows));
    else
        logEnergyExtracted(t) = timeStep*mdot(t)*(enthalpyFluid(1)-enthalpyFluid(numRows));
    end
    
    % 
    if chargeMode
        if useInternalEnergy
            intEnergyFluid(2,1) = getInternalEnergy(tempFluid(t+1,1));
            intEnergyFluid(2,numRows-1:numRows) = intEnergyFluid(2,numRows-2);
		else
			% this is where the unexpected change is happening. Why?
% 			fprintf('Fluid temp used for calculateTempFluid inlet enthalpy: %g\n',tempFluid(t+1,1))
            intEnergyFluid(2,1) = getEnthalpy(tempFluid(t+1,1));
            intEnergyFluid(2,numRows) = intEnergyFluid(2,numRows-1);
% 			if t == 2
% 				return
% 			end
        end
    else
        if useInternalEnergy
            intEnergyFluid(2,numRows) = getInternalEnergy(tempFluid(t+1,numRows));
        else
            intEnergyFluid(2,numRows) = getEnthalpy(tempFluid(t+1,numRows));
        end
	end
	
    intEnergyStorage(2,1) = intEnergyStorage(2,2);
    intEnergyStorage(2,numRows) = intEnergyStorage(2,numRows-1);
    
%     if includePCM
%     logEnergyStored(t) = (sum(intEnergyStorage(2:rowCutoffPCM,2)) - sum(intEnergyStorage(2:rowCutoffPCM,1)))*densityPCMEff*rowVolume(i)*voidFractionInverse(i) ...
%         + (sum(intEnergyStorage(rowStartRock:numRows-1,2)) - sum(intEnergyStorage(rowStartRock:numRows-1,1)))*densityRocks*rowVolume(i)*voidFractionInverse(i);
%     else
%         logEnergyStored(t) = (sum(intEnergyStorage(2:numRows-1,2)) - sum(intEnergyStorage(2:numRows-1,1)))*densityRocks*rowVolume(i)*voidFractionInverse(i);
%     end
       
    try             % Calculate temperatures or throw error
        for i = 1:1:rowCutoffPCM
            if tempStorage(t,i) > TmeltEnd
                cPCM(1,i) = ceffPCMLiquid;
            elseif tempStorage(t,i) > TmeltStart
                cPCM(1,i) = ceffPCMPhaseChange;
            else
                cPCM(1,i) = ceffPCMSolid;
            end
        end
%         if tempStorage(t,4) > TmeltEnd && TmeltEnd > 500
%             disp(cPCM)
%             return
%         end
        
        % storage temps don't depend on int energy vs enthalpy
        if includePCM
            tempStorage(t+1,1:rowCutoffPCM) = calculateTempPCM(tempStorage(t,1:rowCutoffPCM),intEnergyStorage(2,1:rowCutoffPCM),intEnergyStorage(1,1:rowCutoffPCM),cPCM(1,1:rowCutoffPCM));
        end
        tempStorage(t+1,rowStartRock:numRows) = calculateTempRock(tempStorage(t,rowStartRock:numRows),intEnergyStorage(2,rowStartRock:numRows),intEnergyStorage(1,rowStartRock:numRows),getRockSpecificHeat);
%         disp(intEnergyFluid(2,1:5)-intEnergyFluid(1,1:5))
        tempFluid(t+1,:) = calculateTempFluid(intEnergyFluid(2,:),intEnergyFluid(1,:), tempFluid(t,:), specificHeatFluid,getEnthalpy,getSpecificHeat);
%         tempfluid = calculateTempFluid(intEnergyFluid(2,:),intEnergyFluid(1,:), tempFluid(t,:), specificHeatFluid,getEnthalpy,getSpecificHeat)
%         disp(tempFluid(t+1,1:5))
    catch MEx
        disp(t)
        disp(MEx.identifier)
        disp(MEx.message)
        fprintf(['Ended at time ' datestr(clock,'yyyy/mm/dd--HH:MM:SS') '\n'])
        return
    end
    
    logFluidTempChangeActual(t) = tempFluid(t+1,10) - tempFluid(t,10);
    if calculatePressureDrop
        if includePCM   % Pressure drop calculations
            pressureDrop(t) = pressureDrop(t) + getPressureDrop(rowLengthPCM,dEncap, voidFraction(2:rowCutoffPCM), voidFractionInverse(2:rowCutoffPCM), tempFluid(t+1,2:rowCutoffPCM)-tempFluid(t,2:rowCutoffPCM), tempFluid(t,2:rowCutoffPCM),viscosityFluid(2:rowCutoffPCM),densityFluid(2:rowCutoffPCM),G);
        end
        if includeRock
            pressureDrop(t) = pressureDrop(t) + getPressureDrop(rowLengthRock,rockAverageDiameter, voidFraction(rowStartRock:numRows-1), voidFractionInverse(rowStartRock:numRows-1), tempFluid(t+1,rowStartRock:numRows-1)-tempFluid(t,rowStartRock:numRows-1), tempFluid(t,rowStartRock:numRows-1),viscosityFluid(rowStartRock:numRows-1),densityFluid(rowStartRock:numRows-1),G);
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
            break
        end
    end
    t=t+1;
end

disp('Transient analysis complete')
toc
[includeConduction,includeRadiation,includeBulkHeatLoss,includePCMWallExchange,includeCapRadiation]
if breakAfterLogTime
	return
end
%% Energy logger calculations
% {
if printEnergyFlows
    if breakAfterLogTime == false
        
        energyInput = sum(logEnergyInput);
        energyExtracted = sum(logEnergyExtracted);
        
        energyInputKWH = energyInput/3.6e6;
        energyExtractedKWH = energyExtracted/3.6e6;
        if strcmp(conditionSelection,'Zanganeh2015') || strcmp(conditionSelection,'Trahan')
            if totalTime > chargeTime+1
                if includePCM
                    totalEnergyRocksEndDischarge = sum(intEnergyStorage(2,rowStartRock:numRows-1));
                    totalEnergyPCMEndDischarge =  sum(intEnergyStorage(2,2:rowCutoffPCM));
                else
                    totalEnergyRocksEndDischarge = sum(intEnergyStorage(2,2:numRows-1));
                end
                energyStoredRocks = (totalEnergyRocksEndCharge - totalEnergyRocksStart)*densityRocks*rowVolumeRock*voidFractionRocksInverse;
                energyRemovedRocks = (totalEnergyRocksEndCharge - totalEnergyRocksEndDischarge)*densityRocks*rowVolumeRock*voidFractionRocksInverse;
            else % prevents code from breaking if breakAfterLogTime=true
                totalEnergyRocksEndCharge = sum(intEnergyStorage(1,2:numRows-1));
                energyStoredRocks = (totalEnergyRocksEndCharge - totalEnergyRocksStart)*densityRocks*rowVolumeRock*voidFractionRocksInverse;
                energyRemovedRocks = 0;
            end
            energyStoredRocksKWH = energyStoredRocks/3.6e6;
            energyRemovedRocksKWH = energyRemovedRocks/3.6e6;
            energyWallLossRocksChargingKWH = sum(logRockLoss(1:chargeTime))/3.6e6;
            energyWallLossRocksDischargingKWH = sum(logRockLoss(chargeTime+1:totalTime))/3.6e6;
        else
            energyStoredRocksKWH = 0;
            energyWallLossRocksChargingKWH =sum(logRockLoss(1:chargeTime))/3.6e6;
            energyWallLossRocksDischargingKWH=0;
        end

        if includePCM
    %         exist totalEnergyPCMEndCharge
            if ~exist('totalEnergyPCMEndCharge')
                if strcmp(conditionSelection,'TrahanPCM')
                    totalEnergyPCMEndCharge = sum(intEnergyStorage(1,2:numRows-1));
                else
                    totalEnergyPCMEndCharge = sum(intEnergyStorage(1,2:rowCutoffPCM));
                end
                totalEnergyPCMEndDischarge =  0; %sum(intEnergyStorage(2:rowCutoffPCM,1));
            end

            energyStoredPCM = (totalEnergyPCMEndCharge - totalEnergyPCMStart)*densityPCMEff*rowVolumePCM*voidFractionPCMInverse;
            energyRemovedPCM = (totalEnergyPCMEndCharge - totalEnergyPCMEndDischarge)*densityPCMEff*rowVolumePCM*voidFractionPCMInverse;
            energyStoredPCMKWH = energyStoredPCM/3.6e6;
            energyRemovedPCMKWH = energyRemovedPCM / 3.6e6;
            energyWallLossPCMChargingKWH = sum(logPCMLoss(1:chargeTime))/3.6e6;
            energyWallLossPCMDischargingKWH = sum(logPCMLoss(chargeTime+1:totalTime))/3.6e6;
            energyWallPCMKWHCharge = sum(logQWallPCM(1:chargeTime))/3.6e6;
            energyWallPCMKWHDischarge = sum(logQWallPCM(chargeTime+1:totalTime))/3.6e6;
        else
            energyStoredPCMKWH = 0;
            energyWallLossPCMChargingKWH = 0;
            energyWallLossPCMDischargingKWH = 0;
            energyWallPCMKWHCharge = 0;
            energyWallPCMKWHDischarge = 0;
            energyRemovedPCMKWH = 0;
        end
        if beasley || dischargeTimeHours == 0
            energyRemovedRocksKWH = 0;
        end

        energyRadPlateKWHCharge = sum(logQRadPlate(1:chargeTime))/3.6e6;
        energyRadPlateKWHDischarge = sum(logQRadPlate(chargeTime+1:totalTime))/3.6e6;

        energyBalanceChargingKWH = energyInputKWH-energyStoredRocksKWH-energyStoredPCMKWH+ ...
            energyRadPlateKWHCharge+energyWallLossRocksChargingKWH+energyWallPCMKWHCharge;
    %     if dischargeTimeHours > 0
            energyBalanceDischargingKWH = -energyExtractedKWH+energyRemovedRocksKWH+energyRemovedPCMKWH+ ...
                energyRadPlateKWHDischarge+energyWallLossRocksDischargingKWH+energyWallPCMKWHDischarge;
    %     else
    %         energyBalanceDischargingKWH = 0;
	end
end
%% Error calculations
if compareToValidation && trackError
    % model
    errorLinfFluidModel = max(abs(fluiderrorModel(:,2))); errorLinfR1Model = max(abs(row1errorModel(:,2)));
    errorLinfR2Model = max(abs(row2errorModel(:,2))); errorLinfR3Model = max(abs(row3errorModel(:,2))); errorLinfR4Model = max(abs(row4errorModel(:,2)));
    
    errorL1FluidModel = sum(abs(fluiderrorModel(:,2))); errorL1R1Model = sum(abs(row1errorModel(:,2)));
    errorL1R2Model = sum(abs(row2errorModel(:,2))); errorL1R3Model = sum(abs(row3errorModel(:,2))); errorL1R4Model = sum(abs(row4errorModel(:,2)));
    
    errorL2FluidModel = sum(fluiderrorModel(:,2).^2)^0.5; errorL2R1Model = sum(row1errorModel(:,2).^2)^0.5;
    errorL2R2Model = sum(row2errorModel(:,2).^2)^0.5; errorL2R3Model = sum(row3errorModel(:,2).^2)^0.5; errorL2R4Model = sum(row4errorModel(:,2).^2)^0.5;
    
    rmsfluidModel = 0; rmsR1Model = 0; rmsR2Model = 0; rmsR3Model = 0; rmsR4Model = 0;
    for i = 1:1:length(fluiderrorModel)
        if abs(fluiderrorModel(i,1))>0 rmsfluidModel = rmsfluidModel+fluiderrorModel(i,2)^2; nfluid=i; end
        if abs(row1errorModel(i,1))>0 rmsR1Model = rmsR1Model+row1errorModel(i,2)^2; nR1=i; end
        if abs(row2errorModel(i,1))>0 rmsR2Model = rmsR2Model+row2errorModel(i,2)^2; nR2=i; end
        if abs(row3errorModel(i,1))>0 rmsR3Model = rmsR3Model+row3errorModel(i,2)^2; nR3=i; end
        if abs(row4errorModel(i,1))>0 rmsR4Model = rmsR4Model+row4errorModel(i,2)^2; nR4=i; end
    end
    rmsfluidModel = sqrt(rmsfluidModel/nfluid); rmsR1Model = sqrt(rmsR1Model/nR1); rmsR2Model = sqrt(rmsR2Model/nR2); 
    rmsR3Model = sqrt(rmsR3Model/nR3); rmsR4Model = sqrt(rmsR4Model/nR4); 
    
    % experimental
    errorLinfFluidExp = max(abs(fluiderrorExp(:,2))); errorLinfR1Exp = max(abs(row1errorExp(:,2)));
    errorLinfR2Exp = max(abs(row2errorExp(:,2))); errorLinfR3Exp = max(abs(row3errorExp(:,2))); errorLinfR4Exp = max(abs(row4errorExp(:,2)));
    
    errorL1FluidExp = sum(abs(fluiderrorExp(:,2))); errorL1R1Exp = sum(abs(row1errorExp(:,2)));
    errorL1R2Exp = sum(abs(row2errorExp(:,2))); errorL1R3Exp = sum(abs(row3errorExp(:,2))); errorL1R4Exp = sum(abs(row4errorExp(:,2)));
    
    errorL2FluidExp = sum(fluiderrorExp(:,2).^2)^0.5; errorL2R1Exp = sum(row1errorExp(:,2).^2)^0.5;
    errorL2R2Exp = sum(row2errorExp(:,2).^2)^0.5; errorL2R3Exp = sum(row3errorExp(:,2).^2)^0.5; errorL2R4Exp = sum(row4errorExp(:,2).^2)^0.5;
    
    rmsfluidExp = 0; rmsR1Exp = 0; rmsR2Exp = 0; rmsR3Exp = 0; rmsR4Exp = 0;
    for i = 1:1:length(fluiderrorExp)
        if abs(fluiderrorExp(i,1))>0 rmsfluidExp = rmsfluidExp+fluiderrorExp(i,2)^2; nfluid=i; end
        if abs(row1errorExp(i,1))>0 rmsR1Exp = rmsR1Exp+row1errorExp(i,2)^2; nR1=i; end
        if abs(row2errorExp(i,1))>0 rmsR2Exp = rmsR2Exp+row2errorExp(i,2)^2; nR2=i; end
        if abs(row3errorExp(i,1))>0 rmsR3Exp = rmsR3Exp+row3errorExp(i,2)^2; nR3=i; end
        if abs(row4errorExp(i,1))>0 rmsR4Exp = rmsR4Exp+row4errorExp(i,2)^2; nR4=i; end
    end 
    rmsfluidExp = sqrt(rmsfluidExp/nfluid); rmsR1Exp = sqrt(rmsR1Exp/nR1); rmsR2Exp = sqrt(rmsR2Exp/nR2); 
    rmsR3Exp = sqrt(rmsR3Exp/nR3); rmsR4Exp = sqrt(rmsR4Exp/nR4); 
    
    %  PCM1errorModel fluid1errorModel
    
    if includePCM
        % model
        errorLinfPCM1Model = max(abs(PCM1errorModel(:,2))); errorLinfFluid1Model = max(abs(fluid1errorModel(:,2)));
        errorLinfPCM3Model = max(abs(PCM3errorModel(:,2))); errorLinfFluid3Model = max(abs(fluid3errorModel(:,2))); 

        errorL1PCM1Model = sum(abs(PCM1errorModel(:,2))); errorL1Fluid1Model = sum(abs(fluid1errorModel(:,2)));
        errorL1PCM3Model = sum(abs(PCM3errorModel(:,2))); errorL1Fluid3Model = sum(abs(fluid3errorModel(:,2)));

        errorL2PCM1Model = sum(PCM1errorModel(:,2).^2)^0.5; errorL2Fluid1Model = sum(fluid1errorModel(:,2).^2)^0.5;
        errorL2PCM3Model = sum(PCM3errorModel(:,2).^2)^0.5; errorL2Fluid3Model = sum(fluid3errorModel(:,2).^2)^0.5;
        
        rmsPCM1Model = 0; rmsFluid1Model = 0; rmsPCM3Model = 0; rmsFluid3Model = 0;
        for i = 1:1:length(PCM1errorModel)
            if abs(PCM1errorModel(i,1))>0 rmsPCM1Model = rmsPCM1Model+PCM1errorModel(i,2)^2; nR1=i; end
            if abs(fluid1errorModel(i,1))>0 rmsFluid1Model = rmsFluid1Model+fluid1errorModel(i,2)^2; nR2=i; end
            if abs(PCM3errorModel(i,1))>0 rmsPCM3Model = rmsPCM3Model+PCM3errorModel(i,2)^2; nR3=i; end
            if abs(fluid3errorModel(i,1))>0 rmsFluid3Model = rmsFluid3Model+fluid3errorModel(i,2)^2; nR4=i; end
        end
        rmsPCM1Model = sqrt(rmsPCM1Model/nR1); rmsFluid1Model = sqrt(rmsFluid1Model/nR2); 
        rmsPCM3Model = sqrt(rmsPCM3Model/nR3); rmsFluid3Model = sqrt(rmsFluid3Model/nR4); 
    

        % experimental
        errorLinfPCM1Exp = max(abs(PCM1errorExp(:,2))); errorLinfFluid1Exp = max(abs(fluid1errorExp(:,2)));
        errorLinfPCM3Exp = max(abs(PCM3errorExp(:,2))); errorLinfFluid3Exp = max(abs(fluid3errorExp(:,2))); 

        errorL1PCM1Exp = sum(abs(PCM1errorExp(:,2))); errorL1Fluid1Exp = sum(abs(fluid1errorExp(:,2)));
        errorL1PCM3Exp = sum(abs(PCM3errorExp(:,2))); errorL1Fluid3Exp = sum(abs(fluid3errorExp(:,2)));

        errorL2PCM1Exp = sum(PCM1errorExp(:,2).^2)^0.5; errorL2Fluid1Exp = sum(fluid1errorExp(:,2).^2)^0.5;
        errorL2PCM3Exp = sum(PCM3errorExp(:,2).^2)^0.5; errorL2Fluid3Exp = sum(fluid3errorExp(:,2).^2)^0.5;
        
        rmsPCM1Exp = 0; rmsFluid1Exp = 0; rmsPCM3Exp = 0; rmsFluid3Exp = 0;
        for i = 1:1:length(PCM1errorExp)
            if abs(PCM1errorExp(i,1))>0 rmsPCM1Exp = rmsPCM1Exp+PCM1errorExp(i,2)^2; nR1=i; end
            if abs(fluid1errorExp(i,1))>0 rmsFluid1Exp = rmsFluid1Exp+fluid1errorExp(i,2)^2; nR2=i; end
            if abs(PCM3errorExp(i,1))>0 rmsPCM3Exp = rmsPCM3Exp+PCM3errorExp(i,2)^2; nR3=i; end
            if abs(fluid3errorExp(i,1))>0 rmsFluid3Exp = rmsFluid3Exp+fluid3errorExp(i,2)^2; nR4=i; end
        end
        rmsPCM1Exp = sqrt(rmsPCM1Exp/nR1); rmsFluid1Exp = sqrt(rmsFluid1Exp/nR2); 
        rmsPCM3Exp = sqrt(rmsPCM3Exp/nR3); rmsFluid3Exp = sqrt(rmsFluid3Exp/nR4); 

    end
end


%% Written Output


    if includePCM
        fprintf('Total energy input: %g kWh -- target 18.4\n',energyInputKWH)
        fprintf('Total energy stored in rocks: %g kWh\n',energyStoredRocksKWH)
        fprintf('Total energy stored in PCM: %g kWh\n',energyStoredPCMKWH)
        fprintf('Energy exchange with top plate during charge: %g kWh -- target 0.95\n',energyRadPlateKWHCharge);
        fprintf('Energy exchange with top plate during discharge: %g kWh -- target -0.14\n',energyRadPlateKWHDischarge);
        fprintf('Energy lost from rocks section during charging: %g kWh -- target -0.24\n',energyWallLossRocksChargingKWH)
        fprintf('Energy lost from rocks section during discharging: %g kWh -- target -0.25\n',energyWallLossRocksDischargingKWH)
        fprintf('Energy exchange between PCM section and wall during charging: %g kWh -- target 1.27\n',energyWallPCMKWHCharge)
        fprintf('Energy exchange between PCM section and wall during discharging: %g kWh -- target 0.43\n',energyWallPCMKWHDischarge)
        fprintf('Overall energy balance during charging: %g kWh = %g%% of input -- target 0\n',energyBalanceChargingKWH,energyBalanceChargingKWH/(energyInputKWH+energyRadPlateKWHCharge+energyWallPCMKWHCharge)*100)
        fprintf('Overall energy balance during discharging: %g kWh = %g%% of input -- target 0\n',energyBalanceDischargingKWH,energyBalanceDischargingKWH/(energyRemovedRocksKWH+energyRemovedPCMKWH)*100)
        fprintf('Total energy output: %g kWh\n',energyExtractedKWH)
        fprintf('Total energy removed from rocks: %g kWh\n',energyRemovedRocksKWH)
        fprintf('Total energy removed from PCM: %g kWh\n',energyRemovedPCMKWH)
        
        energyFlowOutputs = [energyInputKWH;energyStoredRocksKWH;energyStoredPCMKWH;energyRadPlateKWHCharge;energyRadPlateKWHDischarge;
            energyWallLossRocksChargingKWH;energyWallLossRocksDischargingKWH;energyWallPCMKWHCharge;energyWallPCMKWHDischarge;energyBalanceChargingKWH];
    else
        fprintf('Total energy input: %g kWh -- target 19.3\n',energyInputKWH)
        fprintf('Total energy stored in rocks: %g kWh\n',energyStoredRocksKWH)
        fprintf('Energy exchange with top plate during charge: %g kWh -- target 1.01\n',energyRadPlateKWHCharge);
        fprintf('Energy exchange with top plate during discharge: %g kWh -- target +0.42\n',energyRadPlateKWHDischarge);
        fprintf('Energy lost from rocks section during charging: %g kWh -- target -0.3\n',energyWallLossRocksChargingKWH)
        fprintf('Energy lost from rocks section during discharging: %g kWh -- target -0.3\n',energyWallLossRocksDischargingKWH)
        fprintf('Overall energy balance during charging: %g kWh = %g%% of input -- target 0\n',energyBalanceChargingKWH,energyBalanceChargingKWH/energyInputKWH*100)
        fprintf('Overall energy balance during discharging: %g kWh = %g%% of input -- target 0\n',energyBalanceDischargingKWH,energyBalanceDischargingKWH/(energyRemovedRocksKWH)*100)
        fprintf('Total energy output: %g kWh\n',energyExtractedKWH)
        fprintf('Total energy removed from rocks: %g kWh\n',energyRemovedRocksKWH)
        energyFlowOutputs = [energyInputKWH;energyStoredRocksKWH;0; energyRadPlateKWHCharge;energyRadPlateKWHDischarge;
            energyWallLossRocksChargingKWH;energyWallLossRocksDischargingKWH;0;0;energyBalanceChargingKWH];
    end
    
    if compareToValidation && trackError && printLNorms
        if beasley
            fprintf('\n----------------Errors vs model-----------------\n')
            fprintf('\t\tLinf\t\tL1\t\t\tL2\t\tRMS\n')
            fprintf('R1:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR1Model,errorL1R1Model,errorL2R1Model)
            fprintf('R2:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR2Model,errorL1R2Model,errorL2R2Model)
            fprintf('\n----------------Errors vs exp-----------------\n')
            fprintf('\t\tLinf\t\tL1\t\t\tL2\t\tRMS\n')
            fprintf('R1:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR1Exp,errorL1R1Exp,errorL2R1Exp)
            fprintf('R2:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR2Exp,errorL1R2Exp,errorL2R2Exp)
        else
            fprintf('\n----------------Rock Section--------------------\n')
            fprintf('----------------Errors vs model-----------------\n')
            fprintf('\t\tLinf\t\tL1\t\t\tL2\t\tRMS\n')
            fprintf('fluid:\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfFluidModel,errorL1FluidModel,errorL2FluidModel,rmsfluidModel)
            fprintf('R1:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR1Model,errorL1R1Model,errorL2R1Model,rmsR1Model)
            fprintf('R2:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR2Model,errorL1R2Model,errorL2R2Model,rmsR2Model)
            fprintf('R3:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR3Model,errorL1R3Model,errorL2R3Model,rmsR3Model)
            fprintf('R4:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR4Model,errorL1R4Model,errorL2R4Model,rmsR4Model)
            fprintf('----------------Errors vs exp-----------------\n')
            fprintf('\t\tLinf\t\tL1\t\t\tL2\t\tRMS\n')
            fprintf('fluid:\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfFluidExp,errorL1FluidExp,errorL2FluidExp,rmsfluidExp)
            fprintf('R1:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR1Exp,errorL1R1Exp,errorL2R1Exp,rmsR1Exp)
            fprintf('R2:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR2Exp,errorL1R2Exp,errorL2R2Exp,rmsR2Exp)
            fprintf('R3:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR3Exp,errorL1R3Exp,errorL2R3Exp,rmsR3Exp)
            fprintf('R4:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfR4Exp,errorL1R4Exp,errorL2R4Exp,rmsR4Exp)
            if includePCM
                fprintf('\n----------------PCM Section-----------------\n')
                fprintf('----------------Errors vs model-----------------\n')
                fprintf('\t\t\tLinf\t\tL1\t\t\tL2\t\tRMS\n')
                fprintf('PCM1:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfPCM1Model,errorL1PCM1Model,errorL2PCM1Model,rmsPCM1Model)
                fprintf('Fluid1:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfFluid1Model,errorL1Fluid1Model,errorL2Fluid1Model,rmsFluid1Model)
                fprintf('PCM3:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfPCM3Model,errorL1PCM3Model,errorL2PCM3Model,rmsPCM3Model)
                fprintf('Fluid3:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfFluid3Model,errorL1Fluid3Model,errorL2Fluid3Model,rmsFluid3Model)
                fprintf('----------------Errors vs exp-----------------\n')
                fprintf('\t\t\tLinf\t\tL1\t\t\tL2\t\tRMS\n')
                fprintf('PCM1:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfPCM1Exp,errorL1PCM1Exp,errorL2PCM1Exp,rmsPCM1Exp)
                fprintf('Fluid1:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfFluid1Exp,errorL1Fluid1Exp,errorL2Fluid1Exp,rmsFluid1Exp)
                fprintf('PCM3:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfPCM3Exp,errorL1PCM3Exp,errorL2PCM3Exp,rmsPCM3Exp)
                fprintf('Fluid3:\t\t%g\t\t%g\t\t%g\t\t%g\n',errorLinfFluid3Exp,errorL1Fluid3Exp,errorL2Fluid3Exp,rmsFluid3Exp)
            end
        end
    end
% disp([energyInputKWH; (energyStoredPCMKWH+energyWallLossRocksChargingKWH)]')
% end

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
    t=1;
    while t*timeStepInverse < totalTime
        logTempFluid(:,t) = tempFluid(t*timeStepInverse,:);
        logTempStorage(:,t) = tempStorage(t*timeStepInverse,:);
        logPressure(t) = pressureDrop(t*timeStepInverse);
        logmdot(t) = mdot(t*timeStepInverse);
        logvdot(t) = vdot(t*timeStepInverse);
        logAverageTemp(:,t) = (logTempFluid(:,t)+logTempStorage(:,t))/2;
        logFluidEnergyChangeOutput(:,t) = logFluidEnergyChange(t*timeStepInverse);
        logFluidTempChangeTheoreticalOutput(:,t) = logFluidTempChangeTheoretical(t*timeStepInverse);
        logFluidTempChangeActualOutput(:,t) = logFluidTempChangeActual(t*timeStepInverse);

        t=t+1;
    end
    if strcmp('Zanganeh2015',conditionSelection)
        if includePCM
            parameters = {'Condition: ',conditionSelection;'Fluid: ',fluidName; ...
                'Reference fluid: ',referenceFluid; 'PCM Included?', includePCM; ...
                'PCM Height', heightPCM;'IntEnergy?',useInternalEnergy; ...
                'Conduction',includeConduction; ...
                'Radiation',includeRadiation; 'Heat loss',includeBulkHeatLoss; ...
                'Wall loss', includePCMWallExchange; 'Cap radiation', includeCapRadiation; ...
				
				'Energy input kWh', energyInputKWH; 'Energy stored in rocks kWh', energyStoredRocksKWH; ...
				'Energy stored PCM kWh',energyStoredPCMKWH; 'Qradplate charge', energyRadPlateKWHCharge;
				'Qradplate discharge',energyRadPlateKWHDischarge; 'Qwallrock charge', energyWallLossRocksChargingKWH;
				'Qwallrock discharge',energyWallLossRocksDischargingKWH; 'Qwallpcm charge', energyWallPCMKWHCharge;
				'Qwallpcm discharge',energyWallPCMKWHDischarge';'Energy balance charging kWh',energyBalanceChargingKWH;	
				};
    %         parameters = {'Condition: ',conditionSelection,'Fluid: ',fluidName; 'Reference fluid: ',referenceFluid; 'PCM Included?', includePCM; 'PCM Height', heightPCM;'Ein',energyInputKWH;'Erocks',energyStoredRocksKWH;'EPCM',energyStoredPCMKWH;'Eradplatech',energyRadPlateKWHCharge;'Eradplatedisch',energyRadPlateKWHDischarge;'EwallPCMch',energyWallPCMKWHCharge;'EwallPCMdisch',energyWallPCMKWHDischarge;'Ewallrocksch',energyWallLossRocksChargingKWH;'Ewallrocksdisch',energyWallLossRocksDischargingKWH};		
		else
            parameters = {'Condition: ',conditionSelection;'Fluid: ',fluidName; ...
                'Reference fluid: ',referenceFluid; 'PCM Included?', includePCM; ...
                'PCM Height', 0;'IntEnergy?',useInternalEnergy; ...
                'Conduction',includeConduction; ...
                'Radiation',includeRadiation; 'Heat loss',includeBulkHeatLoss; ...
                'Wall loss', includePCMWallExchange; 'Cap radiation', includeCapRadiation
				
				'Energy input kWh', energyInputKWH; 'Energy stored in rocks kWh', energyStoredRocksKWH; ...
				'Qradplate charge', energyRadPlateKWHCharge; 'Qradplate discharge',energyRadPlateKWHDischarge; 
				'Qwallrock charge', energyWallLossRocksChargingKWH; 'Qwallrock discharge',energyWallLossRocksDischargingKWH;
				'Energy balance charging kWh',energyBalanceChargingKWH;	
				};
    %         parameters = {'Condition: ',conditionSelection,'Fluid: ',fluidName; 'Reference fluid: ',referenceFluid; 'PCM Included?', includePCM; 'PCM Height', 'n/a'; 'Ein',energyInputKWH; 'Erocks',energyStoredRocksKWH;'EPCM', 'n/a';'Eradplatech',energyRadPlateKWHCharge;'Eradplatedisch',energyRadPlateKWHDischarge;'EwallPCMch', 'n/a';'EwallPCMdisch','n/a'; 'Ewallrocksch',energyWallLossRocksChargingKWH;'Ewallrocksdisch',energyWallLossRocksDischargingKWH};   
        end
    else
        parameters = {'Condition: ',conditionSelection;'Fluid: ',fluidName; 'Reference fluid: ',referenceFluid; 'PCM Included?', includePCM; 'PCM Height', heightPCM};
	end
    workSpaceName = ['.\Data Files\Excel Dumps\' conditionSelection '--' num2str(includePCM) '--' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.xlsx'];
    % writematrix(A,workSpaceName)
    % type(workSpaceName)
    warning('off','MATLAB:xlswrite:AddSheet');
    try
        disp('Attempting to write temperature histories to file...')
        % parameters
        xlswrite(workSpaceName,parameters,'Parameters','A1:B26');
        disp('Wrote model parameters to file')
        
        % Main tank, both PCM and rock cases
		if exist('R5')
			mainTankOutput = [logTempFluid(1,:); 
				logTempFluid(R1,:); 
				logTempStorage(R1,:);
				logTempFluid(R2,:); 
				logTempStorage(R2,:);
				logTempFluid(R3,:); 
				logTempStorage(R3,:);
				logTempFluid(R4,:); 
				logTempStorage(R4,:);
				logTempFluid(R5,:);
				logTempStorage(R5,:);
				logmdot(:)' ]';
		else
			mainTankOutput = [logTempFluid(1,:); 
				logTempFluid(R1,:); 
				logTempStorage(R1,:);
				logTempFluid(R2,:); 
				logTempStorage(R2,:);
				logTempFluid(R3,:); 
				logTempStorage(R3,:);
				logTempFluid(R4,:); 
				logTempStorage(R4,:);
				logmdot(:)' ]';
		end
			
			
        xlswrite(workSpaceName,mainTankOutput,'Tank Temps',['A2:L' num2str(1+length(logTempFluid))]);
        disp('Wrote tank temperatures to file')
        
        % PCM
        if includePCM && strcmp(conditionSelection,'Zanganeh2015')
            pcmSectionOutput = [ 
            logTempFluid(2,:); 
            logTempFluid(4,:); 
            logTempStorage(2,:);
            logTempStorage(4,:) ]';
            xlswrite(workSpaceName,pcmSectionOutput,'PCM Temps',['A1:D' num2str(1+length(logTempFluid))]);
            disp('Wrote PCM top temperatures to file')
        end
        
%         xlswrite(workSpaceName,logEnergyFlowJ,'logEnergyFlowJ','A1:M72');
        xlswrite(workSpaceName,logEnergyFlow,'logEnergyFlow','A1:M72');
        disp('Wrote logEnergyFlow to file')
%         xlswrite(workSpaceName,logmdot','mdot',['A1:A' num2str(length(logmdot))]);
%         disp('Wrote logmdot to file')
        
    %     disp('Wrote fluid temp change to file')
        disp(['---> Successfully wrote temperature histories to file: ' workSpaceName])


    catch
        disp('Error writing to file')
    end
end
%}


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
plot(mdot,'LineWidth',2)
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
end
%}
%% Average temperature in tank Plot
if plotAverage
figure
hold on

plot(tempFluid(plotTimeStart:plotTimeEnd,R0),'Color',color1,'LineWidth',3)
plot((tempStorage(plotTimeStart:plotTimeEnd,R1)+tempFluid(plotTimeStart:plotTimeEnd,R1))/2,'Color',color2,'LineWidth',2)
plot((tempStorage(plotTimeStart:plotTimeEnd,R2)+tempFluid(plotTimeStart:plotTimeEnd,R2))/2,'Color',color3,'LineWidth',2)
plot((tempStorage(plotTimeStart:plotTimeEnd,R3)+tempFluid(plotTimeStart:plotTimeEnd,R3))/2,'Color',color4,'LineWidth',2)
plot((tempStorage(plotTimeStart:plotTimeEnd,R4)+tempFluid(plotTimeStart:plotTimeEnd,R4))/2,'Color',color5,'LineWidth',2)

% plot(TmeltStart*ones(1,plotTimeEnd),'LineWidth',2)

if compareToValidation
    % Zanganeh model
    plot(timeStepInverse*3600*validationDataRocks(:,1),validationDataRocks(:,2),':','LineWidth',2,'Color',color1)
    plot(timeStepInverse*3600*validationDataRocks(:,3),validationDataRocks(:,4),':','LineWidth',2,'Color',color2)
    plot(timeStepInverse*3600*validationDataRocks(:,5),validationDataRocks(:,6),':','LineWidth',2,'Color',color3)
    plot(timeStepInverse*3600*validationDataRocks(:,7),validationDataRocks(:,8),':','LineWidth',2,'Color',color4)
    plot(timeStepInverse*3600*validationDataRocks(:,9),validationDataRocks(:,10),':','LineWidth',2,'Color',color5)
    % Zanganeh experimental
    errorbarfluidinlet = 0.197530864 + (0.039407407 + 0.001986762*(validationDataRocks(:,12)-20)) + max(2.2,0.0075*validationDataRocks(:,12));
    errorbarR1 = 0.197530864 + (0.039407407 + 0.001986762*(validationDataRocks(:,14)-20)) + max(2.2,0.0075*validationDataRocks(:,14));
    errorbarR2 = 0.197530864 + (0.039407407 + 0.001986762*(validationDataRocks(:,16)-20)) + max(2.2,0.0075*validationDataRocks(:,16));
    errorbarR3 = 0.197530864 + (0.039407407 + 0.001986762*(validationDataRocks(:,18)-20)) + max(2.2,0.0075*validationDataRocks(:,18));
    errorbarR4 = 0.197530864 + (0.039407407 + 0.001986762*(validationDataRocks(:,20)-20)) + max(2.2,0.0075*validationDataRocks(:,20));
    % max offset error (NI 9218) + max gain error (at given temp)(NI 9218) + K type thermocouple error
    
    errorbar(timeStepInverse*3600*validationDataRocks(:,11),validationDataRocks(:,12),errorbarfluidinlet,'.','LineWidth',2,'Color',color1)
    errorbar(timeStepInverse*3600*validationDataRocks(:,13),validationDataRocks(:,14),errorbarR1,'.','LineWidth',2,'Color',color2)
    errorbar(timeStepInverse*3600*validationDataRocks(:,15),validationDataRocks(:,16),errorbarR2,'.','LineWidth',2,'Color',color3)
    errorbar(timeStepInverse*3600*validationDataRocks(:,17),validationDataRocks(:,18),errorbarR3,'.','LineWidth',2,'Color',color4)
    errorbar(timeStepInverse*3600*validationDataRocks(:,19),validationDataRocks(:,20),errorbarR4,'.','LineWidth',2,'Color',color5)
end
line([chargeTime chargeTime],[0 650])

ylim(ylimTemps);
% ylim([0 700]);
xlabel('Time (hours)')
ylabel(['Temperature (' char(176) 'C)'])
% uistack(h,'bottom');

% yyaxis right
% plot(mdot,'LineWidth',2)
% if includePCM
%     ylim(ylimFlow)
% else
%     ylim(ylimFlow);
% end
% ylabel('Mass Flow (kg/s)')

       
hold off
    
xticks(timeStepInverse*3600*[1 2 3 4 5 6 7 8])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
ylabel('Temperature (C)')

legendEntry0 = '20 mm';
% legend(legendEntry1,legendEntry2,legendEntry3,legendEntry4,legendEntry5,legendEntry6,legendEntry0)

if includePCM
    title('Average temperature in tank - PCM + rocks')
else
%     title('Average temperature in tank - Rocks only')
%     title([num2str(energyRadPlateKWHCharge) '---' num2str(energyRadPlateKWHDischarge) '---' num2str(energyRadPlateKWHDischarge) '---' num2str(energyWallLossRocksDischargingKWH)])
end
   
xlim(3600*timeStepInverse*xlimHours)
% xlim(3600*timeStepInverse*[3 3.25])
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

    if saveFigures
        figureNameAveragePNG = ['..\Model Outputs\Averages\MatrixModel\averageTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
        saveas(fig,figureNameAveragePNG)
        fprintf('Wrote to file: %s\n',figureNameAveragePNG)
    end
end
%}

%% Fluid temperature in tank Plot
if plotFluid
figure
hold on

% plot((tempMatrixStorage(R0,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R0,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)
if strcmp(conditionSelection,'Trahan')
    plot(1:chargeTime,tempFluid(1:chargeTime,1),'LineWidth',2,'Color',color1)
    plot(chargeTime+1:totalTime,tempFluid(chargeTime+1:totalTime,numRows-1),'LineWidth',2,'Color',color1)
else
    plot(tempFluid(plotTimeStart:plotTimeEnd,2),'LineWidth',2,'Color',color1)
end
switch conditionSelection
    case 'Trahan'
        plot((tempFluid(plotTimeStart:plotTimeEnd,R1)+tempFluid(plotTimeStart:plotTimeEnd,R1a))/2,'LineWidth',2,'Color',color2)
        plot((tempFluid(plotTimeStart:plotTimeEnd,R2)+tempFluid(plotTimeStart:plotTimeEnd,R2a))/2,'LineWidth',2,'Color',color3)
        plot((tempFluid(plotTimeStart:plotTimeEnd,R3)+tempFluid(plotTimeStart:plotTimeEnd,R3a))/2,'LineWidth',2,'Color',color4)
        plot((tempFluid(plotTimeStart:plotTimeEnd,R4)+tempFluid(plotTimeStart:plotTimeEnd,R4a))/2,'LineWidth',2,'Color',color5)
        plot((tempFluid(plotTimeStart:plotTimeEnd,R5)+tempFluid(plotTimeStart:plotTimeEnd,R5a))/2,'LineWidth',2,'Color',color6)
    case 'TrahanPCM'
%         plot((tempFluid(R1,plotTimeStart:plotTimeEnd)+tempFluid(R1a,plotTimeStart:plotTimeEnd))/2,'LineWidth',2,'Color',color2)
        plot((tempFluid(plotTimeStart:plotTimeEnd,R2)+tempFluid(plotTimeStart:plotTimeEnd,R2a))/2,'LineWidth',2,'Color',color3)
        plot((tempFluid(plotTimeStart:plotTimeEnd,R3)+tempFluid(plotTimeStart:plotTimeEnd,R3a))/2,'LineWidth',2,'Color',color4)
        plot((tempFluid(plotTimeStart:plotTimeEnd,R4)+tempFluid(plotTimeStart:plotTimeEnd,R4a))/2,'LineWidth',2,'Color',color5)
%         plot((tempFluid(R5,plotTimeStart:plotTimeEnd)+tempFluid(R5a,plotTimeStart:plotTimeEnd))/2,'LineWidth',2,'Color',color6)
end

if compareToValidation
    switch conditionSelection
        case 'Zanganeh2015'
            plot(timeStepInverse*3600*validationDataRocks(:,1),validationDataRocks(:,2),'--','LineWidth',2,'Color',color1)
            plot(timeStepInverse*3600*validationDataRocks(:,5),validationDataRocks(:,6),'--','LineWidth',2,'Color',color2)
            plot(timeStepInverse*3600*validationDataRocks(:,9),validationDataRocks(:,10),'--','LineWidth',2,'Color',color3)
            plot(timeStepInverse*3600*validationDataRocks(:,13),validationDataRocks(:,14),'--','LineWidth',2,'Color',color4)
            plot(timeStepInverse*3600*validationDataRocks(:,17),validationDataRocks(:,18),'--','LineWidth',2,'Color',color5)
        case 'Trahan' % experimental data?
            plot(timeStepInverse*3600*validationDataRocks(:,1),validationDataRocks(:,2),'-.','LineWidth',1,'Color',color1)
            plot(timeStepInverse*3600*validationDataRocks(:,3),validationDataRocks(:,4),'o','LineWidth',1,'Color',color2)
            plot(timeStepInverse*3600*validationDataRocks(:,5),validationDataRocks(:,6),'o','LineWidth',1,'Color',color3)
            plot(timeStepInverse*3600*validationDataRocks(:,7),validationDataRocks(:,8),'o','LineWidth',1,'Color',color4)
            plot(timeStepInverse*3600*validationDataRocks(:,9),validationDataRocks(:,10),'o','LineWidth',1,'Color',color5)
            plot(timeStepInverse*3600*validationDataRocks(:,11),validationDataRocks(:,12),'o','LineWidth',1,'Color',color6)
        case 'TrahanPCM' % experimental data
            uppererrorR2 = validationDataRocks(:,2) - validationDataRocks(:,8);
            lowererrorR2 = validationDataRocks(:,14) - validationDataRocks(:,2);
            
            uppererrorR3 = validationDataRocks(:,4) - validationDataRocks(:,10);
            lowererrorR3 = validationDataRocks(:,16) - validationDataRocks(:,4);
            
            uppererrorR4 = validationDataRocks(:,6) - validationDataRocks(:,12);
            lowererrorR4 = validationDataRocks(:,18) - validationDataRocks(:,6);
            
            errorbar(timeStepInverse*60*validationDataRocks(:,1),validationDataRocks(:,2),lowererrorR2,uppererrorR2,'o','LineWidth',1,'Color',color3)
            errorbar(timeStepInverse*60*validationDataRocks(:,3),validationDataRocks(:,4),lowererrorR3,uppererrorR3,'o','LineWidth',1,'Color',color4)
            errorbar(timeStepInverse*60*validationDataRocks(:,5),validationDataRocks(:,6),lowererrorR4,uppererrorR4,'o','LineWidth',1,'Color',color5)
            
%             errorbar(timeStepInverse*3600*validationDataPCM(:,9),validationDataPCM(:,10),errorbarpcm3,'o','LineWidth',2,'Color',color3)
        otherwise 
    end
    
end
    
if includePCM
plot(TmeltStart*ones(1,plotTimeEnd),'LineWidth',2)
end

ylim(ylimTemps);
ylabel(['Temperature (' char(176) 'C)'])
% uistack(h,'bottom');

hold off

switch conditionSelection
    case 'TrahanPCM'
        xlabel('Time (min)')
        xticks(timeStepInverse*60*[0:20:220])
        xticklabels({'0' '20' '40' '60' '80' '100' '120' '140' '160' '180' '200' '220'})
    otherwise
        xlabel('Time (hours)')
        xticks(1*timeStepInverse*3600*[1 2 3 4 5 6 7 8])
        xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
end
%ylabel('Temperature (C)')

legendEntry0 = '20 mm';
% legend(legendEntry1,legendEntry2,legendEntry3,legendEntry4,legendEntry5,legendEntry6,legendEntry0)

if includePCM
    title('Fluid temperature in tank - PCM + rocks')
else
    title('Fluid temperature in tank - Rocks only')
%     title([num2str(energyRadPlateKWHCharge) '---' num2str(energyRadPlateKWHDischarge) '---' num2str(energyRadPlateKWHDischarge) '---' num2str(energyWallLossRocksDischargingKWH)])
end
   
xlim(3600*timeStepInverse*xlimHours)
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

if saveFigures
    figureNameFluidPNG = ['..\Model Outputs\Fluid\fluidTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
    saveas(fig,figureNameFluidPNG)
    fprintf('Wrote to file: %s\n',figureNameFluidPNG)
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

end
%}



%% Storage temperature in tank
% {
if plotStorage
figure
hold on

% plot((tempMatrixStorage(R0,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R0,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)

switch conditionSelection
    case 'Trahan'
        plot((tempStorage(plotTimeStart:plotTimeEnd,R1)+tempStorage(plotTimeStart:plotTimeEnd,R1a))/2,'LineWidth',2,'Color',color2)
        plot((tempStorage(plotTimeStart:plotTimeEnd,R2)+tempStorage(plotTimeStart:plotTimeEnd,R2a))/2,'LineWidth',2,'Color',color3)
        plot((tempStorage(plotTimeStart:plotTimeEnd,R3)+tempStorage(plotTimeStart:plotTimeEnd,R3a))/2,'LineWidth',2,'Color',color4)
        plot((tempStorage(plotTimeStart:plotTimeEnd,R4)+tempStorage(plotTimeStart:plotTimeEnd,R4a))/2,'LineWidth',2,'Color',color5)
        plot((tempStorage(plotTimeStart:plotTimeEnd,R5)+tempStorage(plotTimeStart:plotTimeEnd,R5a))/2,'LineWidth',2,'Color',color6)
    case 'Zanganeh2015'
        plot(tempStorage(plotTimeStart:plotTimeEnd,R1),'LineWidth',2,'Color',color2)
        plot(tempStorage(plotTimeStart:plotTimeEnd,R2),'LineWidth',2,'Color',color3)
        plot(tempStorage(plotTimeStart:plotTimeEnd,R3),'LineWidth',2,'Color',color4)
        plot(tempStorage(plotTimeStart:plotTimeEnd,R4),'LineWidth',2,'Color',color5)
    case 'Beasley1989-high'
        plot(tempStorage(plotTimeStart:plotTimeEnd,R1),'LineWidth',2,'Color',color2)
        plot(tempStorage(plotTimeStart:plotTimeEnd,R2),'LineWidth',2,'Color',color3)
    case 'Beasley1989-low'
        plot(tempStorage(plotTimeStart:plotTimeEnd,R1),'LineWidth',2,'Color',color2)
        plot(tempStorage(plotTimeStart:plotTimeEnd,R2),'LineWidth',2,'Color',color3)
end

if compareToValidation
    switch conditionSelection
        case 'Trahan'
            plot(timeStepInverse*3600*validationDataRocks(:,13),validationDataRocks(:,14),'o','LineWidth',1,'Color',color2)
            plot(timeStepInverse*3600*validationDataRocks(:,15),validationDataRocks(:,16),'o','LineWidth',1,'Color',color3)
            plot(timeStepInverse*3600*validationDataRocks(:,17),validationDataRocks(:,18),'o','LineWidth',1,'Color',color4)
            plot(timeStepInverse*3600*validationDataRocks(:,19),validationDataRocks(:,20),'o','LineWidth',1,'Color',color5)
            plot(timeStepInverse*3600*validationDataRocks(:,21),validationDataRocks(:,22),'o','LineWidth',1,'Color',color6)
%             plot(timeStepInverse*3600*validationDataRocks(:,11),validationDataRocks(:,12),'o','LineWidth',1,'Color',color6)
    case 'Beasley1989-high'
            plot(timeStepInverse*3600*validationDataRocks(:,1),validationDataRocks(:,2),'o','LineWidth',1,'Color',color2)
            plot(timeStepInverse*3600*validationDataRocks(:,3),validationDataRocks(:,4),'o','LineWidth',1,'Color',color3)
    case 'Beasley1989-low'
            plot(timeStepInverse*3600*validationDataRocks(:,1),validationDataRocks(:,2),'o','LineWidth',1,'Color',color2)
            plot(timeStepInverse*3600*validationDataRocks(:,3),validationDataRocks(:,4),'o','LineWidth',1,'Color',color3)
    end
    
end

ylim(ylimTemps);
xlabel('Time (hours)')
ylabel(['Temperature (' char(176) 'C)'])
% uistack(h,'bottom');
       
hold off
    
xticks(1*timeStepInverse*3600*[1 2 3 4 5 6 7 8])
xticklabels({'1' '2' '3' '4' '5' '6' '7' '8'})
%ylabel('Temperature (C)')

legendEntry0 = '20 mm';
% legend(legendEntry1,legendEntry2,legendEntry3,legendEntry4,legendEntry5,legendEntry6,legendEntry0)

title('Storage temperature in tank - Rocks only')
   
xlim(3600*timeStepInverse*xlimHours)
fig = gcf;
set(fig, 'Units', 'Pixels', 'OuterPosition', [650, 200, 630, 395]);
fig.PaperPositionMode = 'auto';

if saveFigures
    figureNameStoragePNG = ['..\Model Outputs\Storage\storageTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
    saveas(fig,figureNameStoragePNG)
    fprintf('Wrote to file: %s\n',figureNameStoragePNG)
end
end
%}

%% PCM temperature in tank
% {
if plotPCM && strcmp('Zanganeh2015',conditionSelection)
    if includePCM
        figure
        hold on

        % plot((tempMatrixStorage(R0,plotTimeStart:plotTimeEnd)+tempMatrixFluid(R0,plotTimeStart:plotTimeEnd))/2,'LineWidth',2)
        green = [0.4660 0.6740 0.1880]; blue = [0 0.4470 0.7410];
        red = [0.8500 0.3250 0.0980]; yellow = [0.9290 0.6940 0.1250];
        plot(tempStorage(plotTimeStart:plotTimeEnd,2),'LineWidth',2,'Color',green)
        plot(tempFluid(plotTimeStart:plotTimeEnd,2),'LineWidth',2,'Color',blue)
        plot(tempStorage(plotTimeStart:plotTimeEnd,4),'LineWidth',2,'Color',red)
        plot(tempFluid(plotTimeStart:plotTimeEnd,4),'LineWidth',2,'Color',yellow)
        
        if compareToValidation
            % Zanganeh model
%             plot(timeStepInverse*3600*validationDataPCM(:,1),validationDataPCM(:,2),':','LineWidth',2,'Color',green)
%             plot(timeStepInverse*3600*validationDataPCM(:,3),validationDataPCM(:,4),':','LineWidth',2,'Color',blue)
%             plot(timeStepInverse*3600*validationDataPCM(:,5),validationDataPCM(:,6),':','LineWidth',2,'Color',red)
%             plot(timeStepInverse*3600*validationDataPCM(:,7),validationDataPCM(:,8),':','LineWidth',2,'Color',yellow)

            % Zanganeh experimental
            errorbarpcm3 = 0.197530864 + (0.039407407 + 0.001986762*(validationDataPCM(:,10)-20)) + max(2.2,0.0075*validationDataPCM(:,10));
            errorbarfluid3 = 0.197530864 + (0.039407407 + 0.001986762*(validationDataPCM(:,12)-20)) + max(2.2,0.0075*validationDataPCM(:,12));
            errorbarpcm5 = 0.197530864 + (0.039407407 + 0.001986762*(validationDataPCM(:,14)-20)) + max(2.2,0.0075*validationDataPCM(:,14));
            errorbarfluid5 = 0.197530864 + (0.039407407 + 0.001986762*(validationDataPCM(:,16)-20)) + max(2.2,0.0075*validationDataPCM(:,16));
            
            errorbar(timeStepInverse*3600*validationDataPCM(:,9),validationDataPCM(:,10),errorbarpcm3,'o','LineWidth',2,'Color',green)
            errorbar(timeStepInverse*3600*validationDataPCM(:,11),validationDataPCM(:,12),errorbarfluid3,'o','LineWidth',2,'Color',blue)
            errorbar(timeStepInverse*3600*validationDataPCM(:,13),validationDataPCM(:,14),errorbarpcm5,'o','LineWidth',2,'Color',red)
            errorbar(timeStepInverse*3600*validationDataPCM(:,15),validationDataPCM(:,16),errorbarfluid5,'o','LineWidth',2,'Color',yellow)
        end
        line([chargeTime chargeTime],[0 650])


        ylim([500 700]);
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

        if saveFigures
            figureNamePCMPNG = ['..\Model Outputs\PCM\PCMTemp' datestr(fileNameTime,'yyyy-mm-dd-HH-MM-SS') '.png'];
            saveas(fig,figureNamePCMPNG)
            fprintf('Wrote to file: %s\n',figureNamePCMPNG)
        end
    end
end
%}

%% Plot temp gradient at each time sequentially
%{
% transVis = true; %transVis;
% transVis = input('Figure loaded. Press 1 to run transient visualization, or any other number to skip: ');
if transVis

    figure
    tic
    startRow = 3;

    endRow = numRows-1;
    for t = 1:timeStepInverse*60*0.25:chargeTime
        % % f
%         t= logTime

        clf
%         figure
        hold on
        yyaxis right
        ylim([0 100])
        plot(abs(tempFluid(startRow:endRow,t)-tempStorage(startRow:endRow,t)),'LineWidth',3)

        yyaxis left
        ylim(ylimTemps)
        plot(tempFluid(startRow:endRow,t),'LineWidth',3)
        plot(tempStorage(startRow:endRow,t),'LineWidth',3)
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

%% Plot mass flow gradient at each time sequentially
%{
% transVis = true; %transVis;
% transVis = input('Figure loaded. Press 1 to run transient visualization, or any other number to skip: ');
if transVisMdot

    figure
    tic
    startRow = 3;

    endRow = numRows-1;
    for t = 1:timeStepInverse*60*0.25:chargeTime
        % % f
%         t= logTime

        clf
%         figure
        hold on
        yyaxis left
        plot(mdot(startRow:endRow,t),'LineWidth',3)
        ylim([min(mdot(:,t))-0.00001 max(mdot(:,t))+0.00001])
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

%% notify of ending
fprintf(['\nEnded at time ' datestr(clock,'yyyy/mm/dd--HH:MM:SS') '\n'])