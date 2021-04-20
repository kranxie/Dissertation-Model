clear all
addpath(genpath('fluid properties'))

conditionSelection = 'TrahanPCM'; % or 'Trahan'
includePCM = true;

switch conditionSelection
    case 'Trahan'
        fileName = 'inputsTrahan.mat';
        chargeTimeHours = 1.5;                % 3 is standard\
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 1.5;             % 5 standard
    case 'TrahanPCM'
        fileName = 'inputsTrahanPCM.mat';
        chargeTimeHours = 1; %3 + 2/3;                % 3 is standard\
        chargeTimePartialHour = 0;        % 1/6 for Zanganeh, 0 for Beasley
        dischargeTimeHours = 0; %3+1/3;             % 5 standard
    otherwise
        disp('Invalid condition selection')
        return
end


%% ============== Tank dimensions

switch conditionSelection
    case 'Trahan'
        tankHeight = 0.889;
        tankDiameter = 0.445;
        tankRadius = tankDiameter / 2;
        rockHeight = tankHeight;   
    case 'TrahanPCM'
        tankHeight = 0.254;
        tankDiameter = 0.254;
        tankRadius = tankDiameter / 2;
        rockHeight = 0;
        heightPCM = tankHeight;
end

tankArea = pi*(tankRadius)^2;               % m2

%% ============== Insulation parameters
switch conditionSelection
    case 'Trahan'
        tInsulationInside = 0.133; %m
        tInsulationOutside = 50.8/1000; %m
        tWall = 3.175/1000; %m
        
        kInsulation = 0.06; %W/mK
        kWall = 43; %W/mK based on 1% carbon steel
        
    case 'TrahanPCM'
        tInsulationInside = 0; %m
        tInsulationOutside = 0.1524; %m
        tWall = 0.00635; %m
        
        kInsulation = 0.06; %W/mK
        kWall = 43.84; %W/mK based on 1% carbon steel
    otherwise
        lossInsulationTerms(1:numRows) = 0;
end

%% ============== Rock properties
% {
switch conditionSelection
    case 'Trahan'
        getkRocks = @(x) 1.5;
        densityRocks = 3200;                                        % kg/m3
        rockAverageDiameter = 0.0425;                                % m
        rockAverageRadius = rockAverageDiameter/2;
        rockAverageSurfaceArea = 4*pi*(rockAverageRadius)^2;    % m2
        rockAverageVolume = 4/3*pi*(rockAverageRadius)^3;       % m3
        
        voidFractionRocks = 0.51;                                    % .4
        voidFractionRocksInverse = 1-voidFractionRocks;
    case 'TrahanPCM'
        getkRocks = @(x) 0.52;
        densityRocks = 0;                                        % kg/m3
        PCMAverageDiameter = 0.02653;                                % m
        rockAverageDiameter = 0.02653 + 2*0.00045;                                % m
        rockAverageRadius = rockAverageDiameter/2;
        dEncap = rockAverageDiameter;
        rEncap = rockAverageDiameter/2;
        rockAverageSurfaceArea = 4*pi*(rockAverageRadius)^2;    % m2
        rockAverageVolume = 4/3*pi*(rockAverageRadius)^3;       % m3
        surfaceAreaPCM = rockAverageSurfaceArea;
        
        storageWallHXareaTerms = 0;
        fcontPCM = 1;
        
        voidFractionRocks = 1;                              % .4
        voidFractionRocksInverse = 1-voidFractionRocks;
        voidFractionPCM = 0.348;                              % .4
        voidFractionPCMInverse = 1-voidFractionPCM;
end
emissivityRocks = 0.83;
hvRocksCoeff = (1/rockAverageDiameter)^0.92 * 824;


%% ============== Metal Plate properties
% {
voidFractionPlate = 0.142;
voidFractionPlateInverse = 1- voidFractionPlate;
thicknessPlate = 0.02;          % m
kPlate = 38;                    % W/mK 41-35
plateArea = tankArea;
emissivityPlate = 0.7;
densityPlate = 7930;
surfaceAreaPlate = tankArea*(1-voidFractionPlate)*2+221*pi*0.01*0.02;
%}


%% PCM/encapsulation Properties

switch conditionSelection
    case 'Trahan'
    case 'TrahanPCM'
        densityPCMLiquid = 2125;
        densityPCMSolid = 1908;
%         densityPCMEff = (densityPCMLiquid+densityPCMSolid)/2; % * (PCMAverageDiameter/rockAverageDiameter)^2;
        densityPCMEff = densityPCMLiquid; % * (PCMAverageDiameter/rockAverageDiameter)^2;
        
        Tmelt = 306;
        dTMelt = 2;
        TmeltStart = Tmelt - dTMelt/2;
        TmeltEnd = Tmelt + dTMelt/2;
%         TmeltStart = 299.84;
%         TmeltEnd = 310.68;
%         dTMelt = TmeltEnd - TmeltStart;
        TsolidifyStart = 295.2;
        TsolidifyEnd = 302.34;
        
        cPCMSolid = 1835.4;
        cPCMLiquid = 1655;
        enthalpyFusion = 172000;
        cPCMPhaseChange = enthalpyFusion / dTMelt + (cPCMSolid+cPCMLiquid)/2;
        cPCMMeltStart = @(T) (6.138410916 * T.^4) + (8.184093887e3*T.^3) - (4.054467609e6 *T.^2) + (8.859274628e8*T) - 7.212555798e10; % 299.84 to 306.35 C
        cPCMMeltEnd = @(T) (40.885751986 *T.^4) - (5.1536031432e4 *T.^3) + (2.4355700734e7 *T.^2) - (5.1148076097e9 *T) - 4.0272892760e11; % 306.35 to 310.68 C
        cPCMSolidify = @(T) (-11.7992675546556 *T.^6) + (2.10409802521218e4 *T.^5) - (1.5633494286549e7 *T.^4) + (6.19491471878624e9 *T.^3) - (1.38078998166774e12 *T.^2) + (1.64138018646662e14 *T) - 8.12961689535658e15; % 295.2 to 302. 34 C
                
        kPCMLiquid = 0.54;
        kPCMSolid = 0.50;
        kEncaps = 0.22;
        
        ceffPCMLiquid = cPCMLiquid;
        ceffPCMSolid = cPCMSolid;
        ceffPCMPhaseChange = cPCMPhaseChange;
        pcmLowEnergyCutoff = TmeltStart*ceffPCMSolid;
        pcmHighEnergyCutoff = pcmLowEnergyCutoff + ceffPCMPhaseChange*dTMelt;
end

%% Input fluid and top of tank conduction
switch conditionSelection
    case 'Trahan'
        tempProfile = @tempProfileTrahan;
        getmdot = @getmdotTrahan;
        getStorageTopConduction = @getRockConduction;
    case 'TrahanPCM'
        getStorageTopConduction = @getRockConduction;
        tempProfile = @tempProfileTrahanPCM;
        getmdot = @getmdotTrahanPCM;
        calculatehvPCM = @calculateHvPCMBeasley;
end
%%

% if includePCM
%     thicknessEncapsulation = 0; %0.001; % m
%     emissivityEncapsulation = 0.7;

%     Crow = 0.95;
%     densityEncap = 7930;                    % kg/m3
%     enthalpyFusion = 197.71*1000; %466*1000;                   % J/kg
%     TmeltStart = 47.13; %573 ; %+ fluidTempIncreaseAmount;
%     dTMelt = 5; %4;                             % K
%     TmeltEnd = TmeltStart + dTMelt;
    
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
switch conditionSelection
    case 'Trahan'
        dischargeTemp = 192;                                 % C - temperature of HTF at bottom inlet during discharge
        baseTemp = 25;                                      % C - starting temperature of rock bed
        ambientTemp = 25;
        maxTemp = 500;
        avgTemp = (baseTemp + maxTemp)/2;
    case 'TrahanPCM'
        dischargeTemp = 286;                                 % C - temperature of HTF at bottom inlet during discharge
        baseTemp = 285;                                      % C - starting temperature of rock bed
        ambientTemp = 25;
        maxTemp = 326;
        avgTemp = (baseTemp + maxTemp)/2;
    otherwise
        baseTemp = 26.5;
        ambientTemp = 25;
end

save(fileName)
fprintf('Saved workspace to file: %s\n',fileName)
