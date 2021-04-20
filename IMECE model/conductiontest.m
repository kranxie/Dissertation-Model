clc
clear all
close all

baseTemp = 25;
topTemp = 650;
intEnergyRocks(1:70,2) = 0;
timeStepInverse = 40;
totalTime = 3*3600*timeStepInverse;
onePercentTime = round(0.01*totalTime,0);
tempMatrixStorage(1:70,totalTime) = 25;
timeStep = 1/timeStepInverse;
kFluid = 0.0365;
voidFractionRocks = 0.4;
voidFractionRocksInverse = 1-voidFractionRocks;
rockAverageDiameter = 0.032;
numRows = 70;
emissivityRocks = 0.83;
rowSurfaceAreaRocks = 0.2743;
densityRocks = 2635;
rowLengthRock = 0.02;
rowVolumeRock = 0.0024;
% intEnergyMatrixFluid(70,2) = 0;
% intEnergyMatrixFluid(:,1) = getEnthalpyDryAir(30);
% tempMatrixFluid(70,1) = 0;
% tempMatrixFluid(:) = 30;
getEnthalpy = @getEnthalpyDryAir;
getEnergyDerivative = @getEnergyDerivativeDryAir;
getkRocks = @(x) 5+(1-5)*(x-25)/(700-25);
getRockEnergy = @(temp) 747.0995*temp+0.2838*temp.^2;
calculateTemps = @calculateTempsRocksNoFluid;

% intEnergyRocks(1:numRows,1) = getRockEnergy(baseTemp);
% intEnergyRocks(1:2,1) = getRockEnergy(topTemp);
for i = 2:numRows-1
    intEnergyRocks(i,1) = getRockEnergy(topTemp+(i-1)/(numRows-2)*(baseTemp-topTemp));
end
logConduction(totalTime) = 0;

tic
for t = 1:totalTime
%     tempMatrixStorage(1:2,t) = topTemp;
%     intEnergyRocks(1:2,1) = getRockEnergy(topTemp);
    
    kRocks = getkRocks(tempMatrixStorage(1:numRows,t));
    [keffRocks,~] = calculatekeff(kRocks,kFluid,tempMatrixStorage(:,t),tempMatrixStorage(:,t),voidFractionRocks,rockAverageDiameter,numRows,emissivityRocks);
    for i = 2:69
        conduction = getRockConduction(keffRocks(i),tempMatrixStorage(i-1,t),tempMatrixStorage(i,t),tempMatrixStorage(i+1,t),timeStep,rowSurfaceAreaRocks,voidFractionRocksInverse,densityRocks,rowLengthRock,rowVolumeRock);
        intEnergyRocks(i,2) = intEnergyRocks(i,1) + conduction;
        logConduction(t) = logConduction(t) + conduction;
    end
    intEnergyRocks(1,2) = intEnergyRocks(2,2);
    intEnergyRocks(70,2) = intEnergyRocks(69,2);
        
    [tempMatrixStorage(:,t+1)] = calculateTemps(intEnergyRocks(:,2),tempMatrixStorage(:,t));
    intEnergyRocks(:,1) = intEnergyRocks(:,2);
    
    if mod(t,onePercentTime) == 0
        clc
        timeSoFar = toc;
        fprintf('Progress: %g %%\n',round(t/totalTime*100,0))
        fprintf('Time elapsed so far: %g min\n',timeSoFar/60)
        fprintf('Estimated time to completion: %g min\n',timeSoFar*(totalTime/t-1)/60)
    end
end


%% Plot temp gradient at each time sequentially
figure
transVis = input('Figure loaded. Press 1 to run transient visualization, or any other number to skip: ');
if transVis == 1
    tic
    startRow = 1;
    endRow = 70; %numRows
    for t = 1:timeStepInverse*60*1:totalTime
        clf
        hold on
        yyaxis right
        ylim([0 100])
%         plot(abs(tempMatrixFluid(startRow:endRow,t)-tempMatrixStorage(startRow:endRow,t)),'LineWidth',3)

        yyaxis left
        ylim([25 705])
%         plot(tempMatrixFluid(startRow:endRow,t),'LineWidth',3)
        plot(tempMatrixStorage(startRow:endRow,t),'LineWidth',3)
%         plot([2 2], [0 655])
%         plot([3 3], [0 655])
        plot([6 6], [0 705])


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