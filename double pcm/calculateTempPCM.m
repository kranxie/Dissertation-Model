function [tempPCM] = calculateTempPCM(currentTemp,intEnergyStorageCurr,intEnergyStoragePrev,lowEnergyCutoff,highEnergyCutoff, TmeltStart, TmeltEnd, cEffSolid, cEffPhaseChange, cEffLiquid,i,~,~)
   % PCM
   storageEnergyChange = intEnergyStorageCurr-intEnergyStoragePrev;
    if currentTemp > TmeltEnd       % Liquid
%         tempPCM = TmeltEnd +(intEnergyStorage-highEnergyCutoff)/cEffLiquid;
        tempPCM = currentTemp + storageEnergyChange/cEffLiquid;
    elseif currentTemp > TmeltStart      % phase change
%         tempPCM = TmeltStart+(intEnergyStorage-lowEnergyCutoff)/cEffPhaseChange;
        tempPCM = currentTemp + storageEnergyChange / cEffPhaseChange;
        if tempPCM > TmeltEnd
            tempPCM = TmeltEnd+(tempPCM-TmeltEnd)*cEffPhaseChange/cEffLiquid;
        end
    else                                    % solid
%         tempPCM = (intEnergyStorage)/cEffSolid;
        tempPCM = currentTemp + storageEnergyChange/cEffSolid;
        if tempPCM > TmeltStart
            tempPCM = TmeltStart+(tempPCM-TmeltStart)*cEffSolid/cEffPhaseChange;
        end
    end
    
%     if imag(tempPCM)
%         MEx = MException('calculateTemps:StorageTempImaginary',...
%             'Calculated temp imaginary');
%         throw(MEx)
%     end
end