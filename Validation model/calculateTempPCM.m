function [tempPCM] = calculateTempPCM(currentTemp,intEnergyStorageCurr,intEnergyStoragePrev, cPCM)
   % PCM
   tempPCM = currentTemp + (intEnergyStorageCurr-intEnergyStoragePrev)./cPCM;
   
    
    if imag(tempPCM)
        MEx = MException('calculateTemps:StorageTempImaginary',...
            'Calculated temp imaginary');
        throw(MEx)
    end
end