function [tempRock] = calculateTempRock(currentTemp, intEnergyStorageCurr,intEnergyStoragePrev,~,~, ~, ~, ~, ~, ~)
   % PCM
%     tempRock = (-747.0995+sqrt(747.0995^2+4*0.2838*intEnergyStorage))/2/0.2838;
%     disp(currentTemp)
%     disp(storageEnergyChange)
%     disp((747.0995+2*0.2838*currentTemp))
%     disp(storageEnergyChange/(747.0995+2*0.2838*currentTemp))
%     disp(' ')
    storageEnergyChange = intEnergyStorageCurr-intEnergyStoragePrev;
    tempRock = currentTemp + storageEnergyChange/(747.0995+0.5676*currentTemp);
%     0.2838*currentTemp);
       
    if imag(tempRock)
        MEx = MException('calculateTemps:StorageTempImaginary',...
            'Calculated temp imaginary');
        throw(MEx)
    end
end