function [tempRock] = calculateTempRockZang2012CP(currentTempC, intEnergyCurr,intEnergyPrev,useZang2012Rocks, ~, ~,i)
    storageEnergyChange = intEnergyCurr-intEnergyPrev;
    if useZang2012Rocks
        currTempK = currentTempC + 273; % T for this calc should be in Kelvin, per Kelley 1960
        cRock = 705+0.43287*currTempK-13606500*currTempK.^(-2);
    else
        cRock = (747.0995+0.5676*currentTempC);
    end
%     disp(crock)
    tempRock = currentTempC + storageEnergyChange/cRock;
%     (747.0995+0.5676*currentTemp);
       
    if imag(tempRock)
        MEx = MException('calculateTemps:StorageTempImaginary',...
            'Calculated temp imaginary at row %g',i);
        throw(MEx)
    end
end

