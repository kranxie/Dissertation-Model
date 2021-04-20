function [tempRock] = calculateTempRockQuadratic(~, intEnergyCurr,~,~, ~,~, ~, ~, ~, ~,~,~)
%     storageEnergyChange = intEnergyCurr-intEnergyPrev;
%     currTempK = currentTempC + 273; % T for this calc should be in Kelvin, per Kelley 1960
%     cRock = getRockSpecificHeat(currTempK);
    
    tempRock = (-747.0995+sqrt(747.0995^2+4*0.2838*intEnergyCurr))/2/0.2838;
    if imag(tempRock)
        MEx = MException('calculateTemps:StorageTempImaginary',...
            'Calculated temp imaginary');
        throw(MEx)
    end
end

