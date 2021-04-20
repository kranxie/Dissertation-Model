function [tempRock] = calculateTempRockCP(currentTempC, intEnergyCurr,intEnergyPrev,getRockSpecificHeat)
    currTempK = currentTempC + 273; % T for this calc should be in Kelvin, per Kelley 1960
    cRock = getRockSpecificHeat(currTempK);
    
    tempRock = currentTempC + (intEnergyCurr-intEnergyPrev)./cRock;
       
%     if imag(tempRock)
%         MEx = MException('calculateTemps:StorageTempImaginary',...
%             'Calculated temp imaginary at row %g',i);
%         throw(MEx)
%     end
end

