function [tempRock] = calculateTempRockZang2012NR(currentTempC,intEnergy,~, ~,energyFunction,specificHeatFunction,i)
%     storageEnergyChange = intEnergyStorageCurr-intEnergyStoragePrev;
%     if useZang2012Rocks
%         currTempK = currentTempC + 273; % T for this calc should be in Kelvin, per Kelley 1960
%         cRock = 705+0.43287*currTempK-13606500*currTempK.^(-2);
%     else
%         cRock = (747.0995+0.5676*currentTempC);
%     end

    tempK = currentTempC + 273;
    testEnergy = energyFunction(tempK) -intEnergy;
    energyDerivative = specificHeatFunction(tempK);
    tempK = tempK - testEnergy/energyDerivative;
    error = abs(testEnergy);
    while error > 0.00001
        testEnergy = energyFunction(tempK) -intEnergy;
        energyDerivative = specificHeatFunction(tempK);
        tempK = tempK - testEnergy/energyDerivative;
        error = abs(testEnergy);
        if error > 1e6
            MEx = MException('calculateTemps:FluidEnergyOutOfRange',...
                'Temp solver error too high at row %g\n   for rockEnergy: %g J/kg',i,intEnergy);
            throw(MEx)
            return
        end
    end
        
    
    tempRock = tempK - 273;
       
    if imag(tempRock)
        MEx = MException('calculateTemps:StorageTempImaginary',...
            'Calculated temp imaginary');
        throw(MEx)
    end
end