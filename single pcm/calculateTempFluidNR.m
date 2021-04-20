function [tempFluid] = calculateTempFluidNR(intEnergyFluid,~, tempFluid, ~,energyFunction,specificHeatFunction)
% intEnergyMatrix,includePCM,fluidPressure,RH
        
    % air temps
    if ~isreal(intEnergyFluid)
        disp(intEnergyFluid)
        MEx = MException('calculateTemps:FluidEnergyImaginary',...
            'Fluid energy imaginary (%g J/kg)',intEnergyFluid);
        throw(MEx)
    end
        
    testEnergy = energyFunction(tempFluid) - intEnergyFluid;
    energyDerivative = specificHeatFunction(tempFluid);
    tempFluid = tempFluid - testEnergy/energyDerivative;
    error = abs(testEnergy);
    
    while error > 0.00001
% 		disp(error)
%             disp('Error in energy solver > 1000, running extra cycle. Row:')
%             disp(i)
        testEnergy = energyFunction(tempFluid) -intEnergyFluid;
        energyDerivative = specificHeatFunction(tempFluid);
        tempFluid = tempFluid - testEnergy/energyDerivative;
        error = abs(testEnergy);
%             disp(i)
%             disp(error)
%             if error > 1e6
%                 MEx = MException('calculateTemps:FluidEnergyOutOfRange',...
%                     'Temp solver error too high at row %g\n   for airEnergy: %g J/kg',i,intEnergyMatrixFluid(i));
%                 throw(MEx)
%                 return
%             end
    end
    if tempFluid < 24.95
        MEx = MException('calculateTemps:FluidTempOutOfRange',...
            'Calculated Temp out of bounds (%g C)\n   for energy %g J/kg',tempFluid,intEnergyFluid);
%             'Calculated Temp out of bounds (%g C)\n   at row %g for energy %g J/kg',temp,i,intEnergyFluid(i));
        throw(MEx)
    end
    if ~isreal(tempFluid)
        disp(tempFluid)
        MEx = MException('calculateTemps:FluidTempImaginary',...
            'Calculated temp imaginary (%g C) for energy %g',tempFluid,intEnergyFluid);
%             'Calculated temp imaginary (%g C) at row %g for energy %g',temp,i,intEnergyFluid(i));
        throw(MEx)
    end
    
end