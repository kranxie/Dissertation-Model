function [tempsFluid] = calculateTempFluidNR(intEnergyFluid,~, tempsFluid, ~,energyFunction,specificHeatFunction)
% intEnergyMatrix,includePCM,fluidPressure,RH
    numRows = length(intEnergyFluid);
        
    temp = tempsFluid(1);
    % air temps
    for i = 1:1:numRows
        if ~isreal(intEnergyFluid(i))
            disp(intEnergyFluid(i))
            MEx = MException('calculateTemps:FluidEnergyImaginary',...
                'Fluid energy imaginary (%g J/kg) at row %g',intEnergyFluid(i),i);
            throw(MEx)
        end
        
        testEnergy = energyFunction(temp) -intEnergyFluid(i);
        energyDerivative = specificHeatFunction(temp);
        temp = temp - testEnergy/energyDerivative;
        error = abs(testEnergy);
        while error > 0.00001
%             disp('Error in energy solver > 1000, running extra cycle. Row:')
%             disp(i)
            testEnergy = energyFunction(temp) -intEnergyFluid(i);
            energyDerivative = specificHeatFunction(temp);
            temp = temp - testEnergy/energyDerivative;
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
        if temp < 24.95
            MEx = MException('calculateTemps:FluidTempOutOfRange',...
                'Calculated Temp out of bounds (%g C)\n   at row %g for energy %g J/kg',temp,i,intEnergyMatrixFluid(i));
            throw(MEx)
        end
        if ~isreal(temp)
            disp(temp)
            MEx = MException('calculateTemps:FluidTempImaginary',...
                'Calculated temp imaginary (%g C) at row %g for energy %g',temp,i,intEnergyMatrixFluid(i));
            throw(MEx)
        end
        tempsFluid(i) = temp;
    end
    
end