function calculateFluidTempSingleRow(energy,fluidTestEnergy,fluidEnergyDerivative)
    i = 1;
    temp = 650;
    testEnergy = fluidTestEnergy(temp) - energy;
    energyDerivative = fluidEnergyDerivative(temp);
    temp = temp - testEnergy/energyDerivative;
    error = abs(testEnergy);
    while error > 1e1
        i = i+1;
%             disp('Error in energy solver > 1000, running extra cycle. Row:')
        testEnergy = fluidTestEnergy(temp) -energy;
        energyDerivative = fluidEnergyDerivative(temp);
        temp = temp - testEnergy/energyDerivative;
        error = abs(testEnergy);
        if error > 1e6
            MEx = MException('calculateTemps:FluidEnergyOutOfRange',...
                'Temp solver error too high for airEnergy: %g J/kg',energy);
            throw(MEx)
            return
        end
    end
    if temp < 24.95
        MEx = MException('calculateTemps:FluidTempOutOfRange',...
            'Calculated Temp out of bounds (%g C)\n   at row %g',temp);
        throw(MEx)
        return
    end
    fprintf('Energy: %g J/kg \nTemp: %g C \nCycles to converge: %g\n',energy,temp,i)
end