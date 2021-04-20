function [tempsFluid] = calculateTempFluidCP(intEnergyFluid,intEnergyFluidPrev,tempsFluid, specificHeat,~,~)
% intEnergyMatrix,includePCM,fluidPressure,RH
%     numRows = length(intEnergyFluid);
        
    %-------- air temps
    tempsFluid = tempsFluid+(intEnergyFluid - intEnergyFluidPrev)./specificHeat;
    
%     for i = 1:1:numRows
%         
%         if ~isreal(intEnergyFluid(i))
%             disp(intEnergyFluid(i))
%             MEx = MException('calculateTemps:FluidEnergyImaginary',...
%                 'Fluid energy imaginary (%g J/kg) at row %g',intEnergyFluid(i),i);
%             throw(MEx)
%         end
%         
%         energyChange = intEnergyFluid(i) - intEnergyFluidPrev(i);
% %         disp(energyChange)
% %         disp(specificHeat(i))
%         temp = tempsFluid(i)+energyChange/specificHeat(i);    
%           
%         if temp < 24.95
%             MEx = MException('calculateTemps:FluidTempOutOfRange',...
%                 'Calculated Temp out of bounds (%g C)\n   at row %g for energy %g J/kg',temp,i,intEnergyFluid(i));
%             throw(MEx)
%         end
%         if ~isreal(temp)
%             disp(temp)
%             MEx = MException('calculateTemps:FluidTempImaginary',...
%                 'Calculated temp imaginary (%g C) at row %g for energy %g',temp,i,intEnergyFluid(i));
%             throw(MEx)
%         end
%         tempsFluid(i) = temp;
%     end
%     disp(tempMatrixFluid)
end