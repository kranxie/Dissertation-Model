function [tempsFluid] = calculateTempFluidCP(intEnergyFluid,intEnergyFluidPrev,tempsFluid, specificHeat,~,~)
% intEnergyMatrix,includePCM,fluidPressure,RH
tempsFluid = tempsFluid + (intEnergyFluid - intEnergyFluidPrev)./specificHeat;
end