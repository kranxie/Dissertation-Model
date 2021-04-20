function [exergy] = getExergy(temp,deadTemp,getEnthalpy,getEntropy)

    exergy = getEnthalpy(temp) - getEnthalpy(deadTemp) - (deadTemp+273)*(getEntropy(temp)-getEntropy(deadTemp));
end