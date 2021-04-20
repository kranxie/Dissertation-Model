function [vdot] = volumeFlowRate(t,chargeTime,timeStepInverse,bypassFractionInverse,temp,getmdot,getRefDensity)
    mdot = getmdot(t,chargeTime,timeStepInverse,bypassFractionInverse);
    density = getRefDensity(temp);
    vdot = mdot / density;              % m3/s
end