function [mdot] = getmdotWithPCM(t,chargeTime,timeStepInverse,bypassFractionInverse)
    if t > 3*3600*timeStepInverse
        mdot = bypassFractionInverse*0.0138;                   % 0.014
    elseif  t > 3600*timeStepInverse
        mdot = bypassFractionInverse*(0.013 + (0.0138-0.013)*(t-3600*timeStepInverse)/(2*3600*timeStepInverse));
    elseif t > 1200*timeStepInverse
        mdot = bypassFractionInverse*(0.009 + (0.013-0.009)*t/(3600*timeStepInverse));
    else
        mdot = bypassFractionInverse*(0.008 + (0.010333-0.008)*t/(1200*timeStepInverse));
    end
end