% function [mdot] = getmdotNoPCM(~,~,~,~)
%     mdot = 0.01761;
%     
% end

function [mdot] = getmdotNoPCM(t,chargeTime,timeStepInverse,bypassFractionInverse)
    if t < chargeTime+0.25*3600*timeStepInverse && t > chargeTime
        mdot = bypassFractionInverse*(0.625570726660953*(t/3600/timeStepInverse)^4 - 8.35993367504099*(t/3600/timeStepInverse)^3 + 41.8884361826403*(t/3600/timeStepInverse)^2 - 93.2692029047124*(t/3600/timeStepInverse) + 77.8803399131827);
    elseif t > 3*3600*timeStepInverse
        mdot = bypassFractionInverse*0.0138;                   % 0.014
    elseif  t > 3600*timeStepInverse
        mdot = bypassFractionInverse*(0.013 + (0.0138-0.013)*(t-3600*timeStepInverse)/(2*3600*timeStepInverse));
    elseif t > 1200*timeStepInverse
        mdot = bypassFractionInverse*(0.009 + (0.013-0.009)*t/(3600*timeStepInverse));
    else
        mdot = bypassFractionInverse*(0.008 + (0.010333-0.008)*t/(1200*timeStepInverse));
    end
end