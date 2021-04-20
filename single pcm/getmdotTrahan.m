function [mdot] = getmdotTrahan(t,chargeTime,timeStep)
    if t*timeStep < chargeTime
        mdot = 0.044015;
    else
        mdot = 0.04484;
    end
end