function [mdot] = getmdotBeasley(highFlow)
    if highFlow
        mdot = 0.0328;      % From Zanganeh 2014
    else
        mdot = 0.008792;    %360/970 * 0.0328
    end
end

