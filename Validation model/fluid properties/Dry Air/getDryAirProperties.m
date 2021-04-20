function [specificheat, conductivity, viscosity,density,enthalpy,enthalpyDerivative] = getDryAirProperties(airPressure)
    enthalpy = @(temp)(0.00000000025104/5*temp.^5-0.0000007284/4*temp.^4+0.00066686/3*temp.^3-0.010204/2*temp.^2+1006*temp+23294.44473);
    enthalpyDerivative = @(temp)(0.00000000025104*temp.^4-0.0000007284*temp.^3+0.00066686*temp.^2-0.010204*temp+1006);
    specificheat = @(temp)0.00000000025104*temp.^4-0.0000007284*temp.^3+0.00066686*temp.^2-0.010204*temp+1006;
    conductivity = @(temp)-1.526e-8*temp.^2 + 6.976e-5*temp+2.495e-2;
    viscosity = @(temp)-1.248e-11*temp.^2+4.375e-8*temp+1.746e-5;
    density = @(temp)airPressure/287.1/(temp+273.15);
%     density = @(temp) 3.899e-12*temp.^4 - 8.854e-9*temp.^3 + 7.946e-6*temp.^2 -3.828e-3*temp + 1.238;
    % the 4th order density is the one used for 50% RH moist air, I was
    % just testing something
end