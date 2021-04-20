function [specificheat, conductivity, viscosity,density,enthalpy] = getCO2Properties(temp,pressure,RH)

% Properties valid for pressures of approximately 12.38 MPa
% Temperatures in degC
% pressure not used but included for future iterations

    enthalpy = -0.4590*temp^2+1666*temp+295300;
    specificheat = -0.9180*temp+1666;
    if temp > 100
        conductivity = 5.788e-13*temp^4-1.163e-9*temp^3+8.113e-7*temp^2-1.604e-4*temp+0.04444;
        viscosity = 2.364e-16*temp^4-4.731e-13*temp^3+3.247e-10*temp^2-5.675e-8*temp+2.681e-5;
        density = 5.589e-9*temp^4-1.149e-5*temp^3+8.529e-3*temp^2-2.830*temp+458;
    else
        conductivity = 4.260e-8*temp^3+5.432e-6*temp^2-2.195e-3*temp+0.1584;
        viscosity = -3.747e-10*temp^3+9.747e-8*temp^2-8.502e-6*temp+2.743e-4;
        density = -2.115e-3*temp^3+0.6298*temp^2-64.56*temp+2527;
    end
end