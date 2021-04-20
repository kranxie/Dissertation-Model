function specificheat = getCpCO2HP(temp)
% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in J/kgK
% Derivative of enthalpy
specificheat = 0.0053523*temp^2 -5.4414*temp+ 2503.6;
% 0.0053523*temp^2 - 5.4414*temp^2 + 2503.6;
end