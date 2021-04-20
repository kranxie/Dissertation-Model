function cp = getCpCO2HP(temp)
% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in J/kgK
% Derivative of enthalpy
cp = 0.00040137*temp.^2 - 0.3572*temp + 1329.1;
end