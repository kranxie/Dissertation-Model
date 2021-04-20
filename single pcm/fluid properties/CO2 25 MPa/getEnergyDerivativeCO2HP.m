function enthalpyDerivative = getEnergyDerivativeCO2HP(temp)
% **** DELETE THIS WHEN PROPERTY IS UPDATED FOR 25 MPA ***


% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in J/kgK
enthalpyDerivative = -0.9180*temp+1666;
end