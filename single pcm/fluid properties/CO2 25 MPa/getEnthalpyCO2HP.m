function enthalpy = getEnthalpyCO2HP(temp)
% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in J/kgK
% R^2 = 1. Max error = 0.03% (vs refprop from 320-800 C).
    enthalpy = 0.00013379*temp.^3 - 0.1786*temp.^2 + 1329.1*temp + 332710;

end