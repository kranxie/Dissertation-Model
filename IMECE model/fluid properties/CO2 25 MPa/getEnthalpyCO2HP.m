function enthalpy = getEnthalpyCO2HP(temp)
% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in J/kg
% R^2 = 0.99971. Max error = 1.56% (vs refprop from 50-800 C).
    enthalpy = 0.0017841*temp.^3 - 2.7207*temp.^2 + 2503.6*temp + 177350;
end