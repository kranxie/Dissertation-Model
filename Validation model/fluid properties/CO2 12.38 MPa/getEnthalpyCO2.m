function enthalpy = getEnthalpyCO2(temp)
% Properties valid for pressures of approximately 12.38 MPa
% Temperatures in degC
    enthalpy = -0.4590*temp.^2+1666*temp+295300;
end