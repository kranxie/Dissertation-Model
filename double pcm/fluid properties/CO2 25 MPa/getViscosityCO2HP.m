function viscosity = getViscosityCO2HP(temp)% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in Pa-s
% R^2 = 0.99997. Max error 0.10% for temp 320-800 C
viscosity = -0.000000000000023395*temp.^3 + 0.000000000041902*temp.^2 + 0.0000000022433*temp + 0.000028131;
end

