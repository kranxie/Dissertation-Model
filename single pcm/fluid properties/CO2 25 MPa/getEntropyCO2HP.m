function entropy = getEntropyCO2HP(temp)
% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in J/kgK
% R^2 = 1. Max error = 0.03% (vs refprop from 320-800 C).
    entropy = 9.0786E-07*temp.^3 - 2.4191E-03*temp.^2 + 3.3603E+00*temp + 1.4263E+03;
end