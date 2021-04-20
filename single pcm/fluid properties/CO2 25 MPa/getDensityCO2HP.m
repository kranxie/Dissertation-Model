function density = getDensityCO2HP(temp)
% Properties valid for pressures of appro*temp.^imately 25 MPa
% Temperatures in degC, output in kg/m3
% R^2 = 0.99999. Max error = 0.08% (vs refprop from 320-800 C).
    density = 0.0000000011588*temp.^4 - 0.0000032168*temp.^3 + 0.0035117*temp.^2 - 1.9174*temp + 575.18;
end
