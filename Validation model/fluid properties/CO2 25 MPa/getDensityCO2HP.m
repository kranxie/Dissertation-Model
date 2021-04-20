function density = getDensityCO2HP(temp)
% Properties valid for pressures of appro*temp.^imately 25 MPa
% Temperatures in degC, output in kg/m3
% R^2 = 0.99918 (low T), =0.99996 (high T). Max error = 1.12% (vs refprop from 50-800 C).
    if temp > 125
        density = 0.000000010027*temp.^4 - 0.000021306*temp.^3 + 0.016564*temp.^2 - 5.8341*temp + 984.15;
    else
        density = 0.013752*temp.^2 - 7.0248*temp + 1156.5 ;
    end
end
