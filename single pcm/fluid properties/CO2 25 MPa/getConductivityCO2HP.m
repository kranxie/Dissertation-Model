function conductivity = getConductivityCO2HP(temp)
% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in W/mK
% R^2=0.99999; max error = 0.09% from 320-800 C
    conductivity = -0.000000000039554*temp.^3 + 0.000000071063*temp.^2 + 0.000019831*temp + 0.038129;
end
