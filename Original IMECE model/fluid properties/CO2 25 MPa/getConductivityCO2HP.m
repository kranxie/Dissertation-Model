function conductivity = getConductivityCO2HP(temp)
% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in W/mK
% R^2=0.99995 (low T), =0.99978(high T)
    if temp > 150
        conductivity = 0.0000000000013325*temp^4 - 0.0000000027511*temp^3 + 0.0000020201*temp^2 - 0.00056114*temp + 0.098069;
    else
        conductivity = 0.000003328*temp^2 - 0.0011091*temp + 0.14244;
    end
end
