function viscosity = getViscosityCO2(temp)
% Properties valid for pressures of approximately 12.38 MPa
% Temperatures in degC
    if temp > 100
        viscosity = 2.364e-16*temp.^4-4.731e-13*temp.^3+3.247e-10*temp.^2-5.675e-8*temp+2.681e-5;
    else
        viscosity = -3.747e-10*temp.^3+9.747e-8*temp.^2-8.502e-6*temp+2.743e-4;
    end
end