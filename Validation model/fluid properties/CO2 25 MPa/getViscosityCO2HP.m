function viscosity = getViscosityCO2HP(temp)% Properties valid for pressures of approximately 25 MPa
% Temperatures in degC, output in Pa-s
% R^2 = 0.99957 (low T), =0.99973 (high T). Max error 1.12%
    if temp > 150
        viscosity = 7.7411E-16*temp^4 - 0.0000000000016093*temp^3 + 0.0000000011933*temp^2 - 0.00000034619*temp + 0.000064945;
    else
        viscosity = 0.000000003939*temp^2 - 0.0000012074*temp + 0.00012765;
    end
end

