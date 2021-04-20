function conductivity = getConductivityCO2(temp)
% Properties valid for pressures of approximately 12.38 MPa
% Temperatures in degC
    if temp > 100
        conductivity = 5.788e-13*temp.^4-1.163e-9*temp.^3+8.113e-7*temp.^2-1.604e-4*temp+0.04444;
    else
        conductivity = 4.260e-8*temp.^3+5.432e-6*temp.^2-2.195e-3*temp+0.1584;
    end
end