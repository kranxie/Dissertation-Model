function density = getDensityCO2(temp)
% Properties valid for pressures of appro*temp.^imately 12.38 MPa
% Temperatures in degC
    if temp > 100
%         density = 5.589e-9*temp.^4-1.149e-5*temp.^3+8.529e-3*temp.^2-2.830*temp+458;
        density = 5.589E-09*temp.^4 - 1.149E-05*temp.^3 + 8.529E-03*temp.^2 - 2.830E+00*temp + 4.580E+02;
    else
        density = -2.115e-3*temp.^3+0.6298*temp.^2-64.56*temp+2527;
    end
%     density = 1.^1;
end