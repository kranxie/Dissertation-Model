function density = getDensityDryAir(temp)
% density = 101.3e3./287.1./(temp+273.15);
% density = 25e6./287.1./(temp+273.15);

density =6.4418097E-15*temp.^6 - 0.000000000018641417*temp.^5 + 0.000000021959079*temp.^4 - 0.000013767197*temp.^3 + 0.0051668217*temp.^2 - 1.3267948*temp + 304.39537;

end