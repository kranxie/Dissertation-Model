function density = getDensityDryAirAtmo(temp)
% density = 101.3e3./287.1./(temp+273.15);
% density = 25e6./287.1./(temp+273.15);

density = -7.6184155E-15*temp.^5 + 1.9832525E-11*temp.^4 - 2.0929446E-08*temp.^3 + 1.1980056E-05*temp.^2 - 4.4247409E-03*temp + 1.2878778E+00;

end