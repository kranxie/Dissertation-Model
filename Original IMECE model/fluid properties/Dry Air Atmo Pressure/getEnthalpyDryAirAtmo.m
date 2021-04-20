function enthalpy = getEnthalpyDryAirAtmo(temp)

%    enthalpy = (0.00000000025104/5*temp..^5-0.0000007284/4*temp..^4+0.00066686/3*temp..^3-0.010204/2*temp..^2+1006*temp+23294.44473);
enthalpy = -8.620318E-08*temp.^4 + 1.580843E-04*temp.^3 + 1.267703E-02*temp.^2 + 1.004398E+03*temp + 2.733495E+05;
end