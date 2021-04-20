function enthalpy = getEnthalpyDryAir(temp)

%    enthalpy = (0.00000000025104/5*temp..^5-0.0000007284/4*temp..^4+0.00066686/3*temp..^3-0.010204/2*temp..^2+1006*temp+23294.44473);
enthalpy = -0.0000004949198*temp.^4 + 0.001033177*temp.^3 - 0.7092211*temp.^2 + 1310.608*temp + 226503.8;
end