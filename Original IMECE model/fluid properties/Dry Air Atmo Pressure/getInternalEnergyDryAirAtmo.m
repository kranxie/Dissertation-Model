function energy = getInternalEnergyDryAirAtmo(temp)

%    enthalpy = (0.00000000025104/5*temp..^5-0.0000007284/4*temp..^4+0.00066686/3*temp..^3-0.010204/2*temp..^2+1006*temp+23294.44473);
energy =-0.00000008532913*temp.^4+0.000156176*temp.^3+0.01428662*temp.^2+716.5967*temp+194965.1;
end