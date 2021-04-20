function enthalpyDerivative = getEnergyDerivativeDryAir(temp)
% enthalpyDerivative = (0.00000000025104*temp..^4-0.0000007284*temp..^3+0.00066686*temp..^2-0.010204*temp+1006);
enthalpyDerivative =-0.0000019796792*temp.^3 + 0.003099531*temp.^2 - 1.4184422*temp + 1310.608;
end