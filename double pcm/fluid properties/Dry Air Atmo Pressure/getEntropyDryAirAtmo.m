function enthalpy = getEntropyDryAirAtmo(temp)

enthalpy = -0.0012*temp.^2 + 2.6762*temp + 6823.3;
end