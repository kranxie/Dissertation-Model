function specificheat = getSpecificHeatDryAir(temp)
specificheat =-0.0000019796792*temp.^3 + 0.003099531*temp.^2 - 1.4184422*temp + 1310.608;
end