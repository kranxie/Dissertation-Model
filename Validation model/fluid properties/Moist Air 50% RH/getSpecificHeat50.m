function specificheat = getSpecificHeat50(temp)
specificheat = (-3.425e-7*temp.^3 + 4.940e-4*temp.^2 + 0.01043e-2*temp + 1014);
end