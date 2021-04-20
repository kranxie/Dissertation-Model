function cv = getCvDryAirAtmo(temp)
cv = -3.413165E-07 * temp.^3 + 4.685280E-04 * temp.^2 + 2.857324E-02 * temp + 7.165967E+02;
end