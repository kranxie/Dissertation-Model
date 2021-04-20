function cp = getCpDryAirAtmo(temp)
cp = -3.448127E-07*temp.^3 + 4.742529E-04*temp.^2 + 2.535406E-02*temp + 1.004398E+03;
end