function conductivity = getConductivity50(temp)
conductivity = -1.387e-8*temp.^2 + 6.819e-5*temp + 2.438e-2;
end