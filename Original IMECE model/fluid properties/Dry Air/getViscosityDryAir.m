function viscosity = getViscosityDryAir(temp)
% viscosity = -1.248e-11*temp..^2+4.375e-8*temp+1.746e-5;
viscosity =4.2455806E-22*temp.^6 - 1.196575E-18*temp.^5 + 1.3509552E-15*temp.^4 - 7.8319866E-13*temp.^3 + 0.00000000024422398*temp.^2 - 0.0000000093727191*temp + 0.000026240582;
end