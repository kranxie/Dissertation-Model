function viscosity = getViscosityTrahan(temp)
% viscosity = -1.248e-11*temp..^2+4.375e-8*temp+1.746e-5;
viscosity = (6.10504e-10)*(temp^3)-(2.13036e-6)*(temp^2)+(4.71398e-3)*temp+1.67555e-5;
% kg/ms = Pa*s, theoretically, but it seems there might be a typo in the
% dissertation somewhere here. I'm going to use my REFPROP-derived
% equations instead for the time being.

end