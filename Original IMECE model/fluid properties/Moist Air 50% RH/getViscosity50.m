function viscosity = getViscosity50(temp)
viscosity = -1.189e-11*temp.^2 + 4.411e-8*temp + 1.755e-5;
end