function viscosity = getViscosityDryAirAtmo(temp)
% viscosity = -1.248e-11*temp..^2+4.375e-8*temp+1.746e-5;
viscosity = -6.7498152E-24*temp.^6 + 2.3857650E-20*temp.^5 - 3.8746019E-17*temp.^4 + 4.1307105E-14*temp.^3 - 3.6335761E-11*temp.^2 + 5.0156235E-08*temp + 1.7258181E-05;
end