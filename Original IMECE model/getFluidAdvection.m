function advection = getFluidAdvection(mdot,fluiddensity,hin,hout,timeStep,rowVolumeRock,voidFractionRocks) 
advection = timeStep*mdot/fluiddensity/rowVolumeRock/voidFractionRocks*(hin-hout);
end