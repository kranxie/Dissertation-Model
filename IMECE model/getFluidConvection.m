function conv = getFluidConvection(hveff,fluiddensity,ts,tf,timeStep,voidFractionRocks)
conv = timeStep*hveff/voidFractionRocks/fluiddensity*(ts-tf);

end