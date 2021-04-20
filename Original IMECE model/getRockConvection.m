function conv = getRockConvection(hveff,tf,ts,timeStep,fluidVolumeFraction,densityRocks) 
conv = timeStep*hveff/fluidVolumeFraction/densityRocks*(tf-ts);
end
