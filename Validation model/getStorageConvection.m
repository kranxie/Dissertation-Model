function conv = getStorageConvection(hveff,tf,ts,timeStep,rowVolumeStorage,~,~) 
conv = timeStep*hveff*rowVolumeStorage*(tf-ts);
end
