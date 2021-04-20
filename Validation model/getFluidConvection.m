function conv = getFluidConvection(hveff,ts,tf,timeStep,rowVolumeStorage,~,~)
conv = timeStep*hveff*rowVolumeStorage*(ts-tf);

end