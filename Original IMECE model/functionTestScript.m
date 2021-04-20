clc

fluidStepConvection = 0;

tempStorage = 650;
tempFluid = 600;
hvEffRocks = 2500;
densityFluid = 0.3715;
timeStep = 40;
voidFractionRocks = 0.4;

getFluidConvectionHandle = @(hveff,fluiddensity,ts,tf) timeStep*hveff/voidFractionRocks/fluiddensity*(ts-tf);

disp('Anonymous function:')
tic
for i=1:1:79968000
    fluidStepConvection = getFluidConvectionHandle(hvEffRocks,densityFluid,tempStorage,tempFluid);
end
toc

disp('Separate function:')
tic
for i=1:1:79968000
    fluidStepConvection = getFluidConvectionFunc(hvEffRocks,densityFluid,tempStorage,tempFluid,timeStep,voidFractionRocks);
end
toc