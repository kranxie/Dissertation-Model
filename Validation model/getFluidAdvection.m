function advection = getFluidAdvection(mdot,~,hin,hout,timeStep,~,~) 
    advection = timeStep*mdot*(hin-hout);
end