function cond = getPCMConductionBeasley(includeConduction,tsprev,ts,tsnext,~,~,timeStep,~,~,voidFractionInverse,densityPCMEff,rowLength)
% function cond = getPCMConductionBeasley(includeConduction,tsprev,ts,tsnext,timeStep,tankArea,fcont,voidFractionInverse,densityPCMEff,rowVolumePCM,rowLength)
if includeConduction
%     kcurr = 0.24;
%     kprev = 0.24;

%     cond = timeStep*tankArea*fcont*(kprev*tsprev+kcurr*tsnext-(kcurr+kprev)*ts)/voidFractionInverse/densityPCMEff/rowVolumePCM;
    cond = timeStep/voidFractionInverse/densityPCMEff/rowLength^2*0.24*(tsprev+tsnext-2*ts);
else
    cond = 0;
end
end