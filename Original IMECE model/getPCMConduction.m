function cond = getPCMConduction(includeConduction,tsprev,ts,tsnext,timeStep,tankArea,fcont,fluidVolumeFraction,densityPCMEff,rowVolumePCM)
if includeConduction
    if ts < 573
        kencaps = 17.7;
    elseif ts < 577
        kencaps = 22;
    else
        kencaps = 23.1;
    end

    cond = timeStep*tankArea*fcont*kencaps*(tsprev+tsnext-2*ts)/fluidVolumeFraction/densityPCMEff/rowVolumePCM;
else
    cond = 0;
end
end