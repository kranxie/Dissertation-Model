function cond = getPCMConductionZanganeh(includeConduction,tsprev,ts,tsnext,~,~,timeStep,fcont,tankArea,~,~,rowLength)
if includeConduction
    if ts < 573
        kcurr = 17.7;
    elseif ts < 577
        kcurr = 22;
    else
        kcurr = 23.1;
    end
    if tsprev < 573
        kprev = 17.7;
    elseif tsprev < 577
        kprev = 22;
    else
        kprev = 23.1;
    end

    cond = timeStep*tankArea*fcont*(kcurr*(tsnext-ts)-kprev*(ts-tsprev))/rowLength;
%     cond =
%     timeStep*tankArea*fcont*(kprev*tsprev+kcurr*tsnext-(kcurr+kprev)*ts)/voidFractionInverse/densityPCMEff/rowVolumePCM;1
else
    cond = 0;
end
end