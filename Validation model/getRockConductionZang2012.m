% function cond = getRockConduction(includeConduction,keffprev,keff,tsprev,ts,tsnext,timeStep,rowSurfaceAreaRocks,voidFractionRockInverse,densityRocks,rowLengthRock,rowVolumeRock)
function cond = getRockConductionZang2012(includeConduction,tsprev,ts,tsnext,keffprev,keffcurr,timeStep,tankAreaPrev,tankAreaCurr,rowLength)
if includeConduction
    cond = timeStep*(keffcurr*tankAreaCurr*(tsnext-ts)-keffprev*tankAreaPrev*(ts-tsprev))/rowLength;

%     cond = timeStep*rowSurfaceAreaRocks/voidFractionRockInverse/densityRocks/rowLengthRock/rowVolumeRock*(keffprev*tsprev+keff*tsnext-(keff+keffprev)*ts);
else
    cond=0;
end