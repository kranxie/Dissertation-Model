% function cond = getRockConduction(includeConduction,keffprev,keff,tsprev,ts,tsnext,timeStep,rowSurfaceAreaRocks,voidFractionRockInverse,densityRocks,rowLengthRock,rowVolumeRock)
function cond = getRockConduction(includeConduction,tsprev,ts,tsnext,keffprev,keff,timeStep,tankArea,rowLength)
if includeConduction
    cond = timeStep*tankArea*(keff*(tsnext-ts)-keffprev*(ts-tsprev))/rowLength;

%     cond = timeStep*rowSurfaceAreaRocks/voidFractionRockInverse/densityRocks/rowLengthRock/rowVolumeRock*(keffprev*tsprev+keff*tsnext-(keff+keffprev)*ts);
else
    cond=0;
end