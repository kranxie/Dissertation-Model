function cond = getRockConduction(includeConduction,keffprev,keff,tsprev,ts,tsnext,timeStep,rowSurfaceAreaRocks,voidFractionRockInverse,densityRocks,rowLengthRock,rowVolumeRock)
if includeConduction
    cond = timeStep*rowSurfaceAreaRocks/voidFractionRockInverse/densityRocks/rowLengthRock/rowVolumeRock*(keffprev*tsprev+keff*tsnext-(keff+keffprev)*ts);
else
    cond=0;
end