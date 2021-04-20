function loss = getQWallPCM(Tpcm,Twall,Ac,fcontwall,rin,tenc,VFPCMInverse,densityPCM,rowVolumePCM,timeStep,includeWallLoss)
if includeWallLoss
    if Tpcm < 573
        kenc = 17.7;
    elseif Tpcm < 577
        kenc = 22;
    else
        kenc = 23.1;
    end
    
    loss = timeStep*kenc*Ac*fcontwall/(rin-tenc)/log(rin/(rin-tenc))*(Twall-Tpcm)/VFPCMInverse/densityPCM/rowVolumePCM;
else
    loss = 0;
end
end

