function loss = getQWallPCM(Tpcm,Twall,areaTerms,~,~,~,timeStep,includeWallLoss)
if includeWallLoss
    if Tpcm < 573
        kenc = 17.7;
    elseif Tpcm < 577
        kenc = 22;
    else
        kenc = 23.1;
    end
    loss = timeStep*kenc*areaTerms*(Twall-Tpcm);
%     areaTerms = Ac*fcontwall/(rin-tenc)/log(rin/(rin-tenc))
else
    loss = 0;
end
end

