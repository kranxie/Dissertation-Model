function loss = getFluidHeatLoss(includeHeatLoss,i,alphaRad,tf,tinf,lossInsulationTerms,fluiddensity,kFluid,Re,Pr,voidFractionRocks,rockAverageDiameter,rowVolumeRock,rowWallSurfaceAreaRock,timeStep)
if includeHeatLoss
    if i < 28   % Adjust these values for actual conditions
                lossInsulationTerm = lossInsulationTerms(1);
        elseif i == 28
            lossInsulationTerm = (lossInsulationTerms(1)+lossInsulationTerms(2))/2;
        elseif i < 49
                lossInsulationTerm = lossInsulationTerms(2);
        else
                lossInsulationTerm = lossInsulationTerms(3);
        end
        % Calculate uInside
    %     Re = G*rockAverageDiameter/viscosityFluid(i);
    %     Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
        alphaConv = kFluid/rockAverageDiameter*(3.22*(Re*Pr)^(0.333)+0.117*Re^0.8*Pr^0.4);
    %     alphaRad = 0;
        uInside = alphaConv + alphaRad;

        uWall = 1/(1/uInside + lossInsulationTerm);
        loss = timeStep*uWall*rowWallSurfaceAreaRock*(tf-tinf)/voidFractionRocks/fluiddensity/rowVolumeRock;
else
    loss = 0;
end