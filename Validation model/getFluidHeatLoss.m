function loss = getFluidHeatLoss(includeHeatLoss,alphaConv,alphaRad,tf,tinf,lossInsulationTerm,rowWallSurfaceArea,timeStep)
    if includeHeatLoss
        % Calculate uInside
        %     Re = G*rockAverageDiameter/viscosityFluid(i);
        %     Pr = specificHeatFluid(i)*viscosityFluid(i)/kFluid(i);
        
        % alphaConv has been verified
%         alphaConv = kFluid/rockAverageDiameter*(3.22*(Re*Pr)^(0.333333)+0.117*Re^0.8*Pr^0.4);
        
        % alphaRad = heffWall(i)
        % heffWall comes from calculatekeff
        
        % alphaConv has been verified, assuming alphaRad is correct
        uInside = alphaConv + alphaRad;

        % uWall has been verified
        uWall = 1/(1/uInside + lossInsulationTerm);

%         loss = Qwallrocks 
        % rowWallSurfaceAreaRock in 2012 is Cn = row dodecagon surface
        % area, because 2012 uses a dodecagon-shaped cross section instead of a true
        % cone (with circular cross sections). Perhaps I should adjust my
        % geometry accordingly, but that's not in here.
        loss = timeStep*uWall*rowWallSurfaceArea*(tinf-tf);
    else
        loss = 0;
    end
end