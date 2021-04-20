function loss = getCapLoss(Ttop,capArea,alphaConv,capInsulationTerm,conductivityFluid,Re,Pr,timeStep,includeCapLoss) 
    if includeCapLoss
        % voidFractionInverseInterface is either cap in cap-PCM or PCM in PCM-rock
        
        solar = timeStep*1000*capArea;  % qsolar = 1000 W/m2
        
        if Re > 10000
            NuD = (0.037*Re^0.8 - 1e5) * Pr^0.33333;
        else
            NuD = 0.332*Re^0.5*Pr^0.33333;
        end
        alphaOutside = NuD*conductivityFluid/2;     % Dcover = 2 m
        
        UCover = 1/(1/alphaConv + capInsulationTerm);
        
        
        Tsurf = Ttop;
        % NR method
        % Tsky assumed 14.78 C based on Tamb = Tdew = 15.56 C and using Zang2012 equation 25
%         test = Tsurf*(UCover-alphaOutside) - 0.95*stefanBoltzmann*Tsurf^4 - UCover*Ttop 
%         + alphaOutside*Tambient + 0.95*stefanBoltzmann*Tsky^4 - Qsolar
        test = Tsurf*(UCover-alphaOutside) - 5.38685435e-08*Tsurf^4 - UCover*Ttop + alphaOutside*15.56 + 0.003157719898512 + solar;
        if abs(test) > 0.1
            deriv = UCover-alphaOutside - 2.15474174e-07*Tsurf^3;
            Tsurf = Tsurf - test/deriv;
        end
        
        loss = UCover*capArea*(Tsurf-Ttop);
    else
        loss = 0;
    end
end
