function rad = getPCMRadiation(tempPrevRow,tempCurrentRow,tempNextRow,tankArea,~,~,emissivityEncaps,timeStep,voidFractionInverse,includeRadiation) 
    if includeRadiation
        rad = timeStep*5.67e-8*tankArea*voidFractionInverse/(2/emissivityEncaps-1)*((tempNextRow+273.15)^4+(tempPrevRow+273.15)^4-2*(tempCurrentRow+273.15)^4);
        % * voidFractionInverse / voidFractionInverse removed for
        % redundancy sake
        % unit: watts
    else
        rad = 0;
    end
end
