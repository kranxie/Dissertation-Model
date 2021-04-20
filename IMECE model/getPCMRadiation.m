function rad = getPCMRadiation(tempPrevRow,tempCurrentRow,tempNextRow,tankArea,densityPCM,rowVolumePCM,emissivityEncaps,timeStep,includeRadiation) 
    if includeRadiation
        rad = timeStep*5.67e-8*tankArea/(2/emissivityEncaps-1)*((tempNextRow+273.15)^4+(tempPrevRow+273.15)^4-2*(tempCurrentRow+273.15)^4)/densityPCM/rowVolumePCM;
        % unit: watts
    else
        rad = 0;
    end
end
