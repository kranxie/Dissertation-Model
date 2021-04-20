function rad = getInterfaceRadiation(tempPCM,tempOtherMaterial,tankArea,voidFractionInverse,emissitivtyPCM,emissivityOtherMaterial,timeStep,voidFractionInverseMaterial,densityMaterial,rowVolumeMaterial,includeRadiation) 
    if includeRadiation
        rad = timeStep*5.67e-8*tankArea*voidFractionInverse/(1/emissitivtyPCM+1/emissivityOtherMaterial-1)*((tempOtherMaterial+273.15)^4-(tempPCM+273.15)^4)/voidFractionInverseMaterial/densityMaterial/rowVolumeMaterial;
        % unit: watts
    else
        rad = 0;
    end
end
