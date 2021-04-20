function rad = getInterfaceRadiation(tempCurrentRow,tempInterfacedRow,tankArea,voidFractionInverseInterface,emissitivtyCurrentRow,emissivityInterfacedRow,timeStep,~,~,~,includeRadiation) 
    if includeRadiation
        % voidFractionInverseInterface is either cap in cap-PCM or PCM in PCM-rock
        rad = timeStep*5.67e-8*tankArea*voidFractionInverseInterface/(1/emissitivtyCurrentRow+1/emissivityInterfacedRow-1)*((tempInterfacedRow+273.15)^4-(tempCurrentRow+273.15)^4);
    else
        rad = 0;
    end
end
