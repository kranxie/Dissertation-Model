function [tempMatrixStorage,tempMatrixFluid] = calculateTempsWithPCMrefprop(intEnergyMatrixStorage,intEnergyMatrixFluid,tempMatrixStorage,tempMatrixFluid,lowCutoff,highCutoff, TmeltStart, TmeltEnd, cEffSolid, cEffPhaseChange, cEffLiquid, rowCutoffPCM, rowStartRock, ~,~,fluid,pressure)
% intEnergyMatrix,includePCM,fluidPressure,RH
    numRows = length(intEnergyMatrixStorage);
    
    % Rock temps
    tempMatrixStorage(rowStartRock:numRows) = (-747.0995+sqrt(747.0995^2+4*0.2838*intEnergyMatrixStorage(rowStartRock:numRows)))/2/0.2838;
    
    % PCM
    for i = 1:1:rowCutoffPCM
        if intEnergyMatrixStorage(i) > highCutoff       % Liquid
            tempMatrixStorage(i) = TmeltEnd +(intEnergyMatrixStorage(i)-highCutoff)/cEffLiquid;
        elseif intEnergyMatrixStorage(i) > lowCutoff      % phase change
            tempMatrixStorage(i) = TmeltStart+(intEnergyMatrixStorage(i)-lowCutoff)/cEffPhaseChange;
        else                                    % solid
            tempMatrixStorage(i) = (intEnergyMatrixStorage(i))/cEffSolid;
        end
    end  
    
    if any(imag(tempMatrixStorage(:)))
            MEx = MException('calculateTemps:StorageTempImaginary',...
                'Calculated temp imaginary');
        throw(MEx)
    end
    temp = tempMatrixFluid(1);
    % air temps
    for i = 1:1:numRows
        
        if ~isreal(intEnergyMatrixFluid(i))
            disp(intEnergyMatrixFluid(i))
            MEx = MException('calculateTemps:FluidEnergyImaginary',...
                'Fluid energy imaginary (%g J/kg) at row %g',intEnergyMatrixFluid(i),i);
            throw(MEx)
        end
        
        temp = refpropm('t','p',pressure,'h',intEnergyMatrixFluid(i),fluid)-273;
        if temp < 24.95
            MEx = MException('calculateTemps:FluidTempOutOfRange',...
                'Calculated Temp out of bounds (%g C)\n   at row %g for energy %g J/kg',temp,i,intEnergyMatrixFluid(i));
            throw(MEx)
        end
        if ~isreal(temp)
            disp(temp)
            MEx = MException('calculateTemps:FluidTempImaginary',...
                'Calculated temp imaginary (%g C) at row %g for energy %g',temp,i,intEnergyMatrixFluid(i));
            throw(MEx)
        end
        tempMatrixFluid(i) = temp;
    end
%     disp(tempMatrixFluid)
end