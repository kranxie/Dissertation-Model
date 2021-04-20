function [tempMatrixStorage,tempMatrixFluid] = calculateTempsNoRocksrefprop(intEnergyMatrixStorage,intEnergyMatrixFluid,tempMatrixStorage,tempMatrixFluid,lowCutoff,highCutoff, TmeltStart, TmeltEnd, cEffSolid, cEffPhaseChange, cEffLiquid, ~, ~, ~,~,fluid,pressure)
% intEnergyMatrix,includePCM,fluidPressure,RH
    numRows = length(intEnergyMatrixFluid);
    
    for i = 1:1:numRows
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
end