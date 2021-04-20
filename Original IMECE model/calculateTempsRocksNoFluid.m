function [tempMatrixStorage] = calculateTempsRocksNoFluid(intEnergyMatrixStorage,tempMatrixStorage)
% intEnergyMatrix,includePCM,fluidPressure,RH
    numRows = length(intEnergyMatrixStorage);
    
    % Rock temps
    tempMatrixStorage(1:numRows) = (-747.0995+sqrt(747.0995^2+4*0.2838*intEnergyMatrixStorage(1:numRows)))/2/0.2838;
    
    if any(imag(tempMatrixStorage(:)))
            MEx = MException('calculateTemps:StorageTempImaginary',...
                'Calculated temp imaginary');
        throw(MEx)
    end
end