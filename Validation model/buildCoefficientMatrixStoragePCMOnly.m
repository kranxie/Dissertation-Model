function [coefficientMatrix] = buildCoefficientMatrixFluidPCMOnly(baseMatrix,tempFluid,tempStorage,rockProperties,fluidProperties,hvEff)
% (kRocks,kFluid,fluidTemps,rockTemps,voidFractionRocks,d,numRows,epsRocks)
rows = length(baseMatrix);
for i = 2:1:rows-1
% fluid
% advection -- can't really vectorize this properly while using variable
% properties. Will probably need to put this on the RHS
%     coefficientMatrix(i,i) = baseMatrix + 0;

% convection
    coefficientMatrix(i,i) = baseMatrix - hvEff*rowVolume;
    
% heat loss
    coefficientMatrix(i,i) = baseMatrix - hWall(i)*rowSurfaceArea;

% storage
    coefficientMatrix(i,i) = baseMatrix - 
% conduction / radiation
% convection

end
% i = 1

% i = rows

end