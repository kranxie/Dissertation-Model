function press = getPressureDrop(rowLength,diameterRock, voidFractionRock, voidFractionInverse,dTemp, temp, viscosityFluid,densityFluid,g)
% press = 27465.82*viscosityF/densityfluid*g+17.0898*g^2/densityfluid;
press = rowLength * g^2 ./ densityFluid / diameterRock.*(217*voidFractionInverse.^2./voidFractionRock.^3/0.6^2.*viscosityFluid/g/diameterRock ...
    + 1.83.*voidFractionInverse./voidFractionRock.^3/0.6) + densityFluid*9.81.*rowLength.*dTemp./temp;

% press = rowLength
end