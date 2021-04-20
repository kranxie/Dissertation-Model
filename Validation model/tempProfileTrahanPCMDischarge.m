function [temp] = tempProfileTrahanPCMDischarge(t,timeStep,chargeTime)
    sec = (t-chargeTime)*timeStep;
    
% need to update
if sec > 1500
    temp = 286;
elseif sec > 893
    temp = -2.95378E-15*sec.^5 + 0.0000000000324813*sec.^4 - 0.000000139235*sec.^3 + 0.00029049*sec.^2 - 0.295397*sec + 404.018;
elseif sec > 115
    temp = -4.09184E-05*sec.^2 + 1.40631E-02*sec + 3.14781E+02;
else
    temp = 4.34869E-04*sec.^2 - 1.36970E-01*sec + 3.25968E+02;
end
    
%     disp([sec temp])
    
    
    if isnan(temp)
        MEx = MException('tempProfileTrahan:TempIsNotANumber',...
            'Fluid energy NaN at t= %g with timestep=%g and t*timeStep=\n',t,timeStep,t*timeStep);
        throw(MEx)
    end
end