function [temp] = tempProfileTrahan(t,timeStep,~,~,~,~,includeDischarge)
    sec = t*timeStep;
    
% need to update
if includeDischarge
    if sec > 6175
        temp = 192;
    elseif sec > 5600
        temp = -5.31905E-12*sec.^5 + 1.58588E-07*sec.^4 - 1.89112E-03*sec.^3 + 1.12743E+01*sec.^2 - 3.36031E+04*sec + 4.00574E+07;
    elseif sec > 5536
        temp = 1.96205E-06*sec.^4 - 4.37788E-02*sec.^3 + 3.66296E+02*sec.^2 - 1.36208E+06*sec + 1.89926E+09;
    elseif sec > 5400
        temp = -9.75008E-10*sec.^6 + 3.20521E-05*sec.^5 - 4.39016E-01*sec.^4 + 3.20695E+03*sec.^3 - 1.31769E+07*sec.^2 + 2.88751E+10*sec - 2.63638E+13;
    elseif sec > 4290
        temp = 2.65021E-06*sec.^2 - 2.68113E-02*sec + 5.74995E+02;
%         disp('return from 4290')
    elseif sec > 3327.7
%         disp('1')
        temp = -2.33915E-08*sec.^3 + 2.79829E-04*sec.^2 - 1.11539E+00*sec + 1.99027E+03;
    elseif sec > 436
        temp = -2.26464E-12*sec.^4 + 2.21250E-08*sec.^3 - 8.14822E-05*sec.^2 + 1.57211E-01*sec + 3.57176E+02;
    else
%         disp('return from bottom')
        temp = -5.13902E-13*sec.^6 + 9.08409E-10*sec.^5 - 6.45202E-07*sec.^4 + 2.36736E-04*sec.^3 - 4.83670E-02*sec.^2 + 5.68305E+00*sec + 3.95292E+01;
    end
else
    if sec > 5400
        temp = 500;
    elseif sec > 4290
        temp = 2.65021E-06*sec.^2 - 2.68113E-02*sec + 5.74995E+02;
%         disp('return from 4290')
    elseif sec > 3327.7
%         disp('1')
        temp = -2.33915E-08*sec.^3 + 2.79829E-04*sec.^2 - 1.11539E+00*sec + 1.99027E+03;
    elseif sec > 436
        temp = -2.26464E-12*sec.^4 + 2.21250E-08*sec.^3 - 8.14822E-05*sec.^2 + 1.57211E-01*sec + 3.57176E+02;
    else
%         disp('return from bottom')
        temp = -5.13902E-13*sec.^6 + 9.08409E-10*sec.^5 - 6.45202E-07*sec.^4 + 2.36736E-04*sec.^3 - 4.83670E-02*sec.^2 + 5.68305E+00*sec + 3.95292E+01;
    end
end
    if isnan(temp)
        MEx = MException('tempProfileTrahan:TempIsNotANumber',...
            'Fluid energy NaN at t= %g with timestep=%g and t*timeStep=\n',t,timeStep,t*timeStep);
        throw(MEx)
    end
end