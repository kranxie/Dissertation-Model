function [temp] = tempProfileTrahanPCM(t,timeStep,~,~,~,~,~)
    sec = t*timeStep;
    
% need to update
    temp = 326;
    if isnan(temp)
        MEx = MException('tempProfileTrahan:TempIsNotANumber',...
            'Fluid energy NaN at t= %g with timestep=%g and t*timeStep=\n',t,timeStep,t*timeStep);
        throw(MEx)
    end
end