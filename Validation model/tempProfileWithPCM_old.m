function [temp] = tempProfileWithPCM(t,timeStep,tempIncreaseAmount,prevTemp)
    sec = t*timeStep;
    if sec > 7200
        temp = max(prevTemp,tempIncreaseAmount-0.00000000000076251*sec^4+0.000000026076*sec^3-0.00032977*sec^2+1.8379*sec-3221.2);
    elseif sec > 3600
        temp =tempIncreaseAmount+0.00000000000029878*sec^4-0.000000006019*sec^3+0.000030164*sec^2+0.06656*sec;
    else
        temp = tempIncreaseAmount+max(25,-0.0000000065729*sec^3+0.000056327*sec^2-0.013426*sec+25);
    end
end