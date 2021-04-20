% function [temp] = tempProfileNoPCM(~,~,tempIncreaseAmount,~)
%     temp = 500+tempIncreaseAmount;
% end

function [temp] = tempProfileNoPCM(t,timeStep,tempIncreaseAmount,prevTemp,includePCM,tempProfileModifier)
    sec = t*timeStep;
    if includePCM
        if sec > 7200
            temp = tempIncreaseAmount+tempProfileModifier*(-0.00000000000076251*sec^4+0.000000026076*sec^3-0.00032977*sec^2+1.8379*sec-3221.2);
%             temp = max(prevTemp,tempIncreaseAmount+tempProfileModifier*(-0.00000000000076251*sec^4+0.000000026076*sec^3-0.00032977*sec^2+1.8379*sec-3221.2));
        elseif sec > 3600
            temp =tempIncreaseAmount+tempProfileModifier*(0.00000000000029878*sec^4-0.000000006019*sec^3+0.000030164*sec^2+0.06656*sec);
        else
            temp = tempIncreaseAmount+tempProfileModifier*(max(25,-0.0000000065729*sec^3+0.000056327*sec^2-0.013426*sec+25));
        end
    else
        if sec > 5
            temp = max(prevTemp,tempIncreaseAmount+tempProfileModifier*(-3.9993E-17*sec^5+0.0000000000013797*sec^4-0.000000016985*sec^3+0.000080394*sec^2-0.023642*sec+33));
        else
            temp =max(25,tempIncreaseAmount+tempProfileModifier*(-3.9993E-17*sec^5+0.0000000000013797*sec^4-0.000000016985*sec^3+0.000080394*sec^2-0.023642*sec+33));
        end
    end
end