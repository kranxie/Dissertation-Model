% function [temp] = tempProfileNoPCM(~,~,tempIncreaseAmount,~)
%     temp = 500+tempIncreaseAmount;
% end

function [temp] = tempProfileZanganeh(t,timeStep,tempIncreaseAmount,prevTemp,includePCM,tempProfileModifier,~)
    sec = t*timeStep;
    if includePCM
        if sec > 8414
            temp = max(tempIncreaseAmount+tempProfileModifier*(4.760629E-19*sec^6 - 2.891864E-14*sec^5 + 0.0000000007292832*sec^4 - 0.000009773364*sec^3 + 0.07340846*sec^2 - 292.9883*sec + 485996.4),prevTemp);
%             temp = max(prevTemp,tempIncreaseAmount+tempProfileModifier*(-0.00000000000076251*sec^4+0.000000026076*sec^3-0.00032977*sec^2+1.8379*sec-3221.2));
        elseif sec > 3600
            temp =tempIncreaseAmount+tempProfileModifier*(-4.560639E-20*sec^6 + 1.691915E-15*sec^5 - 0.00000000002564325*sec^4 + 0.0000002049408*sec^3 - 0.0009326302*sec^2 + 2.413073*sec - 2394.642);
        else
            temp = tempIncreaseAmount+max(prevTemp,tempProfileModifier*max(25,1.828874E-19*sec^6 - 1.273101E-15*sec^5 + 0.000000000002287858*sec^4 - 0.000000004988416*sec^3 + 0.00004695332*sec^2 - 0.008317037*sec + 25.22574));
        end
    else
        if sec > 5400
            temp = tempIncreaseAmount+tempProfileModifier*(-1.26231E-20*sec^6 + 6.078398E-16*sec^5 - 0.0000000000121383*sec^4 + 0.0000001292506*sec^3 - 0.0007813643*sec^2 + 2.590967*sec - 3153.787);
        elseif sec > 3525
            temp = tempIncreaseAmount+tempProfileModifier*(1.883735E-19*sec^6 - 6.794041E-15*sec^5 + 0.0000000001002798*sec^4 - 0.0000007722869*sec^3 + 0.003247263*sec^2 - 6.929942*sec + 6158.242);
        else
            temp =tempIncreaseAmount+max(prevTemp,tempProfileModifier*max(25,3.799755E-19*sec^6 - 3.320359E-15*sec^5 + 0.00000000001074309*sec^4 - 0.00000002195853*sec^3 + 0.00006102871*sec^2 - 0.002224575*sec + 29.55335));
        end
    end
end