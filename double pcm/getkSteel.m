function [kSteel] = getkSteel(temp)
if temp > 579
    kSteel = 23.1;
elseif temp > 575
    kSteel = 22;
else
    kSteel = 17.7;
end
end

