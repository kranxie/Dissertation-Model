function [densityEncaps,cEncapsAvg,getkEncaps,emissivityEncaps] = getPropertiesEncaps(encapsChoice)
switch encapsChoice
    case 'steel'
        densityEncaps = 7930;
		cEncapsAvg = 558;
%         cEncapsLiquid = 590;
%         cEncapsPhaseChange =  581;
%         cEncapsSolid = 535;
        getkEncaps = @getkSteel;
        emissivityEncaps = 0.7;
    case 'ptfe' % [1] based on https://www.treborintl.com/content/properties-molded-ptfe
        % Tmelt = 327 C so would not actually work in practice, but used to
        % get an order of magnitude approximation.
        densityEncaps = 2200; % [1] 
		cEncapsAvg = 1000;
%         cEncapsLiquid = 1000; % [1] 
%         cEncapsPhaseChange =  1000; % [1] 
%         cEncapsSolid = 1000; % [1] 
        getkEncaps = @(x) 0.25; % [1] 
        emissivityEncaps = 0.92; % [2] based on https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
    case 'ceramic'
end
end

