function [TmeltStart, TmeltEnd, dTMelt, cPCMSolid, cPCMLiquid, enthalpyFusion, kPCMSolid, kPCMLiquid,densityPCMSolid,densityPCMLiquid] = getPropertiesPCM(pcmChoice)

switch pcmChoice
    case 'AlSi12'
        TmeltStart = 575;       dTMelt = 4;         TmeltEnd = TmeltStart + dTMelt;
        cPCMSolid = 1070;       cPCMLiquid = 1170; 	enthalpyFusion = 466000;
        kPCMSolid = 160;        kPCMLiquid = 160;
        densityPCMSolid = 2650; densityPCMLiquid = 2650;
    case 'KCl-NaCl-Tm=657' 
        Tmelt = 657;   % DOI:10.1115/ES2015-49460         
        dTMelt = 4; % assumed
        TmeltStart = Tmelt-dTMelt/2;                    TmeltEnd = TmeltStart + dTMelt;
        
        cPCMSolid = 1835.4;       cPCMLiquid = 1655; % assumed same as NaNO3 from Trahan diss. due to difficulty finding data
        enthalpyFusion = 360000; % DOI:10.1115/ES2015-49460
        kPCMSolid = 0.5;        kPCMLiquid = 0.5; % assumed same as NaNO3 from Trahan diss. due to difficulty finding data
        densityPCMSolid = 1592; densityPCMLiquid = 1592;  % Based on 41.2 KCl-NaCl (mol%) at T=555.2 from ref below
        % Van Artsdalen, E. R., & Yaffe, I. S. (1955). Electrical conductance and density of molten salt systems: KCl–LiCl, KCl–NaCl and KCl–KI. The Journal of Physical Chemistry, 59(2), 118-127.  
    case 'KCl-NaCl-Tm=680'
        Tmelt = 680;   % Increased for the sake of experimentation
        dTMelt = 4; % assumed
        TmeltStart = Tmelt-dTMelt/2;                    TmeltEnd = TmeltStart + dTMelt;
        
        cPCMSolid = 1835.4;       cPCMLiquid = 1655; % assumed same as NaNO3 from Trahan diss. due to difficulty finding data
        enthalpyFusion = 360000; % DOI:10.1115/ES2015-49460
        kPCMSolid = 0.5;        kPCMLiquid = 0.5; % assumed same as NaNO3 from Trahan diss. due to difficulty finding data
        densityPCMSolid = 1592; densityPCMLiquid = 1592;  % Based on 41.2 KCl-NaCl (mol%) at T=555.2 from ref below
        % Van Artsdalen, E. R., & Yaffe, I. S. (1955). Electrical conductance and density of molten salt systems: KCl–LiCl, KCl–NaCl and KCl–KI. The Journal of Physical Chemistry, 59(2), 118-127.  
	case 'KCl-NaCl-Tm=670'
        Tmelt = 670;   % Increased for the sake of experimentation
        dTMelt = 4; % assumed
        TmeltStart = Tmelt-dTMelt/2;                    TmeltEnd = TmeltStart + dTMelt;
        
        cPCMSolid = 1835.4;       cPCMLiquid = 1655; % assumed same as NaNO3 from Trahan diss. due to difficulty finding data
        enthalpyFusion = 360000; % DOI:10.1115/ES2015-49460
        kPCMSolid = 0.5;        kPCMLiquid = 0.5; % assumed same as NaNO3 from Trahan diss. due to difficulty finding data
        densityPCMSolid = 1592; densityPCMLiquid = 1592;  % Based on 41.2 KCl-NaCl (mol%) at T=555.2 from ref below
        % Van Artsdalen, E. R., & Yaffe, I. S. (1955). Electrical conductance and density of molten salt systems: KCl–LiCl, KCl–NaCl and KCl–KI. The Journal of Physical Chemistry, 59(2), 118-127.  
end

pcmProperties = [TmeltStart, TmeltEnd, dTMelt, cPCMSolid, cPCMLiquid, enthalpyFusion, kPCMSolid, kPCMLiquid,densityPCMSolid,densityPCMLiquid];


end