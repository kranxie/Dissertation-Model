function [getSHSEnergy,getSHSSpecificHeat,getkSHS,densitySHS,emissivitySHS] = getPropertiesSHS(SHSChoice)

switch SHSChoice
    case 'Zanganeh2015' % Mixed rocks
        getSHSEnergy = @(tempC) 747.0995*tempC+0.2838*tempC.^2;
        getSHSSpecificHeat = @(tempC) 747.0995 + 0.5676*tempC;
        getkSHS = @(tempC)5+(1-5)*(tempC-25)/(700-25);
        densitySHS = 2635;
        emissivitySHS = 0.83;
    case 'Trahan' % Taconite ore
        % Jamie value
        % raising bed from 25 to 500 C takes input of 60.8 kWh, far higher
        % than the 36 kWh specified by Jamie's dissertation
        getSHSEnergy = @(temp) 1.2953375e-06*(temp+273).^4 - 0.002526691666667*(temp+273).^3 + 2.942740725*(temp+273).^2 - 429.234632*(temp+273);
        getSHSSpecificHeat = @(temp) 5.18135e-6*(temp+273).^3 - 7.580075e-3*(temp+273).^2 + 5.88548145*(temp+273) - 429.234632;
        getkSHS = @(x) 1.5;
        densitySHS = 3200;                                        % kg/m3
        emissivitySHS = 0.83;
    case 'NIST' % Taconite ore
%                 ((temp+273)/1000) = t
        getSHSSpecificHeat = @(temp) 6.262211312058515*(93.43834+108.3577*((temp+273)/1000)-50.86447*((temp+273)/1000).^2+25.58683*((temp+273)/1000).^3-1.61133*((temp+273)/1000).^(-2));
        getSHSEnergy = @(temp) (749.8980157*temp) + 10090488.95*((temp + 273).^(-1)) + (0.270235027*temp.^2) - (0.0000624319*temp.^3) + (0.0000000400575*temp.^4);
        getkSHS = @(x) 1.5;
        densitySHS = 3200;                                        % kg/m3
        emissivitySHS = 0.83;
    case 'Smetana2013' % Taconite ore
        getSHSSpecificHeat = @(temp) 2.42E-11*(temp+273).^4 - 6.71E-08*(temp+273).^3 + 6.96E-05*(temp+273).^2 - 3.15E-02*(temp+273) + 5.78E+00;
        getSHSEnergy = @(temp) 2.42E-11/5*(temp+273).^5 - 6.71E-08/4*(temp+273).^4 + 6.96E-05/3*(temp+273).^3 - 3.15E-02/2*(temp+273)^2 + 5.78E+00*(temp);
        getkSHS = @(x) 1.5;
        densitySHS = 3200;                                        % kg/m3
        emissivitySHS = 0.83;
    otherwise
        getSHSEnergy = 0;
        getSHSSpecificHeat = 0;
        getkSHS = 0;
        densitySHS = 0;
        emissivitySHS = 0;
end


end