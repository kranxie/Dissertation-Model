function [getEnthalpy,getSpecificHeat,getConductivity,getViscosity,getDensity,fluidPressure,getEntropy] = getPropertiesFluid(fluidName)
switch fluidName
    case 'CO2Atmo'
        fluidPressure = 101.3e3;
        getEnthalpy = @getEnthalpyCO2Atmo;
        getConductivity = @getConductivityCO2Atmo;
        getViscosity = @getViscosityCO2Atmo;
        getDensity = @getDensityCO2Atmo;
        getSpecificHeat = @getCpCO2Atmo;
    case 'CO2LP'        
        fluidPressure=12.38e6;
        getEnthalpy = @getEnthalpyCO2LP;
        getSpecificHeat = @getCpCO2LP;
        getConductivity = @getConductivityCO2LP;
        getViscosity = @getViscosityCO2LP;
        getDensity = @getDensityCO2LP;
    case 'CO2HP'
        fluidPressure=25e6;
        getEnthalpy = @getEnthalpyCO2HP;
        getSpecificHeat = @getCpCO2HP;
        getConductivity = @getConductivityCO2HP;
        getViscosity = @getViscosityCO2HP;
        getDensity = @getDensityCO2HP;
		getEntropy = @getEntropyCO2HP;
    case 'DryAirHP'
        fluidPressure = 25e6;
        getEnthalpy = @getEnthalpyDryAir;
        getSpecificHeat = @getCpDryAir;
        getConductivity = @getConductivityDryAir;
        getViscosity = @getViscosityDryAir;
        getDensity = @getDensityDryAir;
    case 'DryAirAtmo'
        fluidPressure = 101.3e3;
        getEnthalpy = @getEnthalpyDryAirAtmo;
        getConductivity = @getConductivityDryAirAtmo;
        getViscosity = @getViscosityDryAirAtmo;
        getDensity = @getDensityDryAirAtmo;
        getSpecificHeat = @getCpDryAirAtmo;
		getEntropy = @getEntropyDryAirAtmo;
    case 'DryAirTrahan'
        fluidPressure = 101.3e3;
        getDensity = @(Tave) ((-5.75399E-16)*(Tave.^5))+((3.02846E-12)*(Tave.^4))-((6.18352E-9)*(Tave.^3))+((6.29927E-6)*(Tave.^2))-((3.5422E-3)*Tave)+1.25079; %density of air (kg/m3) 
        getViscosity = @(Tave) (((6.10504E-10)*(Tave.^3))-((2.13036E-6)*(Tave.^2))+((4.71398E-3)*(Tave))+1.67555)*(10.^-5); %Dynamic viscosity of air (kg/ms) 
        getSpecificHeat = @(Tave) (((1.28806E-13)*(Tave.^4))-((4.46054E-10)*(Tave.^3))+((4.8772E-7)*(Tave.^2))+((1.82754E-5)*Tave)+1.00651)*1000; %Cp of air (J/kg-K) 
        getConductivity = @(Tave) ((-4.44955E-15)*(Tave.^4))+((2.41702E-11)*(Tave.^3))-((4.09601E-8)*(Tave.^2))+((7.91034E-5)*Tave)+.0242006; % thermal conductivity of air (W/mK)
        getEnthalpy = @(Tave) (((1.28806E-13)*(Tave.^5)/5)-((4.46054E-10)*(Tave.^4)/4)+((4.8772E-7)*(Tave.^3)/3)+((1.82754E-5)*Tave.^2)/2+1.00651*Tave)*1000; %Cp of air (J/kg-K) 
    otherwise
        disp('Invalid fluid. Terminating operation.')
        return
end
end

