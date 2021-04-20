function [specificheat, conductivity, viscosity,density,enthalpy,enthalpyDerivative] = getMoistAirProperties(RHpercent)
% temp, airPressure, RHpercent
    switch RHpercent
% ********* enthalpy derivatives haven't been updated for RH > 50 *********
        case 50
            enthalpy = @(temp) 0.1107*temp.^2 + 991.4*temp + 296300;
            enthalpyDerivative = @(temp) 0.2214*temp + 991.4;
            specificheat = @(temp) -3.425e-7*temp.^3 + 4.940e-4*temp.^2 + 0.01043e-2*temp + 1014;
            conductivity = @(temp) -1.387e-8*temp.^2 + 6.819e-5*temp + 2.438e-2;
            viscosity = @(temp) -1.189e-11*temp.^2 + 4.411e-8*temp + 1.755e-5;
            density = @(temp) 3.899e-12*temp.^4 - 8.854e-9*temp.^3 + 7.946e-6*temp.^2 -3.828e-3*temp + 1.238;
        case 60
            enthalpy = @(temp) 0.1111158*temp.^2+993.0278*temp-24151.51;
            enthalpyDerivative = @(temp)  0.2222316*temp+993.0278;
            specificheat = @(temp) -3.475E-07*temp.^3 + 5.014E-04*temp.^2 + 7.995E-03*temp + 1.016E+03;
            conductivity = @(temp)-1.376E-08*temp.^2 + 6.822E-05*temp + 2.436E-02;
            viscosity = @(temp) -1.186E-11*temp.^2 + 4.411E-08*temp + 1.753E-05;
            density = @(temp)3.879E-12*temp.^4 - 8.809E-09*temp.^3 + 7.905E-06*temp.^2 - 3.809E-03*temp + 1.232E+00;
        case 70
            enthalpy = @(temp)1.115403E-01*temp.^2 + 9.946940E+02*temp - 2.418984E+04;
            enthalpyDerivative = @(temp) 2.23081E-01*temp+9.946940E+02;
            specificheat = @(temp)-3.483E-07*temp.^3 + 5.029E-04*temp.^2 + 8.062E-03*temp + 1.017E+03;
            conductivity = @(temp) -1.365E-08*temp.^2 + 6.824E-05*temp + 2.434E-02;
            viscosity = @(temp) -1.184E-11*temp.^2 + 4.410E-08*temp + 1.751E-05;
            density = @(temp)3.859E-12*temp.^4 - 8.764E-09*temp.^3 + 7.864E-06*temp.^2 - 3.789E-03*temp + 1.225E+00;
        case 80
            enthalpy = @(temp)1.119661E-01*temp.^2 + 9.963723E+02*temp - 2.422812E+04;
            enthalpyDerivative = @(temp) 2.23932E-01*temp+9.963723E+02;
            specificheat = @(temp)-3.493E-07*temp.^3 + 5.047E-04*temp.^2 + 8.031E-03*temp + 1.019E+03;
            conductivity = @(temp)-1.354E-08*temp.^2 + 6.827E-05*temp + 2.432E-02;
            viscosity = @(temp) -1.181E-11*temp.^2 + 4.409E-08*temp + 1.749E-05;
            density = @(temp)3.839E-12*temp.^4 - 8.718E-09*temp.^3 + 7.824E-06*temp.^2 - 3.770E-03*temp + 1.219E+00;
        case 90
            enthalpy = @(temp)1.123932E-01*temp.^2 + 9.980629E+02*temp - 2.426633E+04;
            enthalpyDerivative = @(temp) 2.24786E-01*temp+9.980629E+02;
            specificheat = @(temp)-3.505E-07*temp.^3 + 5.067E-04*temp.^2 + 7.888E-03*temp + 1.021E+03;
            conductivity = @(temp) -1.343E-08*temp.^2 + 6.829E-05*temp + 2.431E-02;
            viscosity = @(temp) -1.179E-11*temp.^2 + 4.408E-08*temp + 1.747E-05;
            density = @(temp) 3.819E-12*temp.^4 - 8.673E-09*temp.^3 + 7.783E-06*temp.^2 - 3.750E-03*temp + 1.213E+00;
        case 100
            enthalpy = @(temp)1.128214E-01*temp.^2 + 9.997660E+02*temp - 2.430445E+04;
            enthalpyDerivative = @(temp) 2.25643E-01*temp+9.997660E+02;
            specificheat = @(temp)-3.518E-07*temp.^3 + 5.090E-04*temp.^2 + 7.623E-03*temp + 1.023E+03;
            conductivity = @(temp)-1.332E-08*temp.^2 + 6.832E-05*temp + 2.429E-02;
            viscosity = @(temp)-1.176E-11*temp.^2 + 4.407E-08*temp + 1.746E-05;
            density = @(temp)3.799E-12*temp.^4 - 8.628E-09*temp.^3 + 7.743E-06*temp.^2 - 3.731E-03*temp + 1.206E+00;
        otherwise
            enthalpy = 0;
            enthalpyDerivative = 0;
            specificheat = 0;
            conductivity = 0;
            viscosity = 0;
            density = 0;
    end
end