function density = getDensityTrahan(temp)
% density = 101.3e3./287.1./(temp+273.15);
% density = 25e6./287.1./(temp+273.15);

density = (-5.75399e-16)*(temp^5) + (3.02846e-12)*(temp^4)-(6.18352e-9)*(temp^3)+(6.299277e-6)*(temp^2)-(3.5422e-3)*temp+1.25079;
    % kg/m3, temp in Celcius?

end