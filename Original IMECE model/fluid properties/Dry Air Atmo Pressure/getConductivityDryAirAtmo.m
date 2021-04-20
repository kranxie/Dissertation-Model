function conductivity = getConductivityDryAirAtmo(temp)
conductivity = 9.2603898E-18*temp.^5 - 2.9179175E-14*temp.^4 + 4.2744126E-11*temp.^3 - 4.2635453E-08*temp.^2 + 7.5290320E-05*temp + 2.3997830E-02;
end