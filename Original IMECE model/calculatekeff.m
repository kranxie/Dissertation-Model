function [keff,heff] = calculatekeff(kRocks,kFluid,fluidTemps,rockTemps,voidFractionRocks,d,numRows,emissivityRocks)
%     inputs = [    1.5000,   30.7870 ,   0.4000  ,  0.0320];

    fluidTemps(:) =fluidTemps(:)+273.15;   
    rockTemps(:) =rockTemps(:)+273.15;   
    
    
    beta = 0.9;
    gamma = 0.66666666667;
%     voidFractionRocks1 = 0.476;
%     voidFractionRocks2 = 0.26;
%     VF1 = 0.476;
%     VF2 = 0.26;
    s2t1 = 0.666666666666667;   % 1/.5
    s2t2 = 0.144337567297406;   % 1/4/sqrt(3)
    ct1 = 0.577350269189626;    % cos(asin(sqrt(s2t1)))
    ct2 = 0.925020233671996;    % cos(asin(sqrt(s2t2)))

    % keffinputs = [kRocks, fluidTemps(i), voidFractionRocks, rockAverageDiameter];
    keff(numRows) = 0;
    heff(numRows) = 0;
    for i=1:1:numRows
        kappa = kRocks(i)/kFluid(i);
    
        hrv = (0.1952/(1+voidFractionRocks/2/(1-voidFractionRocks)*(1-emissivityRocks)/emissivityRocks))*(fluidTemps(i)/100)^3;         
        hrs = 0.1952*(emissivityRocks/(2-emissivityRocks))*(rockTemps(i)/100)^3;
%         hrv = (0.1952/(1+voidFractionRocks/2/(1-voidFractionRocks)*(1-0.85)/0.85))*(fluidTemps(i)/100)^3;         
%         hrs = 0.1952*(0.85/(2-0.85))*(fluidTemps(i)/100)^3;

        phi1 = 0.5*(((kappa-1)/kappa)^2*s2t1) / (log(kappa-(kappa-1)*ct1)-(kappa-1)/kappa*(1-ct1))-2/3/kappa;
        phi2 = 0.5*(((kappa-1)/kappa)^2*s2t2) / (log(kappa-(kappa-1)*ct2)-(kappa-1)/kappa*(1-ct2))-2/3/kappa;
%         phi1 = 0.5*(((kappa-1)/kappa)^2*0.6666667) / (log(kappa-(kappa-1)*0.577350269)-(kappa-1)/kappa*(1-0.577350269))-2/3/kappa;
%         phi2 = 0.5*(((kappa-1)/kappa)^2*0.1443376) / (log(kappa-(kappa-1)*0.925020234)-(kappa-1)/kappa*(1-0.925020234))-2/3/kappa;

        voidFractionWall = voidFractionRocks*1.2;
        phi = phi2 + (phi1-phi2)*(voidFractionRocks-0.26)*0.7840;
        phiWall = phi1; % based on assumption that VFwall = 1.2*VFrocks = 0.48 > VF1
        
%         phi = phi2 + (phi1-phi2)*(voidFractionRocks-0.26)/(0.216);
%         phi = phi2 + (phi1-phi2)*(voidFractionRocks-VF2)/(VF1-VF2);

        keff(i) = kFluid(i)*(voidFractionRocks*(1+beta*hrv*d/kFluid(i))+beta*(1-voidFractionRocks)/(1/phi+hrs*d/kFluid(i))+gamma/kappa);
%         keff(i) = kFluid(i)*(voidFractionRocks*(1+0.9*hrv*d/kFluid(i))+0.9*(1-voidFractionRocks)/(1/(1/phi)+hrs*d/kFluid(i))+0.66667*(1/kappa));
    
        heff(i) = kFluid(i) *(voidFractionWall*(2+hrv*d/kFluid(i)) + ((1-voidFractionWall)/(1/(1/phiWall+hrs*d/kFluid(i))+1/kappa/3)));
    end
end