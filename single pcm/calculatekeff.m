function [keff,heff] = calculatekeff(ksolid,kfluid,fluidTemps,rockTemps,voidFraction,d,numRows,emissivityStorage,Re,Pr)
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
	kw0(numRows) = 0;
	hw0(numRows) = 0;
    heff(numRows) = 0;
    for i=1:1:numRows
        kappa = ksolid(i)/kfluid(i);
%     disp(voidFraction(i))
        hrv = (0.1952/(1+voidFraction(i)/2/(1-voidFraction(i))*(1-emissivityStorage(i))/emissivityStorage(i)))*(fluidTemps(i)/100)^3;         
        hrs = 0.1952*(emissivityStorage(i)/(2-emissivityStorage(i)))*(rockTemps(i)/100)^3;
%         hrv = (0.1952/(1+voidFractionRocks/2/(1-voidFractionRocks)*(1-0.85)/0.85))*(fluidTemps(i)/100)^3;         
%         hrs = 0.1952*(0.85/(2-0.85))*(fluidTemps(i)/100)^3;

        phi1 = 0.5*(((kappa-1)/kappa)^2*s2t1) / (log(kappa-(kappa-1)*ct1)-(kappa-1)/kappa*(1-ct1))-2/3/kappa;
        phi2 = 0.5*(((kappa-1)/kappa)^2*s2t2) / (log(kappa-(kappa-1)*ct2)-(kappa-1)/kappa*(1-ct2))-2/3/kappa;
%         phi1 = 0.5*(((kappa-1)/kappa)^2*0.6666667) / (log(kappa-(kappa-1)*0.577350269)-(kappa-1)/kappa*(1-0.577350269))-2/3/kappa;
%         phi2 = 0.5*(((kappa-1)/kappa)^2*0.1443376) / (log(kappa-(kappa-1)*0.925020234)-(kappa-1)/kappa*(1-0.925020234))-2/3/kappa;

        voidFractionWall = voidFraction(i)*1.2;
        phi = phi2 + (phi1-phi2)*(voidFraction(i)-0.26)*0.7840;
        phiWall = phi1; % based on assumption that VFwall = 1.2*VFrocks = 0.48 > VF1
        
%         phi = phi2 + (phi1-phi2)*(voidFractionRocks-0.26)/(0.216);
%         phi = phi2 + (phi1-phi2)*(voidFractionRocks-VF2)/(VF1-VF2);
%         disp(phi)
        keff(i) = kfluid(i)*(voidFraction(i)*(1+beta*hrv*d(i)/kfluid(i))+beta*(1-voidFraction(i))/((1/phi+hrs*d(i)/kfluid(i))+gamma/kappa));
%         keff(i) = kFluid(i)*(voidFractionRocks*(1+0.9*hrv*d/kFluid(i))+0.9*(1-voidFractionRocks)/(1/(1/phi)+hrs*d/kFluid(i))+0.66667*(1/kappa));
    
        kw0(i) = kfluid(i)*(voidFractionWall*(2+hrv*d(i)/kfluid(i)) + ((1-voidFractionWall)/(1/(1/phiWall+hrs*d(i)/kfluid(i))+1/kappa/3)));
		hw0(i) = 1/(d(i)/kw0(i) - 0.5*d(i)/keff(i));
		heff(i)= hw0(i) + 0.0205*Re(i)*Pr(i);
		% alpha_w = 0.5 * alpha, where alpha = 0.0205 for cylindrical cross section
	end
end