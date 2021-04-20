function [hp,hv,Re] = calculateHvPCMBeasley(G,dEncap,viscosity,Pr,~,conductivity,~,~,voidFractionPCMInverse)
    Re = G*dEncap./viscosity;
    Nu = 2 + 2.03*Re.^(0.5).*Pr.^(1/3) + 0.049.*Re.*Pr.^0.5;          % From Beasley 1984 and Beasley 1989
    
    hp = Nu.*conductivity/dEncap;
    hv = 6.*hp.*voidFractionPCMInverse/dEncap;
end

