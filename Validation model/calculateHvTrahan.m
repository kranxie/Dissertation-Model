function [hp,hv,Re] = calculateHvTrahan(G,dParticle,viscosity,Pr,~,conductivityFluid,conductivityStorage,~,voidFractionInverse)
    Re = G*dParticle./viscosity;
    Nu = 2+1.1*(Re.^0.6 .* Pr.^(1/3));
    hp = Nu.*conductivityFluid/dParticle;
    hv = 6*voidFractionInverse/dParticle.*hp;
end