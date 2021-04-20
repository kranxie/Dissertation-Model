function [hp,hv,Re] = calculateHvPCMZanganeh(G,dEncap,viscosity,Pr,PrS,kFluid,~,surfaceAreaPCM,~)
    Re = 3.4658*G*dEncap./viscosity;
    Nu = 0.4845 * sqrt(Re).*Pr.^0.37.*(Pr./PrS).^0.25;
    hp = Nu.*kFluid./dEncap;
    hv = hp.*surfaceAreaPCM;
end

