function conductivity = getConductivityDryAir(temp)
conductivity = 7.7053418E-19*temp.^6 - 2.1725308E-15*temp.^5 + 0.000000000002454112*temp.^4 - 0.0000000014263403*temp.^3 + 0.00000045399548*temp.^2 - 0.000030159064*temp + 0.041669106;
end