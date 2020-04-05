function X=ImIgdv(Yg,Kg,V,expdlt_a,expdlt_v,lam)


X=Igdv(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2j-conj(Igdv(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2j;

end