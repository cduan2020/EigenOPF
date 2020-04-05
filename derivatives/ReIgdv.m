function X=ReIgdv(Yg,Kg,V,expdlt_a,expdlt_v,lam)

X=Igdv(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2+conj(Igdv(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2;

end