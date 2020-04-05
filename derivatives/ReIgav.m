function X=ReIgav(Yg,Kg,V,expdlt_a,expdlt_v,lam)

X=Igav(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2+conj(Igav(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2;

end