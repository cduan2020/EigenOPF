function X=ReIgda(Yg,Kg,V,expdlt_a,expdlt_v,lam)

X=Igda(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2+conj(Igda(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2;

end