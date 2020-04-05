function X=ReIgva(Yg,Kg,V,expdlt_a,expdlt_v,lam)

X=Igva(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2+conj(Igva(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2;

end