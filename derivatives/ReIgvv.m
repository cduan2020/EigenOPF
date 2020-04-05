function X=ReIgvv(Yg,Kg,V,expdlt_a,expdlt_v,lam)

X=Igvv(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2+conj(Igvv(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2;

end