function X=ImIgda(Yg,Kg,V,expdlt_a,expdlt_v,lam)


X=Igda(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2j-conj(Igda(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2j;

end