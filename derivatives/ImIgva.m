function X=ImIgva(Yg,Kg,V,expdlt_a,expdlt_v,lam)


X=Igva(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2j-conj(Igva(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2j;

end