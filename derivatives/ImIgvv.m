function X=ImIgvv(Yg,Kg,V,expdlt_a,expdlt_v,lam)


X=Igvv(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2j-conj(Igvv(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2j;

end