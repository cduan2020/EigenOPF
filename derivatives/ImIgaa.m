function X=ImIgaa(Yg,Kg,V,expdlt_a,expdlt_v,lam)


X=Igaa(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2j-conj(Igaa(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2j;

end