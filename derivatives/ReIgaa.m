function X=ReIgaa(Yg,Kg,V,expdlt_a,expdlt_v,lam)

X=Igaa(Yg,Kg,V,expdlt_a,expdlt_v,lam)/2+conj(Igaa(Yg,Kg,V,expdlt_a,expdlt_v,conj(lam)))/2;

end