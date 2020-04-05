function X=ReEpda(Kg,V,xy2dq,expdlt_a,lam)

X=Epda(Kg,V,xy2dq,expdlt_a,lam)/2+conj(Epda(Kg,V,xy2dq,expdlt_a,conj(lam)))/2;

end