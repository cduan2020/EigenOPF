function X=ImEpda(Kg,V,xy2dq,expdlt_a,lam)

X=Epda(Kg,V,xy2dq,expdlt_a,lam)/2j-conj(Epda(Kg,V,xy2dq,expdlt_a,conj(lam)))/2j;

end