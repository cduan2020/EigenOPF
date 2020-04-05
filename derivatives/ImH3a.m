function X=ImH3a(Cl,Yl,Kg,V,xy2dq,expdlt_a,lam)

X=H3a(Cl,Yl,Kg,V,xy2dq,expdlt_a,lam)/2j-conj(H3a(Cl,Yl,Kg,V,xy2dq,expdlt_a,conj(lam)))/2j;

end