function X=ImH4v(Cl,Yl,Kg,V,lam,expdlt_v,xy2dq)

X=H4v(Cl,Yl,Kg,V,lam,expdlt_v,xy2dq)/2j-conj(H4v(Cl,Yl,Kg,V,conj(lam),expdlt_v,xy2dq))/2j;

end