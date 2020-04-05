function X=ImH3v(Cl,Yl,Kg,V,xy2dq,expdlt_v,lam)

X=H3v(Cl,Yl,Kg,V,xy2dq,expdlt_v,lam)/2j-conj(H3v(Cl,Yl,Kg,V,xy2dq,expdlt_v,conj(lam)))/2j;

end