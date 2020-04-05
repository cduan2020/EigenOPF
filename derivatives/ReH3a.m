function X=ReH3a(Cl,Yl,Kg,V,xy2dq,expdlt_a,lam)

X=H3a(Cl,Yl,Kg,V,xy2dq,expdlt_a,lam)/2+conj(H3a(Cl,Yl,Kg,V,xy2dq,expdlt_a,conj(lam)))/2;

end