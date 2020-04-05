function X=ReH4a(Cl,Yl,Kg,V,lam,expdlt_a,xy2dq)

X=H4a(Cl,Yl,Kg,V,lam,expdlt_a,xy2dq)/2+conj(H4a(Cl,Yl,Kg,V,conj(lam),expdlt_a,xy2dq))/2;

end