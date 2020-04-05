function X=ReH3v(Cl,Yl,Kg,V,xy2dq,expdlt_v,lam)

X=H3v(Cl,Yl,Kg,V,xy2dq,expdlt_v,lam)/2+conj(H3v(Cl,Yl,Kg,V,xy2dq,expdlt_v,conj(lam)))/2;

end