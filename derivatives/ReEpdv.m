function X=ReEpdv(Kg,V,xy2dq,expdlt_v,lam)

X=Epdv(Kg,V,xy2dq,expdlt_v,lam)/2+conj(Epdv(Kg,V,xy2dq,expdlt_v,conj(lam)))/2;

end