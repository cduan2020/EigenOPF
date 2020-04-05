function X=ImEpdv(Kg,V,xy2dq,expdlt_v,lam)

X=Epdv(Kg,V,xy2dq,expdlt_v,lam)/2j-conj(Epdv(Kg,V,xy2dq,expdlt_v,conj(lam)))/2j;

end