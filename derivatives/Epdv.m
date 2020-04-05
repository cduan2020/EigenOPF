function X=Epdv(Kg,V,xy2dq,expdlt_v,lam)

X =  sparsediag(xy2dq)*sparsediag(-1j*lam)*Kg*sparsediag(1j*V)+sparsediag(Kg*V)*sparsediag(-1j*lam)*sparsediag(xy2dq)*(-1j)*expdlt_v;

end