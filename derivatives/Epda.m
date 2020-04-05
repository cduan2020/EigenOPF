function X=Epda(Kg,V,xy2dq,expdlt_a,lam)

X = sparsediag(xy2dq)*sparsediag(-1j*lam)*Kg*sparsediag(1j*V)+sparsediag(Kg*V)*sparsediag(-1j*lam)*sparsediag(xy2dq)*(-1j)*expdlt_a;

end