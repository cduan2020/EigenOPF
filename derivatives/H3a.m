function X=H3a(Cl,Yl,Kg,V,xy2dq,expdlt_a,lam)


X=[sparsediag(xy2dq)*Kg*sparsediag(lam.*V)-sparsediag(Kg*sparsediag(1j*V)*lam)*sparsediag(xy2dq)*(-1j)*expdlt_a; Saa(Cl,Yl,V,lam)];

end