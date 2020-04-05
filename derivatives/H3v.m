function X=H3v(Cl,Yl,Kg,V,xy2dq,expdlt_v,lam)

theta=angle(V);

X=[sparsediag(xy2dq)*Kg*sparsediag(-lam.*(1j*exp(1j*theta)))-sparsediag(Kg*sparsediag(1j*V)*lam)*sparsediag(xy2dq)*(-1j)*expdlt_v; Sav(Cl,Yl,V,lam)];

end