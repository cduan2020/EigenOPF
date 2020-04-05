function X=H4v(Cl,Yl,Kg,V,lam,expdlt_v,xy2dq)

theta=angle(V);

X=[-sparsediag(Kg*sparsediag(exp(1j*theta))*lam)*sparsediag(xy2dq)*(-1j)*expdlt_v; Svv(Cl,Yl,V,lam)];

end