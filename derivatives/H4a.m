function X=H4a(Cl,Yl,Kg,V,lam,expdlt_a,xy2dq)

theta=angle(V);

X=[sparsediag(xy2dq)*Kg*sparsediag((-1j*exp(1j*theta)).*lam)-sparsediag(Kg*sparsediag(exp(1j*theta))*lam)*sparsediag(xy2dq)*(-1j)*expdlt_a; Sva(Cl,Yl,V,lam)];

end