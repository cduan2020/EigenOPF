function X=Igaa(Yg,Kg,V,expdlt_a,expdlt_v,lam)

volt=abs(V);
theta=angle(V);
nb=length(V);
ng=size(Yg,1);
Eg = Kg*V;
xy2dq = 1j*conj(Eg./abs(Eg));


X = sparsediag(xy2dq)*Yg*sparsediag(-lam.*V)+ sparsediag(Yg*(1j*V.*lam))*sparsediag(xy2dq)*(-1j)*expdlt_a;


end