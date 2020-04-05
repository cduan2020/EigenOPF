function X=Igdv(Yg,Kg,V,expdlt_a,expdlt_v,lam)

volt=abs(V);
theta=angle(V);
nb=length(V);
ng=size(Yg,1);
Eg = Kg*V;
xy2dq = 1j*conj(Eg./abs(Eg));


X = sparsediag(xy2dq.*(-1j*lam))*(Yg*sparsediag(exp(1j*theta))+sparsediag(Yg*V)*(-1j)*expdlt_v);

end