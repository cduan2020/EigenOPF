function X=Igda(Yg,Kg,V,expdlt_a,expdlt_v,lam)

volt=abs(V);
theta=angle(V);
nb=length(V);
ng=size(Yg,1);
Eg = Kg*V;
xy2dq = 1j*conj(Eg./abs(Eg));


X = sparsediag(xy2dq.*(-1j*lam))*(Yg*sparsediag(1j*V)+sparsediag(Yg*V)*(-1j)*expdlt_a);

end