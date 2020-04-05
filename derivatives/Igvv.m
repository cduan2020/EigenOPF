function X=Igvv(Yg,Kg,V,expdlt_a,expdlt_v,lam)

volt=abs(V);
theta=angle(V);
nb=length(V);
ng=size(Yg,1);
Eg = Kg*V;
xy2dq = 1j*conj(Eg./abs(Eg));


X = sparsediag(Yg*(exp(1j*theta).*lam))*sparsediag(xy2dq)*(-1j)*expdlt_v;




end