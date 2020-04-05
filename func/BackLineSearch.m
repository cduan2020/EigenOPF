function [t, dxnorm]=BackLineSearch(ps,V0,V1,obj,tol)

nb=length(V0);
x0=[angle(V0);abs(V0)];
x1=[angle(V1);abs(V1)];
dx=x1-x0;
dxnorm=norm(dx,inf);

t=0.5;
while (t*dxnorm>tol)
    xtry=x0+t*dx;
    V=exp(1j*xtry(1:nb)).*xtry(nb+1:2*nb);
    [Asys1, Bsys1, Csys1, Dsys1]=DAEsys(ps,V);
    Afull=Asys1-Bsys1*(Dsys1\Csys1);
    lambda1= eig(full(Afull));
    abassia=max(real(lambda1(abs(lambda1)>10^-6)));
    if abassia<obj
        break;
    else
        t=t/2;
    end
end

end