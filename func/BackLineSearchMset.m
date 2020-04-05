function [t, dxnorm, Mset]=BackLineSearchMset(ps,V0,V1,obj,tol,Mset)

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
    [Ueig1,D,Veig1] = eig(full(Afull));
    lambda1=diag(D);
    
    Inx1=find(abs(lambda1)>10^-6 & imag(lambda1)>=0);
    [abassia,Inx2]=max(real(lambda1(Inx1)));
    maxInx=Inx1(Inx2);
    
    Mset.num=Mset.num+1;
    Mset.x=[Mset.x xtry];
    Mset.abassia=[Mset.abassia abassia];
    ueig=[Ueig1(:,maxInx);-(Dsys1)\(Csys1*Ueig1(:,maxInx))];
    veig=[Veig1(:,maxInx);-(Dsys1')\(Bsys1'*Veig1(:,maxInx))];
    [dsys_da, dsys_dv]=SysGradient(ps,V1,ueig,veig);
    Mset.daba_da=[Mset.daba_da ; real(dsys_da)];
    Mset.daba_dv=[Mset.daba_dv ; real(dsys_dv)];
    
    
    if abassia<obj
        break;
    else
        t=t/2;
    end
end

end