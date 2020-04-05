x0=[angle(V0);abs(V0)];
x1=[angle(V1);abs(V1)];
dx=x1-x0;

[Asys0, Bsys0, Csys0, Dsys0]=DAEsys(ps,V0);
Afull=Asys0-Bsys0*(Dsys0\Csys0);
[Ueig0,D,Veig0] = eig(full(Afull));
lambda=diag(D);

Inx1=find(abs(lambda)>10^-6 & imag(lambda)>=0);
[abassia,Inx2]=max(real(lambda(Inx1)));
maxInx=Inx1(Inx2);

ueig=[Ueig0(:,maxInx);-inv(Dsys0)*Csys0*Ueig0(:,maxInx)];
veig=[Veig0(:,maxInx);-inv(Dsys0)'*Bsys0'*Veig0(:,maxInx)];
 [dsys_da, dsys_dv]=SysGradient(ps,V0,ueig,veig);
 
 
 [real(dsys_da), real(dsys_dv)]*dx