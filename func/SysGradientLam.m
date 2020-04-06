function [dsys_da, dsys_dv]=SysGradientLam(ps,V,ueig,veig)

%% import dynamical parameters
nb=size(ps.bus,1);
ng=size(ps.gen,1);
nl=nb-ng;
baseMVA = ps.baseMVA;  %common system base
H = ps.gen_dyn(:,5);   %inertia constants for generators
D = ps.gen_dyn(:,6);    % damping constants for generators
Ka = ps.gen_dyn(:,9);
if isfield(ps, 'ref_freq')
    omega_R = 2 * pi * ps.ref_freq;
else
    omega_R = 2 * pi * 60;    
end
Ks = ps.gen_dyn(:,10)/omega_R;

Td0p =ps.gen_dyn(:,7);
Tq0p =ps.gen_dyn(:,8);

Xd = ps.gen_dyn(:,2);
Xdp = ps.gen_dyn(:,1);
Xq = ps.gen_dyn(:,4);
Xqp = ps.gen_dyn(:,3);

D = sparse(1:ng, 1:ng, D/omega_R , ng, ng);
M = sparse(1:ng, 1:ng, 2*H/omega_R, ng, ng);
invM=sparse(1:ng, 1:ng, omega_R./(2*H), ng, ng);
T=[speye(4*ng), sparse(4*ng,2*nb); sparse(2*ng+2*nl,4*ng) sparse(2*ng+2*nl,2*nb)];

%% import network parameters

gbus=ps.gen(:,1);
lbus=setdiff(1:nb,gbus);
Ybus=ps.Y;
Yg=Ybus(gbus,:);
Yl=Ybus(lbus,:);
Cg=sparse(1:ng,gbus,ones(ng,1),ng,nb);
Cl=sparse(1:nl,lbus,ones(nl,1),nl,nb);

%% calculat system matrices
V=sparse(V);
volt=abs(V);
theta=angle(V);

%=================== additional useful values =============================
Kg = @(X) Cg+sparsediag(1j*X)*Yg;
xy2dq = 1j*conj(Kg(Xq)*V./abs(Kg(Xq)*V));
Eg = Kg(Xq)*V; 

expdlt_a = spdiags(1j*Eg./abs(Eg),0,ng,ng)*imag(spdiags(1./Eg,0,ng,ng)*Kg(Xq)*spdiags(1j*V,0,nb,nb));
expdlt_v = spdiags(1j*Eg./abs(Eg),0,ng,ng)*imag(spdiags(1./Eg,0,ng,ng)*Kg(Xq)*spdiags(V./abs(V),0,nb,nb));

ueig=sparse(ueig);
veig=sparse(veig);

Index1_R=1:ng;
Index2_R=ng+1:2*ng;
Index3_R=2*ng+1:3*ng;
Index4_R=3*ng+1:4*ng;
Index5_R=4*ng+1:4*ng+nb;
Index6_R=4*ng+nb+1:4*ng+2*nb;

Index1_L=1:ng;
Index2_L=ng+1:2*ng;
Index3_L=2*ng+1:3*ng;
Index4_L=3*ng+1:4*ng;
Index5_L=4*ng+1:5*ng;
Index6_L=5*ng+1:5*ng+nl;
Index7_L=5*ng+nl+1:6*ng+nl;
Index8_L=6*ng+nl+1:6*ng+2*nl;

a31_a= veig(Index3_L)'*(-sparsediag((Xd-Xdp)./Td0p)*ReIgda(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index1_R)));
a31_v= veig(Index3_L)'*(-sparsediag((Xd-Xdp)./Td0p)*ReIgdv(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index1_R)));
a41_a= veig(Index4_L)'*(sparsediag((Xq-Xqp)./Tq0p)*ImIgda(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index1_R)));
a41_v= veig(Index4_L)'*(sparsediag((Xq-Xqp)./Tq0p)*ImIgdv(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index1_R)));

b21_a=-veig(Index2_L)'*invM*ReSaa(Cg,Yg,V,ueig(Index5_R));
b21_v=-veig(Index2_L)'*invM*ReSav(Cg,Yg,V,ueig(Index5_R));
b22_a=-veig(Index2_L)'*invM*ReSva(Cg,Yg,V,ueig(Index6_R));
b22_v=-veig(Index2_L)'*invM*ReSvv(Cg,Yg,V,ueig(Index6_R));

b31_a = veig(Index3_L)'*(-sparsediag((Xd-Xdp)./Td0p)*ReIgaa(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index5_R)));
b31_v = veig(Index3_L)'*(-sparsediag((Xd-Xdp)./Td0p)*ReIgav(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index5_R)));
b32_a = veig(Index3_L)'*(-sparsediag((Xd-Xdp)./Td0p)*ReIgva(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index6_R)));
b32_v = veig(Index3_L)'*(-sparsediag((Xd-Xdp)./Td0p)*ReIgvv(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index6_R)));
b41_a = veig(Index4_L)'*(sparsediag((Xq-Xqp)./Tq0p)*ImIgaa(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index5_R)));
b41_v = veig(Index4_L)'*(sparsediag((Xq-Xqp)./Tq0p)*ImIgav(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index5_R)));
b42_a = veig(Index4_L)'*(sparsediag((Xq-Xqp)./Tq0p)*ImIgva(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index6_R)));
b42_v = veig(Index4_L)'*(sparsediag((Xq-Xqp)./Tq0p)*ImIgvv(Yg,Kg(Xq),V,expdlt_a,expdlt_v,ueig(Index6_R)));

c_11_a = -veig(Index5_L)'*ReEpda(Kg(Xqp),V,xy2dq,expdlt_a,ueig(Index1_R));
c_11_v = -veig(Index5_L)'*ReEpdv(Kg(Xqp),V,xy2dq,expdlt_v,ueig(Index1_R));
c_31_a = -veig(Index7_L)'*ImEpda(Kg(Xdp),V,xy2dq,expdlt_a,ueig(Index1_R));
c_31_v = -veig(Index7_L)'*ImEpdv(Kg(Xdp),V,xy2dq,expdlt_v,ueig(Index1_R)); 



lam=ueig(Index5_R);
d11_a=  veig([Index5_L Index6_L])'*ReH3a(Cl,Yl,Kg(Xqp),V,xy2dq,expdlt_a,lam);
d11_v= veig([Index5_L Index6_L])'*ReH3v(Cl,Yl,Kg(Xqp),V,xy2dq,expdlt_v,lam);


lam=ueig(Index5_R);
d21_a=  veig([Index7_L Index8_L])'*ImH3a(Cl,Yl,Kg(Xdp),V,xy2dq,expdlt_a,lam);
d21_v= veig([Index7_L Index8_L])'*ImH3v(Cl,Yl,Kg(Xdp),V,xy2dq,expdlt_v,lam);


lam=ueig(Index6_R);
d12_a=  veig([Index5_L Index6_L])'*ReH4a(Cl,Yl,Kg(Xqp),V,lam,expdlt_a,xy2dq);
d12_v= veig([Index5_L Index6_L])'*ReH4v(Cl,Yl,Kg(Xqp),V,lam,expdlt_v,xy2dq);


lam=ueig(Index6_R);
d22_a=  veig([Index7_L Index8_L])'*ImH4a(Cl,Yl,Kg(Xdp),V,lam,expdlt_a,xy2dq);
d22_v= veig([Index7_L Index8_L])'*ImH4v(Cl,Yl,Kg(Xdp),V,lam,expdlt_v,xy2dq);


dsys_da=(a31_a+a41_a+b21_a+b22_a+b31_a+b32_a+b41_a+b42_a+c_11_a+c_31_a+d11_a+d21_a+d12_a+d22_a)/(veig'*T*ueig);
dsys_dv=(a31_v+a41_v+b21_v+b22_v+b31_v+b32_v+b41_v+b42_v+c_11_v+c_31_v+d11_v+d21_v+d12_v+d22_v)/(veig'*T*ueig);

end

