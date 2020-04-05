function [Asys, Bsys, Csys, Dsys]=DAEsys(ps,V)

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

%% import network parameters

gbus=ps.gen(:,1);
lbus=setdiff(1:nb,gbus);
Ybus=ps.Y;
Yg=Ybus(gbus,:);
Yl=Ybus(lbus,:);
Cg=sparse(1:ng,gbus,ones(ng,1),ng,nb);
Cl=sparse(1:nl,lbus,ones(nl,1),nl,nb);
Kg=Cg+sparsediag(1j*Xq)*Yg;
Kgp=Cg+sparsediag(1j*Xqp)*Yg;

%% calculat system matrices
V=sparse(V);
volt=abs(V);
theta=angle(V);
Eg = Kg*V;
Egp = Kgp*V;
cosdelta = cos(angle(Eg./abs(Eg)));
sindelta = sin(angle(Eg./abs(Eg)));

dSg_dtheta = sparsediag(Egp)*conj(Yg)*conj(sparsediag(1j*V))+sparsediag(conj(Yg)*conj(V))*Kgp*sparsediag(1j*V);
dSg_dvolt = sparsediag(Egp)*conj(Yg)*conj(sparsediag(exp(1j*theta)))+sparsediag(conj(Yg)*conj(V))*Kgp*sparsediag(exp(1j*theta));
dSl_dtheta = sparsediag(Cl*V)*conj(Yl)*conj(sparsediag(1j*V))+sparsediag(conj(Yl)*conj(V))*Cl*sparsediag(1j*V);
dSl_dvolt = sparsediag(Cl*V)*conj(Yl)*conj(sparsediag(exp(1j*theta)))+sparsediag(conj(Yl)*conj(V))*Cl*sparsediag(exp(1j*theta));

xy2dq = 1j*conj(Eg./abs(Eg));
dIg_delta = sparsediag(Yg*V)*sparsediag(xy2dq)*(-1j);
dIg_theta = sparsediag(xy2dq)*Yg*sparsediag((1j*V));
dIg_volt = sparsediag(xy2dq)*Yg*sparsediag(exp(1j*theta));

% dIg_delta =dIg_delta*0;
% dIg_theta =dIg_theta*0;
% dIg_volt =dIg_volt*0;

% dPg_dE = diag(((Kg*V)./abs(Kg*V)).*(conj(Yg)*conj(V)));


Asys=[sparse(ng,ng) speye(ng) sparse(ng,ng) sparse(ng,ng);
    sparse(ng,ng) -invM*D  sparse(ng,ng) sparse(ng,ng);
    -sparsediag((Xd-Xdp)./Td0p)*real(dIg_delta) sparsediag(Ka.*Ks./Td0p) -sparsediag(1./Td0p) sparse(ng,ng);
    +sparsediag((Xq-Xqp)./Tq0p)*imag(dIg_delta) sparse(ng,ng) sparse(ng,ng) -sparsediag(1./Tq0p)];

Bsys=[sparse(ng,nb) sparse(ng,nb);
    -invM*real(dSg_dtheta) -invM*(real(dSg_dvolt));
    -sparsediag((Xd-Xdp)./Td0p)*real(dIg_theta) -sparsediag((Xd-Xdp)./Td0p)*real(dIg_volt)-sparsediag(Ka./Td0p)*Cg;
    sparsediag((Xq-Xqp)./Tq0p)*imag(dIg_theta) sparsediag((Xq-Xqp)./Tq0p)*imag(dIg_volt)];

% Csys=[sparsediag(real(xy2dq.*Egp).*cosdelta)-sparsediag(imag(xy2dq.*Egp).*sindelta) sparse(ng,ng) sparsediag(cosdelta) sparsediag(sindelta);
%     sparse(nl,ng) sparse(nl,ng) sparse(nl,ng) sparse(nl,ng);
%     sparsediag(imag(xy2dq.*Egp).*cosdelta)+sparsediag(real(xy2dq.*Egp).*sindelta) sparse(ng,ng) sparsediag(sindelta) -sparsediag(cosdelta);
%     sparse(nl,ng) sparse(nl,ng) sparse(nl,ng) sparse(nl,ng)];

Csys=[sparsediag(real(1j*Egp)) sparse(ng,ng) sparsediag(cosdelta) sparsediag(sindelta);
    sparse(nl,ng) sparse(nl,ng) sparse(nl,ng) sparse(nl,ng);
    sparsediag(imag(1j*Egp)) sparse(ng,ng) sparsediag(sindelta) -sparsediag(cosdelta);
    sparse(nl,ng) sparse(nl,ng) sparse(nl,ng) sparse(nl,ng)];

Dsys=[-real(Kgp*sparsediag(1j*V)) -real(Kgp*sparsediag(exp(1j*theta)));
    real(dSl_dtheta) real(dSl_dvolt);
    -imag(Kgp*sparsediag(1j*V)) -imag(Kgp*sparsediag(exp(1j*theta)));
    imag(dSl_dtheta) imag(dSl_dvolt)];


end













