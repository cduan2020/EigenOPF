function [Tm, Sd, Efd] = result2contrl(ps,Vps)

%% import dynamical parameters
nb=size(ps.bus,1);
ng=size(ps.gen,1);
nl=nb-ng;
baseMVA = ps.baseMVA;  %common system base
mBase = ps.gen_dyn(:,4); % machine base
x_d = sparse(ps.gen_dyn(:,1)); %transient reactances for generators

%% import network parameters
gbus=ps.gen(:,1);
lbus=setdiff(1:nb,gbus);
Ybus=ps.Y;
Yg=Ybus(gbus,:);
Yl=Ybus(lbus,:);
Cg=sparse(1:ng,gbus,ones(ng,1),ng,nb);
Cl=sparse(1:nl,lbus,ones(nl,1),nl,nb);
Kg=Cg+diag(1j*x_d)*Yg;

%% Output
Pm = real(diag(Kg*Vps)*conj(Yg)*conj(Vps));
Tm = Pm*baseMVA./mBase; % torque in p.u. machine base
Sd = Cl'*(-diag(Cl*Vps)*conj(Yl)*conj(Vps));
Sd = Sd*ps.baseMVA; % loads in MVA system base
Efd = abs(Kg*Vps);

end
