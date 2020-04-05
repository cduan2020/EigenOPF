function [Asys, Bsys, Csys, Dsys]=DAEsys2(ps,V)

%% import dynamical parameters
nb=size(ps.bus,1);
ng=size(ps.gen,1);
nl=nb-ng;
baseMVA = ps.baseMVA;  %common system base
x_d = sparse(ps.gen_dyn(:,1)); %transient reactances for generators
H = ps.gen_dyn(:,2);   %inertia constants for generators
D = ps.gen_dyn(:,3);    % damping constants for generators
Ka = ps.gen_dyn(:,5);
Ks = ps.gen_dyn(:,6);
if isfield(ps, 'ref_freq')
    omega_R = 2 * pi * ps.ref_freq;
else
    omega_R = 2 * pi * 60;    
end

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
Kg=Cg+diag(1j*x_d)*Yg;


%% calculat system matrices
V=sparse(V);
volt=abs(V);
theta=angle(V);
Eg = Kg*V;

dSg_dtheta = diag(Eg)*conj(Yg)*conj(diag(1j*V))+diag(conj(Yg)*conj(V))*Kg*diag(1j*V);
dSg_dvolt = diag(Eg)*conj(Yg)*conj(diag(exp(1j*theta)))+diag(conj(Yg)*conj(V))*Kg*diag(exp(1j*theta));
dSl_dtheta = diag(Cl*V)*conj(Yl)*conj(diag(1j*V))+diag(conj(Yl)*conj(V))*Cl*diag(1j*V);
dSl_dvolt = diag(Cl*V)*conj(Yl)*conj(diag(exp(1j*theta)))+diag(conj(Yl)*conj(V))*Cl*diag(exp(1j*theta));

dPg_dE = diag(((Kg*V)./abs(Kg*V)).*(conj(Yg)*conj(V)));
Eg = Kg*V;

% H1=[ -1j*diag(Eg); zeros(nl,ng) ];
% H2=[zeros(ng,ng); zeros(nl,ng)];  
% H3=[Kg*diag(1j*V); dSl_dtheta];
% H4=[Kg*diag(exp(1j*theta)); dSl_dvolt];  

H1=[ -1j*diag(Eg); zeros(nl,ng) ];
H2=[-diag(((Eg)./abs(Eg)).*Ka.*Ks/omega_R); zeros(nl,ng)];   % change to take Ka Ks
H3=[Kg*diag(1j*V); dSl_dtheta];
H4=[Kg*diag(exp(1j*theta))+diag(((Eg)./abs(Eg)).*Ka)*Cg; dSl_dvolt];  % change to take Ka Ks


A11=sparse(ng,ng);
A12=sparse(1:ng,1:ng,ones(1,ng),ng,ng);
A21=sparse(ng,ng);
A22=-invM*dPg_dE*diag(Ka.*Ks)/omega_R;  % change to Ka Ks generated D

B11=sparse(ng,nb);
B12=sparse(ng,nb);
B21=-invM*real(dSg_dtheta);
B22=-invM*(real(dSg_dvolt)-dPg_dE*diag(Ka)*Cg); % change to Ka Ks generated D


C11=real(H1);
C12=real(H2);
C21=imag(H1);
C22=imag(H2);

D11=real(H3);
D12=real(H4);
D21=imag(H3);
D22=imag(H4);

Asys=[A11 A12; A21 A22];
Bsys=[B11 B12; B21 B22];
Csys=[C11 C12; C21 C22];
Dsys=[D11 D12; D21 D22];

end












