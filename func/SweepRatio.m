clear;clc;


%========================================================================
% system='test_system_3gen';
% tbd0=0.8;
%========================================================================
% system='test_system_10gen';
% tbd0=0.6;
%========================================================================
system='test_system_50gen';
tbd0=0.2;
%========================================================================
% system='test_system_505gen';
% tbd0=0.8;
%========================================================================
% system='test_system_1445gen';
% tbd0=0.8;
%========================================================================
% system='case_ACTIVSg2000';
% tbd0=1.0;
%========================================================================
loadmode=0;
mpc = load_system(system, loadmode);
ps = mpc2ps(mpc);
nb=size(ps.bus,1);
ng=size(ps.gen,1);

ps.gen(:,10)=0;


TableFinal=zeros(41,1);
for jjj=0:40


%==================== Set the controllable load ==========================
disploc=find(ps.bus(:,2)==1 & ps.bus(:,3)>0);   % locations of controllable load
nd=length(disploc);     % number of controllable load
dispratio=jjj*0.02*ones(nd,1);   % percentage of controllable load
%=========================================================================
 


maxReal=[];
ActivePower=[];
Mset.x=[];
Mset.num=0;
Mset.abassia=[];
Mset.daba_da=[];
Mset.daba_dv=[];


%==================== Get the eig for nomical power flow===================

% [results,success] = runpf(ps);
% V0=results.bus(:,8).*exp(1j*results.bus(:,9)/180*pi);
% results.x=[results.bus(:,9)/180*pi ; results.bus(:,8)];
% [Asys0, Bsys0, Csys0, Dsys0]=DAEsys(ps,V0);
% Afull=Asys0-Bsys0*(Dsys0\Csys0);
% [Ueig0,D,Veig0] = eig(full(Afull));
% lambda=diag(D);
% maxReal=[maxReal max(real(lambda(abs(lambda)>10^-6)))];


%==================== Set up the optimization problem =====================
mpc_disp=load2control(ps,[],disploc,0,dispratio);
theta_IND=1:nb;
volt_IND=nb+1:2*nb;
Pg_IND=2*nb+1:2*nb+ng+nd;
Qg_IND=2*nb+ng+nd+1:2*nb+2*(ng+nd);
gamma_IND=2*nb+2*(ng+nd)+1;
Ctheta=sparse(1:nb,theta_IND,ones(nb,1),nb,2*nb+2*(ng+nd)+1);
Cvolt=sparse(1:nb,volt_IND,ones(nb,1),nb,2*nb+2*(ng+nd)+1);
Cgamma=sparse(1,gamma_IND,1,1,2*nb+2*(ng+nd)+1);


mpc_disp.gencost(:,5:7)=mpc_disp.gencost(:,5:7)*0;
mpopt = mpoption;
mpopt.opf.ac.solver='ipopt';
om = opf_setup(mpc_disp, mpopt);

% [results, success, raw] = opf_execute(om, mpopt); 
% V0=results.bus(:,8).*exp(1j*results.bus(:,9)/180*pi);   

%=========================================================================
if jjj==0
% rng('shuffle')
% cutratio=rand(nd,1);
cutratio=zeros(nd,1);
ps_new=ps;
ps_new.bus(disploc,3)=ps_new.bus(disploc,3).*(1-cutratio.*dispratio);
ps_new.bus(disploc,4)=ps_new.bus(disploc,4).*(1-cutratio.*dispratio);
ps_new.gen(:,2)=ps_new.gen(:,2)*sum(ps_new.bus(:,3))/sum(ps.bus(:,3));
[results,success] = runpf(ps_new); ps.bus=results.bus; ps.gen=results.gen; ps.branch=results.branch;
V0=results.bus(:,8).*exp(1j*results.bus(:,9)/180*pi);
end
%=========================================================================


% [lmax, eigvals] = EN_Lmax( ps );
% results.x=[results.bus(:,9)/180*pi ; results.bus(:,8)];

[Asys0, Bsys0, Csys0, Dsys0]=DAEsys(ps,V0);
Afull=Asys0-Bsys0*(Dsys0\Csys0);
[Ueig0,D,Veig0] = eig(full(Afull));

n=size(Afull,1)/2;
P=-full(Afull(n+1:2*n,1:n));
[UP,DP,VP]=eig(P);

% [Ueig0,D,Veig0] = speigs(Afull,2,'largestreal');

lambda=diag(D);
% [freq,damping_ratio,pfac,mode2gen]=ModalAnalysis(lambda,Ueig0,Veig0);
maxReal=[maxReal max(real(lambda(abs(lambda)>10^-6)))];


plot(lambda, '*');
grid on;
hold on;
xlim([-0.08 0]);
ylim([-20 20]);
xlabel 'Real part'
ylabel 'Imaginary part'
lamIndex=find(abs(lambda)>10^-6 & imag(lambda)>=0 & real(lambda)>min(maxk(real(lambda),20))-0.0);



Q=sparse(2*nb+2*(ng+nd)+1,2*nb+2*(ng+nd)+1);
c=[sparse(2*nb+2*(ng+nd),1); 1];
om = add_vars(om, 'gamma', 1, [], -100, 100);
om = add_costs(om, 'usr', struct('H', Q, 'Cw', c));
om = build_cost_params(om);


iter=1;
fobj=maxReal(1);
tbd=tbd0;
tol=10^-4;
tic;
while (iter<=100)
    
    Aadd=[Ctheta; Cvolt];
    bup=[ones(nb,1)*tbd+angle(V0); ones(nb,1)*tbd+abs(V0)];
    bdn=[-ones(nb,1)*tbd+angle(V0); -ones(nb,1)*tbd+abs(V0)];
    
    for k=1:length(lamIndex)
        num_of_lam=lamIndex(k);
        ueig=[Ueig0(:,num_of_lam);-Dsys0\(Csys0*Ueig0(:,num_of_lam))];
        veig=[Veig0(:,num_of_lam);-(Dsys0')\(Bsys0'*Veig0(:,num_of_lam))];
        [dsys_da, dsys_dv]=SysGradient(ps,V0,ueig,veig);
        Aadd=[Aadd; real(dsys_da)*Ctheta+real(dsys_dv)*Cvolt-Cgamma];
        bup=[bup; real(dsys_da)*angle(V0)+real(dsys_dv)*abs(V0)-real(lambda(num_of_lam))];
        bdn=[bdn; -Inf];
    end
    
    for k=1:Mset.num
        if norm([angle(V0); abs(V0)]-Mset.x(:,k),inf)<=tbd
            Aadd=[Aadd; Mset.daba_da(k,:)*Ctheta+Mset.daba_dv(k,:)*Cvolt-Cgamma];
            bup=[bup; [Mset.daba_da(k,:) Mset.daba_dv(k,:)]*Mset.x(:,k)-Mset.abassia(k)];
            bdn=[bdn; -Inf];
        end
    end
    
    
    om_new=om;
    om_new = add_constraints(om_new, 'EigCon', Aadd, bdn, bup);
    
    
%     [results, success, raw] = opf_execute(om_new, mpopt);
    
    
    mpopt.verbose=1;
    mpopt.opf.init_from_mpc=-1;
%     [results, success, raw] = ktropf_solver_eigen(om_new, mpopt, V0);
    
    [results, success, raw] = ipoptopf_solver_eigen(om_new, mpopt, V0);
    
%     V1=(Cvolt*results.x).*exp(1j*Ctheta*results.x);
    
    V1=results.bus(:,8).*exp(1j*results.bus(:,9)/180*pi);  
    
    [Asys1, Bsys1, Csys1, Dsys1]=DAEsys(ps,V1);
    Afull=Asys1-Bsys1*(Dsys1\Csys1);
    [Ueig1,D,Veig1] = eig(full(Afull));
    lambda1=diag(D);
    
    lamIndex1=find(abs(lambda1)>10^-6 & imag(lambda1)>=0 & real(lambda1)>min(maxk(real(lambda1),20))-0.0);
    
    Inx1=find(abs(lambda1)>10^-6 & imag(lambda1)>=0);
    [abassia,Inx2]=max(real(lambda1(Inx1)));
    maxInx=Inx1(Inx2);
    
    
    if abassia<fobj && success==1
        V0=V1;
        tbd=min(2*tbd,tbd0);
        fobj=abassia;
        Ueig0=Ueig1;
        Veig0=Veig1;
        Asys0=Asys1;
        Bsys0=Bsys1;
        Csys0=Csys1;
        Dsys0=Dsys1;
        maxReal=[maxReal abassia];
        ActivePower=[ActivePower results.x(Pg_IND)];
        lambda=lambda1;
        lamIndex=lamIndex1;
%         plot(lambda, '*');
%         grid on;
%         hold on;
%         xlim([-0.08 0]);
%         ylim([-20 20]);
%         xlabel 'Real part'
%         ylabel 'Imaginary part'
    else
        Mset.num=Mset.num+1;
        Mset.x=[Mset.x results.x(1:2*nb)];
        Mset.abassia=[Mset.abassia abassia];
        ueig=[Ueig1(:,maxInx);-(Dsys1)\(Csys1*Ueig1(:,maxInx))];
        veig=[Veig1(:,maxInx);-(Dsys1')\(Bsys1'*Veig1(:,maxInx))];
        [dsys_da, dsys_dv]=SysGradient(ps,V1,ueig,veig);
        Mset.daba_da=[Mset.daba_da ; real(dsys_da)];
        Mset.daba_dv=[Mset.daba_dv ; real(dsys_dv)];
        
        [t, dxnorm, Mset]=BackLineSearchMset(ps,V0,V1,fobj,tol,Mset);
        log2(1/t)
        
        
        tbd=t*dxnorm;
        if tbd<tol
            break;
        end
    end
    iter=iter+1;
end
toc

maxReal
improve=(maxReal(1)-maxReal(end))/abs(maxReal(1))

results.x(Pg_IND(ng+1:ng+nd))

beta_max= max(ps.gen_dyn(:,3)./(4*ps.gen_dyn(:,2)))

beta_min= min(ps.gen_dyn(:,3)./(4*ps.gen_dyn(:,2)))

TableFinal(jjj+1)=maxReal(end);

end


% (maxReal(5)-maxReal(1))/((maxReal(end)-maxReal(1)))
% maxReal(5)

plot(lambda, '*');
grid on;
hold on;
xlim([-0.08 0]);
ylim([-20 20]);
xlabel 'Real part'
ylabel 'Imaginary part'
