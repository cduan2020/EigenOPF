clear;clc;

addpath(genpath(cd));

% choose the test system
system='test_system_10gen';
% system='test_system_16gen';
% system='test_system_50gen';

% specify the number of load controllability levels
% dispratio: relative increase or decrease allowed for each controllable load
N=11;
maxReal=zeros(1,N+1);
improve=zeros(1,N+1);
best_id=zeros(1,N+1);
Vps1=cell(3,N+1);

% specify the load level
loadability = 1;


%=========================================================
% initialize the power flow solution
loadmode=0;
ps = load_system_ps(system,loadmode);
ps = mpc2ps(ps);
ps.bus(ps.bus(:,8) > ps.bus(:,12),12) = 1.2;
ps.bus(ps.bus(:,8) < ps.bus(:,13),13) = 0.8;
ps.gen = ps.gen(:,1:21);
ps.bus(:,3:4)=ps.bus(:,3:4)*loadability;
ps.gen(:,2)=ps.gen(:,2)*loadability;
[results,success] = runpf(ps); ps.bus=results.bus; ps.gen=results.gen; ps.branch=results.branch;
Vps0=results.bus(:,8).*exp(1j*results.bus(:,9)/180*pi);
%=========================================================

% calculate the initial linearized system and damping ratio and plot the system eigenvalues
[Asys0, Bsys0, Csys0, Dsys0]=DAEsys(ps,Vps0);
Afull=Asys0-Bsys0*(Dsys0\Csys0);
[Ueig0,D,Veig0] = eig(full(Afull));
lambda=diag(D);
figure(1);
hold off;
y=-12:.01:12;
x=-.03*abs(y);
plot(x,y, '-.');
hold on;
grid on;
plot(lambda, '*');
grid on;
hold on;
xlim([-5 5]);
ylim([-25 25]);
xlabel 'Real part'
ylabel 'Imaginary part'

% gradually increase the ratio of controllable load and obtain optimized power flow solution
dispratio=0.0;
load_level=1;
maxReal(1,1)=max(real(lambda(abs(lambda)>10^-6 & imag(lambda)>0.01))./(abs(imag(lambda(abs(lambda)>10^-6 & imag(lambda)>0.01)))));
Vps1{1,1}=Vps0;
for i=1:N
    i
    [Vps1{1,i+1},success,maxReal(1,i+1),improve(1,i+1),best_id(1,i+1)]=OptimizeRatio(ps,Vps0,0.01*(i-1),load_level);
    Vps0 = Vps1{1,i+1};
    maxReal
end

% calculate the optimized linearized system and plot the system eigenvalues
figure(1);
[Asys0, Bsys0, Csys0, Dsys0]=DAEsys(ps,Vps0);
Afull=Asys0-Bsys0*(Dsys0\Csys0);
[Ueig0,D,Veig0] = eig(full(Afull));
lambda=diag(D);
figure(1);
plot(lambda, '*');
grid on;
hold on;
xlim([-5 1]);
ylim([-15 15]);
xlabel 'Real part'
ylabel 'Imaginary part'


