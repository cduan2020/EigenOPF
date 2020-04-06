function mpc=ps2mpc(ps)

% Originally created by Ferenc Molnar.

% mpc = 
% 
%   struct with fields:
% 
%     ref_freq: 60
%      gen_dyn: [3×4 double]
%      version: '2'
%      baseMVA: 100
%          bus: [9×13 double]
%          gen: [3×21 double]
%       branch: [9×13 double]
%      gencost: [3×7 double]

% ps = 
% 
%   struct with fields:
% 
%          gen: [3×21 double]
%      gen_dyn: [3×4 double]
%          bus: [9×13 double]
%      baseMVA: 100
%            Y: [9×9 double]
%     ref_freq: 60
%      version: '2'
%       branch: [9×13 double]
%      gencost: [3×7 double]

[branch, Gs, Bs]=Y2Branch(ps.Y);
mpc=ps;
mpc.branch=branch;
mpc.bus(:,5)=Gs*mpc.baseMVA;
mpc.bus(:,6)=Bs*mpc.baseMVA;
ng=size(mpc.gen,1);
mpc.gencost=ones(ng,1)*[2 0 0 3 0 0 0];

end
