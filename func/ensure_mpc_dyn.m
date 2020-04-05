function mpc = ensure_mpc_dyn( mpc )
%ensure_mpc_dyn Ensures that mpc has proper gen_dyn data. 

% Number of generators:
ngi = size(mpc.gen, 1);
nb = size(mpc.bus, 1);

% preallocate gen_dyn if not present
if ~isfield(mpc,'gen_dyn')        
    mpc.gen_dyn = nan(ngi, 3); %Xd, H, D for generators
end

if ~isfield(mpc,'bus_dyn')
    mpc.bus_dyn = nan(nb, 1); % D for buses
end

%% Estimating missing dynamic parameters

default_max_xd = 1;

% If H_i is not given, use the default method to estimate it.
i = isnan(mpc.gen_dyn(:,2));
mpc.gen_dyn(i,2) = default_H(abs(mpc.gen(i,2)));

% If D_i is not given (for generators), use the default method to estimate it. 
i = isnan(mpc.gen_dyn(:,3));
mpc.gen_dyn(i,3) = default_D(abs(mpc.gen(i,2)));

% If x'_{d,i} is not given, use the default method to estimate it.
i = isnan(mpc.gen_dyn(:,1));    
mpc.gen_dyn(i,1) = default_x_d(abs(mpc.gen(i,2)), default_max_xd);

% If D_i is not given (for buses), use the default method to estimate it.
i = isnan(mpc.bus_dyn(:,1));
mpc.bus_dyn(i,1) = default_D2(abs(mpc.bus(i,3)));

end

function x_d = default_x_d(P, x_d_max)
% The default method for estimating the d-axis transient reactances x'_d.
% Assign the values accoding to the empirical relation used in A.E. Motter,
% S.A. Myers, M. Anghel, & T. Nishikawa, Spontaneous synchrony in
% power-grid networks, Nat Phys 9, 191 (2013).  Impose a maximum value.

% x_d_max = 1;
% x_d_max = 2.2;

%x_d = min([92.8*P.^(-1.3), x_d_max*ones(size(P))], [], 2);
x_d = min([102.30253 * P.^(-1.35742), x_d_max*ones(size(P))], [], 2);
%disp(x_d);
end


function H = default_H(P)
% The default method for estimating the inertia constant H_i.  Assign the
% values accoding to the empirical relation used in A.E. Motter, S.A.
% Myers, M. Anghel, & T. Nishikawa, Spontaneous synchrony in
% power-grid networks, Nat Phys 9, 191 (2013).  Impose a minimum value.

H_min = 0.1;
H = max([0.04*P, H_min*ones(size(P))], [], 2);
end

function D = default_D(P)
% The default method for estimating the damping coefficient D_i. Setting
% D_i = 50 (which is what we did in the 2013 Nat Phys paper) is equivalent
% to setting
%
%   D_m = 0 * P_R/omega_R
%   D_e = 0 * P_R/omega_R
%   R = 0.02 * omega_R/P_R
%
% in the notation of the New J Phys paper.

D = 50*ones(size(P));
end


function D = default_D2(P)
    D = 50 * ones(size(P));
end


