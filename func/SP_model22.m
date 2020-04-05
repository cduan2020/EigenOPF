function [ M, D, K ] = SP_model22( ps )

% 2nd order dynamics on ALL buses and generator internals, linearized
% M \ddot{\delta} + D \dot{\delta} + K \delta = 0

x_d = ps.gen_dyn(:,1); %transient reactances for generators
H = ps.gen_dyn(:,2);   %inertia constants for generators
D = ps.gen_dyn(:,3);    % damping constants for generators
baseMVA = ps.baseMVA;  %common system base

%% calculate the extended bus admittance matrix
ngi = size(ps.gen,1); % number of generators --> generator internal nodes
generator_buses = unique(ps.gen(:,1)); %note, ID == row, because internal indexing
ngt = length(generator_buses);

% Other (load) buses: buses with no generators attached
n = size(ps.bus, 1);
loads = true(n, 1);
loads(generator_buses) = false;
load_buses = ps.bus(loads, 1); %note, ID == now, because internal indexing
nl = length(load_buses);

Y0 = ps.Y; %bus admittance matrix of the original system
Y0gg = Y0(generator_buses,generator_buses);
Y0ll = Y0(load_buses,load_buses);
Y0gl = Y0(generator_buses,load_buses);
Y0lg = Y0(load_buses,generator_buses);

% Compute Y0it (admittance between generator internals and terminals)
% a generator terminal can have multiple generator internal nodes on it,
% so this matrix is not square.
% But the number of nnz's is exactly ngi. Preallocate, fast sparse

rows = 1:ngi;
[~, cols ] = ismember(ps.gen(:,1), generator_buses);
vals = -1 ./ (1j*x_d);

Y0it = sparse(rows, cols, vals, ngi, ngt);
Yd = sparse(1:ngi, 1:ngi, -sum(Y0it,2)); %negative row sum of possibly multiple transient reactances
trsum = sum(Y0it.',2);

% formulate the extended bus admittance matrix
% ordering: generator internal buses, generator terminals, loads
Y_SP = [
    Yd,             Y0it,              sparse(ngi,nl);
    Y0it.',         Y0gg-diag(trsum),  Y0gl;
    sparse(nl,ngi), Y0lg,              Y0ll];


%% Calculate generator internal voltages
% note, for buses, internal voltage = terminal voltage, and x_d=0

% P and Q injected at generator terminals in p.u. on system base MVA.
Pi = ps.gen(:,2) / baseMVA; 
Qi = ps.gen(:,3) / baseMVA;

% Voltage magnitude V and phase angle phi for the generator terminal buses
tb =  ps.gen(:, 1);
V =   ps.bus(tb, 8);
phi = ps.bus(tb, 9) / 180 * pi;

% Compute the complex voltage E at the internal nodes of generators
E = ((V + Qi .* x_d ./ V) + 1j * (Pi .* x_d ./ V)) .* exp(1j * phi);

%% calculate system matrices

% System reference frequency omega_R (in radian):
if isfield(ps, 'ref_freq')
    omega_R = 2 * pi * ps.ref_freq;
else
    omega_R = 2 * pi * 60;    
end

N = ngi + ngt + nl;

D = sparse(1:ngi, 1:ngi, D / omega_R, N, N);
M = sparse(1:ngi, 1:ngi, 2*H / omega_R, N, N);

%B = sparse(1:N, 1:N, D ./ (2*H)); % very large beta for virtual machines

%% linearization

% terminal buses: voltage angle = rotor angle (since x_d = 0 assumed)
Vg = ps.bus(generator_buses, 8);
Vl = ps.bus(load_buses,8);
phi_g = ps.bus(generator_buses,9) / 180.0 * pi;
phi_l = ps.bus(load_buses,9) / 180.0 * pi;

% concat generator internals with the virtual bus internal voltages
E = [E; Vg.*exp(1j*phi_g); Vl.*exp(1j*phi_l)];

DE = sparse(1:N, 1:N, abs(E));
K1 = DE * abs(Y_SP) * DE;  %|Ei * Ek * Yik|
gamma = angle(Y_SP) - pi/2 .* (Y_SP ~= 0);
delta = angle(E);

%% Calculate sparse Pik = -|Ei * Ek * Yik| .* cos(delta(i) - delta(k) - gamma(i,k))

% use the sparsity pattern of K1, which determines the sparsity of K
[rows, cols, k] = find(K1);
DIK = sparse(rows, cols, delta(rows) - delta(cols)); %delta(i) - delta(k)

% sparse cosine!!! otherwise cos(0) = 1 ==> full
% but use the sparsity pattern of K1!!!
arg = DIK - gamma;
vals = arg(K1~=0); %sequence the nonzeros into a vector (same order as k)

vals = -k .* cos(vals); %cos

K = sparse(rows, cols, vals);
rowsum = sum(K, 2); % make K(i,i) = negative row sum
K = K - sparse(1:N, 1:N, rowsum);


end

