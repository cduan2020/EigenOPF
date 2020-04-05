function [Ynn,Ynr,Yrn,Yrr,genbus,loadbus] = MakeExtendedYbus( ps )
% MakeExtendedYbus calculates the extended admittance matrix and returns it
% split to 4 blocks suitable for Kron-reduction. 

% Note: in ps format, one bus may have at most one generator.

%% Get constants
x_d = ps.gen_dyn(:,1); %transient reactances

%% Get indices
nbus = size(ps.bus,1);  % number of buses
ngen = size(ps.gen,1);  % number of generator internal nodes = generator terminal nodes in PS
genbus = ps.gen(:,1);   % generator buses (buses with generator attached)
%genbus = unique(ps.gen(:,1));
%ngt = length(gtb);      %number of generator terminal nodes

% Other (load) buses.
is_load = true(nbus,1);
is_load(genbus) = false;
loadbus = ps.bus(is_load, 1);
nload = size(loadbus, 1);

%% Computing Y
Y0gl = ps.Y(genbus,loadbus);
Y0lg = ps.Y(loadbus,genbus);
Yd = sparse(1:ngen, 1:ngen, 1 ./ (1j * x_d));

% Add the shunt admittances equivalent to the given loads to the diagonal
% of Y0gg and Y0ll to obatin Y0ggt and Y0llt.
Plg = ps.bus(genbus,3) / ps.baseMVA;
Qlg = ps.bus(genbus,4) / ps.baseMVA;
Vg =  ps.bus(genbus,8);
Y0ggt = ps.Y(genbus,genbus) + sparse(1:ngen, 1:ngen, Plg ./ (Vg .^ 2) - 1j * (Qlg ./ (Vg .^ 2)));

Pll = ps.bus(loadbus,3) / ps.baseMVA;
Qll = ps.bus(loadbus,4) / ps.baseMVA;
Vl  = ps.bus(loadbus,8);
Y0llt = ps.Y(loadbus,loadbus) + sparse(1:nload,  1:nload,  Pll ./ (Vl .^ 2) - 1j * (Qll ./ (Vl .^ 2)));

% Parts for Kron reduction (using the notation in our 2013 Nat Phys paper).
Ynn = Yd;
Ynr = sparse(1:ngen, 1:ngen, -1 ./ (1j * x_d), ngen, nbus);
Yrn = Ynr.';
Yrr = [...
    Y0ggt - diag(sum(Ynr(:, 1:ngen), 1)), Y0gl; ...
    Y0lg, Y0llt];

end

