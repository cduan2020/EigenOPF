function [ps, success] = runpf_ps( ps, mpopt )
%runpf_ps Runs AC power flow exactly as MATPOWER, except using ps format.
% Matpower options used: 
% - verbosity is set to zero
% - force Q limits -> fail if unfeasible

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;


%% record the original bus types
origBusTypes = ps.bus(:,2);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(ps.bus, ps.gen);

%% generator info
on = find(ps.gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = ps.gen(on, GEN_BUS);                %% what buses are they at?

%% initial state
V0  = ps.bus(:, VM) .* exp(sqrt(-1) * pi/180 * ps.bus(:, VA));
vcb = ones(size(V0));           %% create mask of voltage-controlled buses
vcb(pq) = 0;                    %% exclude PQ buses
k = find(vcb(gbus));            %% in-service gens at v-c buses
V0(gbus(k)) = ps.gen(on(k), VG) ./ abs(V0(gbus(k))).* V0(gbus(k));

%% prepare for Qlimit forcing iteration
ref0 = ref;                         %% save index and angle of
Varef0 = ps.bus(ref0, VA);          %% original reference bus(es)
limited = [];                       %% list of indices of gens @ Q lims
fixedQg = zeros(size(ps.gen, 1), 1);   %% Qg of gens at Q limits

if (nargin < 2)
    mpopt = mpoption;
    mpopt.verbose = 0;
end

% Qlimit enforcing loop
repeat = 1;
while repeat    
    %% compute complex bus power injections (generation minus load)
    Sbus = makeSbus(ps.baseMVA, ps.bus, ps.gen);

    %% run the power flow (default tolerances)
    [V, success] = newtonpf(ps.Y, Sbus, V0, ref, pv, pq, mpopt);

    %% ----- update data ----- (from pfsoln)
    
    %------ bus voltages
    ps.bus(:, VM) = abs(V);
    ps.bus(:, VA) = angle(V) * 180 / pi;

    %----- update Qg for gens at PV/slack buses and Pg for slack bus(es) -----    
    on = find(ps.gen(:, GEN_STATUS) > 0 & ...  %% which generators are on?
            ps.bus(ps.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and not at PQ buses
    off = find(ps.gen(:, GEN_STATUS) <= 0);    %% which generators are off?
    gbus = ps.gen(on, GEN_BUS);                %% what buses are they at?

    % compute total injected bus powers
    Sbus = V(gbus) .* conj(ps.Y(gbus, :) * V);

    % update Qg for generators at PV/slack buses
    ps.gen(off, QG) = zeros(length(off), 1);   %% zero out off-line Qg

    % don't touch the ones at PQ buses
    ps.gen(on, QG) = imag(Sbus) * ps.baseMVA + ps.bus(gbus, QD); %% inj Q + local Qd
    
    % ... at this point any buses with more than one generator will have
    % the total Q dispatch for the bus assigned to each generator. This
    % must be split between them. We do it first equally, then in proportion
    % to the reactive range of the generator.

    if length(on) > 1
        % build connection matrix, element i, j is 1 if gen on(i) at bus j is ON
        nb = size(ps.bus, 1);
        ngon = size(on, 1);
        Cg = sparse((1:ngon)', gbus, ones(ngon, 1), ngon, nb);

        % divide Qg by number of generators at the bus to distribute equally
        ngg = Cg * sum(Cg)';    %% ngon x 1, number of gens at this gen's bus
        ps.gen(on, QG) = ps.gen(on, QG) ./ ngg;

        % divide proportionally
        Cmin = sparse((1:ngon)', gbus, ps.gen(on, QMIN), ngon, nb);
        Cmax = sparse((1:ngon)', gbus, ps.gen(on, QMAX), ngon, nb);
        Qg_tot = Cg' * ps.gen(on, QG);     %% nb x 1 vector of total Qg at each bus
        Qg_min = sum(Cmin)';            %% nb x 1 vector of min total Qg at each bus
        Qg_max = sum(Cmax)';            %% nb x 1 vector of max total Qg at each bus
        ig = find(Cg * Qg_min == Cg * Qg_max);  %% gens at buses with Qg range = 0
        Qg_save = ps.gen(on(ig), QG);
        ps.gen(on, QG) = ps.gen(on, QMIN) + ...
            (Cg * ((Qg_tot - Qg_min)./(Qg_max - Qg_min + eps))) .* ...
                (ps.gen(on, QMAX) - ps.gen(on, QMIN));    %%    ^ avoid div by 0
        ps.gen(on(ig), QG) = Qg_save;
    end                                             %% (terms are mult by 0 anyway)

    % update Pg for slack gen(s)
    for k = 1:length(ref)
        refgen = find(gbus == ref(k));              %% which is(are) the reference gen(s)?
        ps.gen(on(refgen(1)), PG) = real(Sbus(refgen(1))) * ps.baseMVA ...
                                + ps.bus(ref(k), PD);  %% inj P + local Pd
        if length(refgen) > 1       %% more than one generator at this ref bus
            % subtract off what is generated by other gens at this bus
            ps.gen(on(refgen(1)), PG) = ps.gen(on(refgen(1)), PG) ...
                                    - sum(ps.gen(on(refgen(2:length(refgen))), PG));
        end
    end
    
    %% check for Q limit violations
    mx = find( ps.gen(:, GEN_STATUS) > 0 ...
            & ps.gen(:, QG) > ps.gen(:, QMAX) );
    mn = find( ps.gen(:, GEN_STATUS) > 0 ...
            & ps.gen(:, QG) < ps.gen(:, QMIN) );
        
    if ~isempty(mx) || ~isempty(mn)  % we have some Q limit violations
        % first check for INFEASIBILITY (all remaining gens violating)
        infeas = union(mx', mn')';  %% transposes handle fact that union of scalars is a row vector
        remaining = find( ps.gen(:, GEN_STATUS) > 0 & ...
                        ( ps.bus(ps.gen(:, GEN_BUS), BUS_TYPE) == PV | ...
                          ps.bus(ps.gen(:, GEN_BUS), BUS_TYPE) == REF ));
                      
        if length(infeas) == length(remaining) && all(infeas == remaining)
            %all generators violating == infeasible
            success = 0;
            break;
        end

        % save corresponding limit values
        fixedQg(mx) = ps.gen(mx, QMAX);
        fixedQg(mn) = ps.gen(mn, QMIN);
        mx = [mx;mn];

        % convert to PQ bus
        ps.gen(mx, QG) = fixedQg(mx);      %% set Qg to binding limit
        ps.gen(mx, GEN_STATUS) = 0;        %% temporarily turn off gen,
        for i = 1:length(mx)            %% (one at a time, since
            bi = ps.gen(mx(i), GEN_BUS);   %%  they may be at same bus)
            ps.bus(bi, [PD,QD]) = ...      %% adjust load accordingly,
                ps.bus(bi, [PD,QD]) - ps.gen(mx(i), [PG,QG]);
        end
        
        if length(ref) > 1 && any(ps.bus(ps.gen(mx, GEN_BUS), BUS_TYPE) == REF)
            error('Sorry, MATPOWER cannot enforce Q limits for slack buses in systems with multiple slacks.');
        end
        
        ps.bus(ps.gen(mx, GEN_BUS), BUS_TYPE) = PQ;   %% & set bus type to PQ

        % update bus index lists of each type of bus
        ref_temp = ref;
        [ref, pv, pq] = bustypes(ps.bus, ps.gen);
        
        % previous line can modify lists to select new REF bus
        % if there was none, so we should update bus with these
        % just to keep them consistent
        if ref ~= ref_temp
            ps.bus(ref, BUS_TYPE) = REF;
            ps.bus( pv, BUS_TYPE) = PV;
            %fprintf('Bus %d is new slack bus\n', ref);
        end
        limited = [limited; mx];
    else
        repeat = 0; % no more generator Q limits violated
    end
end %while repeat Qlimit

%% restore injections from limited gens (those at Q limits), if any
if ~isempty(limited)

    ps.gen(limited, QG) = fixedQg(limited);    % restore Qg value,
    for i = 1:length(limited)                  % (one at a time, since
        bi = ps.gen(limited(i), GEN_BUS);      % they may be at same bus)
        ps.bus(bi, [PD,QD]) = ...              % re-adjust load,
            ps.bus(bi, [PD,QD]) + ps.gen(limited(i), [PG,QG]);
    end
    ps.gen(limited, GEN_STATUS) = 1;           % and turn gen back on.
    if ref ~= ref0
        % adjust voltage angles to make original ref bus correct
        ps.bus(:, VA) = ps.bus(:, VA) - ps.bus(ref0, VA) + Varef0;
    end
end

%% restore bus types
ps.bus(:,2) = origBusTypes;

end

