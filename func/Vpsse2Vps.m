function [Vps, Ypsse]=Vpsse2Vps(Vpsse, ps, nb_mpc, nb_ps)

% Vpsse: complex bus voltages in psse format
% ps: ps data structure
% nb_mpc: number of buses after converting three winding generators
% nb_ps: number of buses after adding buses for zero-load generators

% Vps: complex bus voltages in ps format
% Ypsse: psse format admittance matrix constructed from ps data

Ybus=ps.Y;
nb_psse=length(Vpsse);
nadd2=nb_ps-nb_mpc;

starbusind=nb_psse+1:nb_mpc;
oldgenind=ps.branch(end-nadd2+1:end,1).';
Iinj=zeros(nb_ps,1);
Iinj(oldgenind)=-conj((ps.bus(oldgenind,3)+1j*ps.bus(oldgenind,4))/ps.baseMVA./Vpsse(oldgenind));


Y21=Ybus([starbusind oldgenind],1:nb_psse);
Y22=Ybus([starbusind oldgenind], nb_psse+1:nb_ps);

Vadd=Y22\(Iinj([starbusind oldgenind])-Y21*Vpsse);

Vps=[Vpsse; Vadd];


if nargout>1
    Y11=Ybus(1:nb_psse,1:nb_psse);
    Y12=Ybus(1:nb_psse,[starbusind oldgenind]);
    Ypsse=Y11-Y12*(Y22\Y21);
end


end

%% tests
% V0=ps.bus(:,8).*exp(1j*ps.bus(:,9)/180*pi);
% nb_psse=39;
% nb_mpc=40;
% nb_ps=40;
% Vpsse=V0(1:nb_psse);
% Vps=Vpsse2Vps(Vpsse, ps, nb_mpc, nb_ps);
% abs(Vps)

