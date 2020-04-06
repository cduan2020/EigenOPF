function mpc = load_system_ps( casename, convert )

% Originnal created by Ferenc Molnar.

%load_system Loads the specified powergrid sample by case name and
%conversion
coder.extrinsic('eval');

if (convert == 0)       % matpower func
    mpc = eval(casename);
elseif (convert == 1)   %matpower from .mat file
    load(casename);
elseif (convert == 2)   %pst func
    [bus, line, mac_con] = eval(casename);
    fprintf('Using OLD pst conversion (adding fictitious generators)\n');
    mpc = pst2matpower_old(bus, line, mac_con);        
elseif (convert == 3)   %US power grid sample
    fname = [casename '.mat'];
    if exist(fname, 'file')
        load(fname);
    else
        %load nametable
        names = readtable('nametable.txt', 'Delimiter', '\t', 'ReadVariableNames', false);
        [rows, ~] = size(names);
        for i = 1 : rows
            id = cell2mat(table2array(names(i,1)));
            fname = cell2mat(table2array(names(i,2)));

            if (id == casename)
                load(fname); %contains mpc
                return
            end
        end
    end
end

mpc=ext2int(mpc);

% if there are loads at the generator bus, convert it into a pure load bus
% and add a new generator bus connecting the original one so that the
% generator bus always do not have load
for i=1:size(mpc.gen,1)
    if (mpc.bus(mpc.gen(i,1),3)~=0) || (mpc.bus(mpc.gen(i,1),4)~=0)
        newbus=mpc.bus(mpc.gen(i,1),:);
        newbus(1)=size(mpc.bus,1)+1;
        newbus(3)=0; newbus(4)=0;newbus(5)=0; newbus(6)=0;
        mpc.bus = [mpc.bus ; newbus];
        mpc.bus(mpc.gen(i,1),2)=1;   
        newbrh=[mpc.gen(i,1) newbus(1) 0 0.00001 0 Inf Inf Inf 0 0 1 -360 360];
        mpc.branch=[mpc.branch; newbrh];
        mpc.gen(i,1)=newbus(1);
    end
end

end

