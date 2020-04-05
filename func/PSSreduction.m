function [K,T] = PSSreduction(Ka,Tr,Ta,Tb,Tc,Ks,T1,T2,T3,T4)

% sys = zpk([0 -1/T1 -1/T3 -1/Tc],[-1/Tw -1/T2 -1/T4 -1/Ta -1/Tb -1/Td0p],Ks*Ka);


sys1 = zpk([-1/T1 -1/T3 -1/Tc],[-1/T2 -1/T4 -1/Ta -1/Tb],Ks*Ka*T1/T2*T3/T4/Ta*Tc/Tb);

sys2 = zpk([-1/Tc],[-1/Tr -1/Ta -1/Tb],1/Tr*Ka/Ta*Tc/Tb);

bopt = balredOptions('StateElimMethod','Truncate','FreqIntervals',[0.3*pi,2*pi]);

rsys1 = balred(sys1,1,bopt);

rsys2 = balred(sys2,1,bopt);

end