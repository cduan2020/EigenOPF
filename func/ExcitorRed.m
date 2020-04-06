function [Ka_red,Ks_red] = ExcitorRed(Ka,Tr,Ta,Tb,Tc,Ks,T1,T2,T3,T4)

s=1.5*pi*1j;
Tw=10*ones(size(Ka,1),1);



Ka_red = abs((1./(1+s*Tr)).*(Ka./(1+s*Ta)).*((1+s*Tc)./(1+s*Tb)));


Ks_red = abs(Ks.*1./(1+s*Tr).*((s*Tw)./(1+s*Tw)).*((1+s*T1)./(1+s*T2)).*((1+s*T3)./(1+s*T4)).*(Ka./(1+s*Ta)).*((1+s*Tc)./(1+s*Tb)));


end
