function dImS_dav=ImSav(K,Y,V,lam)

dImS_dav=Sav(K,Y,V,lam)/2j-conj(Sav(K,Y,V,conj(lam)))/2j;

end