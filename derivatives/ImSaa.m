function dImS_daa=ImSaa(K,Y,V,lam)

dImS_daa=Saa(K,Y,V,lam)/2j-conj(Saa(K,Y,V,conj(lam)))/2j;

end