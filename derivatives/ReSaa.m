function dReS_daa=ReSaa(K,Y,V,lam)

dReS_daa=Saa(K,Y,V,lam)/2+conj(Saa(K,Y,V,conj(lam)))/2;

end