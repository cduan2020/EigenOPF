function dReS_dav=ReSav(K,Y,V,lam)

dReS_dav=Sav(K,Y,V,lam)/2+conj(Sav(K,Y,V,conj(lam)))/2;

end