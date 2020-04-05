function dReS_dvv=ReSvv(K,Y,V,lam)

dReS_dvv=Svv(K,Y,V,lam)/2+conj(Svv(K,Y,V,conj(lam)))/2;

end