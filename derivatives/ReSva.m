function dReS_dva=ReSva(K,Y,V,lam)

dReS_dva=Sva(K,Y,V,lam)/2+conj(Sva(K,Y,V,conj(lam)))/2;

end