function dS_daa=Saa(K,Y,V,lam)

% S=diag(K*V)*Y*V

volt=abs(V);
theta=angle(V);
nb=length(V);


dS_daa=sparsediag(K*V)*conj(Y)*sparsediag(lam)*conj(sparsediag(-V))+sparsediag(conj(Y)*conj(sparsediag(1j*V))*lam)*K*sparsediag(1j*V)...
    +diag(conj(Y*V))*K*sparsediag(lam)*sparsediag(-V)+diag(K*sparsediag(1j*V)*lam)*conj(Y)*conj(sparsediag(1j*V));

end